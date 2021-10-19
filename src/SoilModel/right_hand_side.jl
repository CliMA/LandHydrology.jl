export make_rhs, make_update_aux, coordinates

"""
    coordinates(cs::Spaces.CenterFiniteDifferenceSpace)::Fields.Field
Returns the `z` coordinates of the space passed as an argument.
"""
coordinates(cs::Spaces.CenterFiniteDifferenceSpace)::Fields.Field =
    getproperty(Fields.coordinate_field(cs), :z)


"""
    zero_field(ft, cs::Spaces.CenterFiniteDifferenceSpace)::Fields.Field
Wrapper function returning a field on the space `cs`,
with all values = 0, of type `ft`.
"""
zero_field(ft, cs::Spaces.CenterFiniteDifferenceSpace)::Fields.Field =
    Fields.zeros(ft, cs)

"""
    make_rhs(model::SoilModel)

A function which takes a model::AbstractModel as argument, 
and returns function which computes the rhs
of a set of ordinary differential equations corresponding to
that model.

Currently, the arguments of the returned rhs function
are configured for use with OrdinaryDiffEq.jl.
For the soil model, the rhs function depends on the type of the 
components of the model (the energy and hydrology models), as well as 
whether additional sources are included.
"""
function make_rhs(model::SoilModel)
    update_aux_en! = make_update_aux(model.energy_model)
    update_aux_hydr! = make_update_aux(model.hydrology_model)
    rhs_soil! = make_rhs(model.energy_model, model.hydrology_model, model)
    function rhs!(dY, Y, Ya, t)
        update_aux_en!(Ya, t)
        update_aux_hydr!(Ya, t)
        rhs_soil!(dY, Y, Ya, t)
        return dY
    end
    return rhs!
end

"""
    make_update_aux(
        energy::PrescribedTemperatureModel,
    )

Returns a function which updates the auxiliary state vector in place by 
modifying the temperature field to the prescribed current value.
"""
function make_update_aux(energy::PrescribedTemperatureModel)
    function update_aux!(Ya, t)
        T = Ya.soil.T
        zc = Ya.zc
        T .= energy.T_profile.(zc, t)
        return Ya
    end
    return update_aux!
end

"""
    make_update_aux(
        hydrology::PrescribedHydrologyModel
    )

Returns a function which updates the auxiliary state vector in place by 
modifying the water content fields to the prescribed current value.
"""
function make_update_aux(hydrology::PrescribedHydrologyModel)
    function update_aux!(Ya, t)
        zc = Ya.zc
        @unpack ϑ_l, θ_i = Ya.soil
        ϑ_l .= hydrology.ϑ_l_profile.(zc, t)
        θ_i .= hydrology.θ_i_profile.(zc, t)
        return Ya
    end
    return update_aux!
end

"""
    make_update_aux(
        mode::AbstractSoilComponentModel
    )

Returns a function which does not update the auxilary state vector.
This is appropriate for models which do not add auxiliary state variables.
"""
function make_update_aux(model::AbstractSoilComponentModel)
    function update_aux!(Ya, t)
        nothing
    end
    return update_aux!
end


"""
    make_rhs(energy::PrescribedTemperatureModel, hydrology::PrescribedHydrologyModel, model::SoilModel)

"""
function make_rhs(
    energy::PrescribedTemperatureModel,
    hydrology::PrescribedHydrologyModel,
    model::SoilModel,
)
    function rhs!(dY, Y, Ya, t)
        nothing
    end
    return rhs!
end

"""
    make_rhs(energy::PrescribedTemperatureModel, hydrology::SoilHydrologyModel, model::SoilModel)

"""
function make_rhs(
    energy::PrescribedTemperatureModel,
    hydrology::SoilHydrologyModel{FT},
    model::SoilModel,
) where {FT}
    function rhs!(dY, Y, Ya, t)
        dϑ_l = dY.soil.ϑ_l
        dθ_i = dY.soil.θ_i
        ϑ_l = Y.soil.ϑ_l
        θ_i = Y.soil.θ_i
        T = Ya.soil.T
        zc = Ya.zc

        cspace = axes(ϑ_l)
        # construct extended state vector
        X = Fields.FieldVector(ϑ_l = ϑ_l, θ_i = θ_i, T = T) # blend of state and prescribed variables needed at boundary
        # boundary conditions and parameters
        faces = model.domain.boundary_tags
        bcs = getproperty.(Ref(model.boundary_conditions), faces)
        fluxes = (;
            zip(
                faces,
                boundary_fluxes.(
                    Ref(X),
                    bcs,
                    faces,
                    Ref(model),
                    Ref(cspace),
                    t,
                ),
            )...,
        )

        sp = model.soil_param_set
        param_set = model.earth_param_set
        hm = hydrology.hydraulic_model
        @unpack θr = hm
        @unpack ν, S_s = sp
        ν_eff = ν .- θ_i
        θ_l = volumetric_liquid_fraction.(ϑ_l, ν_eff)
        # Compute hydraulic head, conductivity
        f_i = θ_i ./ (θ_l .+ θ_i)
        viscosity_f = viscosity_factor.(Ref(hydrology.viscosity_factor), T)
        impedance_f = impedance_factor.(Ref(hydrology.impedance_factor), f_i)

        S = effective_saturation.(ν, ϑ_l, θr)
        K = hydraulic_conductivity.(Ref(hm), S, viscosity_f, impedance_f)

        ψ = pressure_head.(Ref(hm), ϑ_l, ν_eff, S_s)
        h = ψ .+ zc

        # right hand side operators
        interpc2f = Operators.InterpolateC2F()
        gradc2f_water = Operators.GradientC2F()
        divf2c_water = Operators.DivergenceF2C(
            top = Operators.SetValue(
                Geometry.Cartesian3Vector(fluxes.top.fϑ_l),
            ),
            bottom = Operators.SetValue(
                Geometry.Cartesian3Vector(fluxes.bottom.fϑ_l),
            ),
        )

        @. dϑ_l = -(divf2c_water(-interpc2f(K) * gradc2f_water(h))) #Richards equation
        dθ_i .= zero_field(FT, cspace)
        return dY
    end
    return rhs!
end

"""
    make_rhs(energy::SoilEnergyModel, hydrology::PrescribedHydrologyModel, model::SoilModel)

"""
function make_rhs(
    energy::SoilEnergyModel,
    hydrology::PrescribedHydrologyModel,
    model::SoilModel,
)
    function rhs!(dY, Y, Ya, t)
        dρe_int = dY.soil.ρe_int
        ρe_int = Y.soil.ρe_int
        ϑ_l = Ya.soil.ϑ_l
        θ_i = Ya.soil.θ_i

        # Parameters
        sp = model.soil_param_set
        param_set = model.earth_param_set
        @unpack ν, ρc_ds, κ_sat_unfrozen, κ_sat_frozen = sp

        # Compute center values of everything
        ν_eff = ν .- θ_i
        θ_l = volumetric_liquid_fraction.(ϑ_l, ν_eff)

        ρc_s = volumetric_heat_capacity.(θ_l, θ_i, ρc_ds, Ref(param_set))
        T = temperature_from_ρe_int.(ρe_int, θ_i, ρc_s, Ref(param_set))
        κ_dry = k_dry(param_set, sp)
        S_r = relative_saturation.(θ_l, θ_i, ν)
        kersten = kersten_number.(θ_i, S_r, Ref(sp))
        κ_sat =
            saturated_thermal_conductivity.(
                θ_l,
                θ_i,
                κ_sat_unfrozen,
                κ_sat_frozen,
            )
        κ = thermal_conductivity.(κ_dry, kersten, κ_sat)

        # construct extended state vector
        X = Fields.FieldVector(ϑ_l = ϑ_l, θ_i = θ_i, T = T) # blend of state and prescribed variables needed at boundary

        # boundary conditions
        cspace = axes(θ_i)
        faces = model.domain.boundary_tags
        bcs = getproperty.(Ref(model.boundary_conditions), faces)
        fluxes = (;
            zip(
                faces,
                boundary_fluxes.(
                    Ref(X),
                    bcs,
                    faces,
                    Ref(model),
                    Ref(cspace),
                    t,
                ),
            )...,
        )


        # rhs operators
        interpc2f = Operators.InterpolateC2F()
        gradc2f_heat = Operators.GradientC2F()
        divf2c_heat = Operators.DivergenceF2C(
            top = Operators.SetValue(
                Geometry.Cartesian3Vector(fluxes.top.fρe_int),
            ),
            bottom = Operators.SetValue(
                Geometry.Cartesian3Vector(fluxes.bottom.fρe_int),
            ),
        )
        @. dρe_int = -divf2c_heat(-interpc2f(κ) * gradc2f_heat(T))
        return dY
    end
    return rhs!
end

"""
    make_rhs(energy::SoilEnergyModel, hydrology::SoilHydrologyModel, model::SoilModel)

"""
function make_rhs(
    energy::SoilEnergyModel,
    hydrology::SoilHydrologyModel{FT},
    model::SoilModel,
) where {FT}
    function rhs!(dY, Y, Ya, t)
        dϑ_l = dY.soil.ϑ_l
        dθ_i = dY.soil.θ_i
        dρe_int = dY.soil.ρe_int
        ϑ_l = Y.soil.ϑ_l
        θ_i = Y.soil.θ_i
        ρe_int = Y.soil.ρe_int
        zc = Ya.zc

        # parameters
        sp = model.soil_param_set
        param_set = model.earth_param_set
        hm = hydrology.hydraulic_model
        @unpack θr = hm
        @unpack ν, S_s, ρc_ds, κ_sat_unfrozen, κ_sat_frozen = sp

        # Compute center values of everything
        ν_eff = ν .- θ_i
        θ_l = volumetric_liquid_fraction.(ϑ_l, ν_eff)
        ρc_s = volumetric_heat_capacity.(θ_l, θ_i, ρc_ds, Ref(param_set))
        T = temperature_from_ρe_int.(ρe_int, θ_i, ρc_s, Ref(param_set))
        κ_dry = k_dry(param_set, sp)
        S_r = relative_saturation.(θ_l, θ_i, ν)
        kersten = kersten_number.(θ_i, S_r, Ref(sp))
        κ_sat =
            saturated_thermal_conductivity.(
                θ_l,
                θ_i,
                κ_sat_unfrozen,
                κ_sat_frozen,
            )
        κ = thermal_conductivity.(κ_dry, kersten, κ_sat)
        ρe_int_l = volumetric_internal_energy_liq.(T, Ref(param_set))

        f_i = θ_i ./ (θ_l .+ θ_i)
        viscosity_f = viscosity_factor.(Ref(hydrology.viscosity_factor), T)
        impedance_f = impedance_factor.(Ref(hydrology.impedance_factor), f_i)
        S = effective_saturation.(ν, ϑ_l, θr)
        K = hydraulic_conductivity.(Ref(hm), S, viscosity_f, impedance_f)
        ψ = pressure_head.(Ref(hm), ϑ_l, ν_eff, S_s)
        h = ψ .+ zc

        # Construct extended state vector
        X = Fields.FieldVector(ϑ_l = ϑ_l, θ_i = θ_i, T = T) # blend of state and prescribed variables needed at boundary
        # boundary conditions
        cspace = axes(θ_i)
        faces = model.domain.boundary_tags
        bcs = getproperty.(Ref(model.boundary_conditions), faces)
        fluxes = (;
            zip(
                faces,
                boundary_fluxes.(
                    Ref(X),
                    bcs,
                    faces,
                    Ref(model),
                    Ref(cspace),
                    t,
                ),
            )...,
        )

        # rhs operators
        interpc2f = Operators.InterpolateC2F()
        gradc2f_heat = Operators.GradientC2F()
        divf2c_heat = Operators.DivergenceF2C(
            top = Operators.SetValue(
                Geometry.Cartesian3Vector(fluxes.top.fρe_int),
            ),
            bottom = Operators.SetValue(
                Geometry.Cartesian3Vector(fluxes.bottom.fρe_int),
            ),
        )

        gradc2f_water = Operators.GradientC2F()
        divf2c_water = Operators.DivergenceF2C(
            top = Operators.SetValue(
                Geometry.Cartesian3Vector(fluxes.top.fϑ_l),
            ),
            bottom = Operators.SetValue(
                Geometry.Cartesian3Vector(fluxes.bottom.fϑ_l),
            ),
        )

        @. dϑ_l = -divf2c_water(-interpc2f(K) * gradc2f_water(h)) #Richards equation
        dθ_i .= zero_field(FT, cspace)

        @. dρe_int =
            -divf2c_heat(
                -interpc2f(κ) * gradc2f_heat(T) -
                interpc2f(ρe_int_l * K) * gradc2f_water(h),
            )
        return dY
    end
    return rhs!
end
