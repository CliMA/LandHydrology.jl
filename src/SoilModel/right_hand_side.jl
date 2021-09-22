export make_rhs

"""
    coordinates(cs::Spaces.CenterFiniteDifferenceSpace)::Fields.Field

Returns the coordinates of the space passed as an argument. For columns,
this returns the `z` values.
"""
coordinates(cs::Spaces.CenterFiniteDifferenceSpace)::Fields.Field =
    Fields.local_geometry_field(cs).coordinates


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
    rhs_soil! = make_rhs!(model.energy_model, model.hydrology_model, model)
    function rhs!(dY, Y, p, t)
        rhs_soil!(dY, Y, p, t)
        return dY
    end
    return rhs!
end

"""
    make_rhs!(energy::PrescribedTemperatureModel, hydrology::PrescribedHydrologyModel, model::SoilModel)

"""
function make_rhs!(
    energy::PrescribedTemperatureModel,
    hydrology::PrescribedHydrologyModel,
    model::SoilModel,
)
    function rhs!(dY, Y, _, t)
        dϑ_l = dY.ϑ_l
        dθ_i = dY.θ_i
        dρe_int = dY.ρe_int
        ϑ_l = Y.ϑ_l
        θ_i = Y.θ_i
        ρe_int = Y.ρe_int

        #update the state Y with prescribed values at the current time
        ϑ_l = hydrology.ϑ_l_profile.(zc, t)
        θ_i = hydrology.θ_i_profile.(zc, t)
        θ_l = ϑ_l # eventually have a conversion to liquid water content
        T = energy.T_profile.(zc, t)
        ρc_s =
            volumetric_heat_capacity.(
                θ_l,
                θ_i,
                model.soil_param_set.ρc_ds,
                Ref(model.earth_param_set),
            )
        ρe_int =
            volumetric_internal_energy.(
                θ_i,
                ρc_s,
                T,
                Ref(model.earth_param_set),
            )

        # RHS is zero
        cspace = axes(ϑ_l)
        zc = coordinates(cspace)
        FT = eltype(zc)

        dϑ_l = zero_field(FT, cspace)
        dθ_i = zero_field(FT, cspace)
        dρe_int = zero_field(FT, cspace)
        return dY
    end
    return rhs!
end

"""
    make_rhs!(energy::PrescribedTemperatureModel, hydrology::SoilHydrologyModel, model::SoilModel)

"""
function make_rhs!(
    energy::PrescribedTemperatureModel,
    hydrology::SoilHydrologyModel,
    model::SoilModel,
)
    function rhs!(dY, Y, _, t)
        dϑ_l = dY.ϑ_l
        dθ_i = dY.θ_i
        dρe_int = dY.ρe_int
        ϑ_l = Y.ϑ_l
        θ_i = Y.θ_i
        ρe_int = Y.ρe_int

        cspace = axes(ϑ_l)
        zc = coordinates(cspace)
        FT = eltype(zc)

        # Evaluate T at the current time, update ρe_int
        θ_l = ϑ_l
        T = energy.T_profile.(zc, t)
        ρc_s =
            volumetric_heat_capacity.(
                θ_l,
                θ_i,
                model.soil_param_set.ρc_ds,
                Ref(model.earth_param_set),
            )
        ρe_int =
            volumetric_internal_energy.(
                θ_i,
                ρc_s,
                T,
                Ref(model.earth_param_set),
            )

        # boundary conditions and parameters
        faces = model.domain.x3boundary
        bcs = getproperty.(Ref(model.boundary_conditions), faces)
        fluxes = (;
            zip(
                faces,
                boundary_fluxes.(
                    Ref(Y),
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
        @unpack ν, S_s, ρc_ds = sp
        hm = hydrology.hydraulic_model
        @unpack θr = hm

        # Compute hydraulic head, conductivity
        S = effective_saturation.(θ_l; ν = ν, θr = θr)
        K = hydraulic_conductivity.(Ref(hm), S)
        ψ = pressure_head.(Ref(hm), S; ν = ν, S_s = S_s)

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
        dθ_i = zero_field(FT, cspace)
        dρe_int = zero_field(FT, cspace)

        return dY
    end
    return rhs!
end

"""
    make_rhs!(energy::SoilEnergyModel, hydrology::PrescribedHydrologyModel, model::SoilModel)

"""
function make_rhs!(
    energy::SoilEnergyModel,
    hydrology::PrescribedHydrologyModel,
    model::SoilModel,
)
    function rhs!(dY, Y, _, t)
        dϑ_l = dY.ϑ_l
        dθ_i = dY.θ_i
        dρe_int = dY.ρe_int
        ϑ_l = Y.ϑ_l
        θ_i = Y.θ_i
        ρe_int = Y.ρe_int

        cspace = axes(ϑ_l)
        zc = coordinates(cspace)
        FT = eltype(zc)

        # update water content based on prescribed profiles, set RHS to zero.
        ϑ_l = hydrology.ϑ_l_profile.(zc, t)
        θ_i = hydrology.θ_i_profile.(zc, t)
        dϑ_l = zero_field(FT, cspace)
        dθ_i = zero_field(FT, cspace)

        # boundary conditions and parameters
        faces = model.domain.x3boundary
        bcs = getproperty.(Ref(model.boundary_conditions), faces)
        fluxes = (;
            zip(
                faces,
                boundary_fluxes.(
                    Ref(Y),
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
        @unpack ν, ρc_ds, κ_sat_unfrozen, κ_sat_frozen = sp

        # Compute center values of everything
        θ_l = ϑ_l
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
    make_rhs!(energy::SoilEnergyModel, hydrology::SoilHydrologyModel, model::SoilModel)

"""
function make_rhs!(
    energy::SoilEnergyModel,
    hydrology::SoilHydrologyModel,
    model::SoilModel,
)
    function rhs!(dY, Y, _, t)
        dϑ_l = dY.ϑ_l
        dθ_i = dY.θ_i
        dρe_int = dY.ρe_int
        ϑ_l = Y.ϑ_l
        θ_i = Y.θ_i
        ρe_int = Y.ρe_int

        cspace = axes(ϑ_l)
        zc = coordinates(cspace)
        FT = eltype(zc)

        # boundary conditions and parameters
        faces = model.domain.x3boundary
        bcs = getproperty.(Ref(model.boundary_conditions), faces)
        fluxes = (;
            zip(
                faces,
                boundary_fluxes.(
                    Ref(Y),
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
        @unpack ν, S_s, ρc_ds, κ_sat_unfrozen, κ_sat_frozen = sp

        # Compute center values of everything
        θ_l = ϑ_l
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

        S = effective_saturation.(θ_l; ν = ν, θr = θr)
        K = hydraulic_conductivity.(Ref(hm), S)
        ψ = pressure_head.(Ref(hm), S; ν = ν, S_s = S_s)
        h = ψ .+ zc

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
        dθ_i = zero_field(FT, cspace)

        @. dρe_int =
            -divf2c_heat(
                -interpc2f(κ) * gradc2f_heat(T) -
                interpc2f(ρe_int_l * K) * gradc2f_water(h),
            )
        return dY
    end
    return rhs!
end
