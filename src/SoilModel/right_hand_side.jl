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
        @unpack ??_l, ??_i = Ya.soil
        ??_l .= hydrology.??_l_profile.(zc, t)
        ??_i .= hydrology.??_i_profile.(zc, t)
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
        d??_l = dY.soil.??_l
        d??_i = dY.soil.??_i
        ??_l = Y.soil.??_l
        ??_i = Y.soil.??_i
        T = Ya.soil.T
        zc = Ya.zc

        cspace = axes(??_l)
        # construct extended state vector
        X = Fields.FieldVector(??_l = ??_l, ??_i = ??_i, T = T) # blend of state and prescribed variables needed at boundary
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
        @unpack ??r = hm
        @unpack ??, S_s = sp
        ??_eff = ?? .- ??_i
        ??_l = volumetric_liquid_fraction.(??_l, ??_eff)
        # Compute hydraulic head, conductivity
        f_i = ??_i ./ (??_l .+ ??_i)
        viscosity_f = viscosity_factor.(Ref(hydrology.viscosity_factor), T)
        impedance_f = impedance_factor.(Ref(hydrology.impedance_factor), f_i)

        S = effective_saturation.(??, ??_l, ??r)
        K = hydraulic_conductivity.(Ref(hm), S, viscosity_f, impedance_f)

        ?? = pressure_head.(Ref(hm), ??_l, ??_eff, S_s)
        h = ?? .+ zc

        # right hand side operators
        interpc2f = Operators.InterpolateC2F()
        gradc2f_water = Operators.GradientC2F()
        divf2c_water = Operators.DivergenceF2C(
            top = Operators.SetValue(
                Geometry.Cartesian3Vector(fluxes.top.f??_l),
            ),
            bottom = Operators.SetValue(
                Geometry.Cartesian3Vector(fluxes.bottom.f??_l),
            ),
        )

        @. d??_l = -(divf2c_water(-interpc2f(K) * gradc2f_water(h))) #Richards equation
        d??_i .= zero_field(FT, cspace)
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
        d??e_int = dY.soil.??e_int
        ??e_int = Y.soil.??e_int
        ??_l = Ya.soil.??_l
        ??_i = Ya.soil.??_i

        # Parameters
        sp = model.soil_param_set
        param_set = model.earth_param_set
        @unpack ??, ??c_ds, ??_sat_unfrozen, ??_sat_frozen = sp

        # Compute center values of everything
        ??_eff = ?? .- ??_i
        ??_l = volumetric_liquid_fraction.(??_l, ??_eff)

        ??c_s = volumetric_heat_capacity.(??_l, ??_i, ??c_ds, Ref(param_set))
        T = temperature_from_??e_int.(??e_int, ??_i, ??c_s, Ref(param_set))
        ??_dry = k_dry(param_set, sp)
        S_r = relative_saturation.(??_l, ??_i, ??)
        kersten = kersten_number.(??_i, S_r, Ref(sp))
        ??_sat =
            saturated_thermal_conductivity.(
                ??_l,
                ??_i,
                ??_sat_unfrozen,
                ??_sat_frozen,
            )
        ?? = thermal_conductivity.(??_dry, kersten, ??_sat)

        # construct extended state vector
        X = Fields.FieldVector(??_l = ??_l, ??_i = ??_i, T = T) # blend of state and prescribed variables needed at boundary

        # boundary conditions
        cspace = axes(??_i)
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
                Geometry.Cartesian3Vector(fluxes.top.f??e_int),
            ),
            bottom = Operators.SetValue(
                Geometry.Cartesian3Vector(fluxes.bottom.f??e_int),
            ),
        )
        @. d??e_int = -divf2c_heat(-interpc2f(??) * gradc2f_heat(T))
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
        d??_l = dY.soil.??_l
        d??_i = dY.soil.??_i
        d??e_int = dY.soil.??e_int
        ??_l = Y.soil.??_l
        ??_i = Y.soil.??_i
        ??e_int = Y.soil.??e_int
        zc = Ya.zc

        # parameters
        sp = model.soil_param_set
        param_set = model.earth_param_set
        hm = hydrology.hydraulic_model
        @unpack ??r = hm
        @unpack ??, S_s, ??c_ds, ??_sat_unfrozen, ??_sat_frozen = sp

        # Compute center values of everything
        ??_eff = ?? .- ??_i
        ??_l = volumetric_liquid_fraction.(??_l, ??_eff)
        ??c_s = volumetric_heat_capacity.(??_l, ??_i, ??c_ds, Ref(param_set))
        T = temperature_from_??e_int.(??e_int, ??_i, ??c_s, Ref(param_set))
        ??_dry = k_dry(param_set, sp)
        S_r = relative_saturation.(??_l, ??_i, ??)
        kersten = kersten_number.(??_i, S_r, Ref(sp))
        ??_sat =
            saturated_thermal_conductivity.(
                ??_l,
                ??_i,
                ??_sat_unfrozen,
                ??_sat_frozen,
            )
        ?? = thermal_conductivity.(??_dry, kersten, ??_sat)
        ??e_int_l = volumetric_internal_energy_liq.(T, Ref(param_set))

        f_i = ??_i ./ (??_l .+ ??_i)
        viscosity_f = viscosity_factor.(Ref(hydrology.viscosity_factor), T)
        impedance_f = impedance_factor.(Ref(hydrology.impedance_factor), f_i)
        S = effective_saturation.(??, ??_l, ??r)
        K = hydraulic_conductivity.(Ref(hm), S, viscosity_f, impedance_f)
        ?? = pressure_head.(Ref(hm), ??_l, ??_eff, S_s)
        h = ?? .+ zc

        # Construct extended state vector
        X = Fields.FieldVector(??_l = ??_l, ??_i = ??_i, T = T) # blend of state and prescribed variables needed at boundary
        # boundary conditions
        cspace = axes(??_i)
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
                Geometry.Cartesian3Vector(fluxes.top.f??e_int),
            ),
            bottom = Operators.SetValue(
                Geometry.Cartesian3Vector(fluxes.bottom.f??e_int),
            ),
        )

        gradc2f_water = Operators.GradientC2F()
        divf2c_water = Operators.DivergenceF2C(
            top = Operators.SetValue(
                Geometry.Cartesian3Vector(fluxes.top.f??_l),
            ),
            bottom = Operators.SetValue(
                Geometry.Cartesian3Vector(fluxes.bottom.f??_l),
            ),
        )

        @. d??_l = -divf2c_water(-interpc2f(K) * gradc2f_water(h)) #Richards equation
        d??_i .= zero_field(FT, cspace)

        @. d??e_int =
            -divf2c_heat(
                -interpc2f(??) * gradc2f_heat(T) -
                interpc2f(??e_int_l * K) * gradc2f_water(h),
            )
        return dY
    end
    return rhs!
end
