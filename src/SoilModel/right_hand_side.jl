# We may still be able to preserve not using LandModel (make_rhs(soil_model))

function SubComponentModels.make_update_aux(model::SoilModel,lm::LandHydrologyModel)
    update_aux_en! = SubComponentModels.make_update_aux(model.energy_model, lm)
    update_aux_hydr! = SubComponentModels.make_update_aux(model.hydrology_model, lm)
    function update_aux!(Ya, Y, t)
        update_aux_en!(Ya,Y, t)
        update_aux_hydr!(Ya,Y, t)
    end
    return update_aux!
end




"""
            SubComponentModels.make_update_aux(
                energy::PrescribedTemperatureModel,
            )
    
        Returns a function which updates the auxiliary state vector in place by 
        modifying the temperature field to the prescribed current value.
        """
function SubComponentModels.make_update_aux(energy::PrescribedTemperatureModel,_)
    function update_aux!(Ya, Y,t)
        T = Ya.soil.T
        zc = Ya.soil.zc
        T .= energy.T_profile.(zc, t)
        return Ya
    end
    return update_aux!
end

"""
        SubComponentModels.make_update_aux(
        hydrology::PrescribedHydrologyModel
    )

Returns a function which updates the auxiliary state vector in place by 
modifying the water content fields to the prescribed current value.
"""
function SubComponentModels.make_update_aux(hydrology::PrescribedHydrologyModel,_)
    function update_aux!(Ya,Y, t)
        zc = Ya.soil.zc
        @unpack ϑ_l, θ_i = Ya.soil
        ϑ_l .= hydrology.ϑ_l_profile.(zc, t)
        θ_i .= hydrology.θ_i_profile.(zc, t)
        return Ya
    end
    return update_aux!
end


"""
    SubComponentModels.make_update_aux(
        mode::AbstractSoilComponentModel
    )

Returns a function which does not update the auxilary state vector.
This is appropriate for models which do not add auxiliary state variables.
"""
function SubComponentModels.make_update_aux(model::AbstractSoilComponentModel, _)
    function update_aux!(Ya,Y, t)
        nothing
    end
    return update_aux!
end


"""
    SubComponentModels.make_tendency_terms(model::SoilModel{FT, dm, PrescribedTemperatureModel, PrescribedHydrologyModel},) where {FT,dm}

"""
function SubComponentModels.make_tendency_terms(model::SoilModel{FT, dm, PrescribedTemperatureModel, PrescribedHydrologyModel},lm::LandHydrologyModel) where {FT, dm}
    function tendency_terms!(dY, Y, Ya, t)
        nothing
    end
    return tendency_terms!
end

"""
    SubComponentModels.make_tendency_terms(model::SoilModel{FT, dm, PrescribedTemperatureModel, SoilHydrologyModel{FT}},) where {FT, dm}

"""
function SubComponentModels.make_tendency_terms(
    model::SoilModel{FT, dm, PrescribedTemperatureModel, SoilHydrologyModel{FT}},lm::LandHydrologyModel) where {FT, dm}
    function tendency_terms!(dY, Y, Ya, t)
        hydrology = model.hydrology_model
        dϑ_l = dY.soil.ϑ_l
        dθ_i = dY.soil.θ_i
        ϑ_l = Y.soil.ϑ_l
        θ_i = Y.soil.θ_i
        T = Ya.soil.T
        zc = Ya.soil.zc

        cspace = axes(ϑ_l)
        # boundary conditions and parameters
        faces = model.domain.boundary_tags
        bcs = getproperty.(Ref(model.boundary_conditions), faces)
        fluxes = (;
            zip(
                faces,
                boundary_fluxes.(
                    Ref(Y),
                    Ref(Ya),
                    bcs,
                    faces,
                    Ref(lm),
                    Ref(cspace),
                    t,
                ),
            )...,
        )
        sp = model.soil_param_set
        param_set = lm.earth_param_set
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
                Geometry.WVector(fluxes.top.fϑ_l),
            ),
            bottom = Operators.SetValue(
                Geometry.WVector(fluxes.bottom.fϑ_l),
            ),
        )

        @. dϑ_l = -(divf2c_water(-interpc2f(K) * gradc2f_water(h))) #Richards equation
        dθ_i .= zero_field(FT, cspace)
    end
    return tendency_terms!
end

"""
    SubComponentModels.make_tendency_terms(model::SoilModel{FT, dm, SoilEnergyModel, PrescribedHydrologyModel},) where {FT, dm}

"""
function SubComponentModels.make_tendency_terms(model::SoilModel{FT, dm, SoilEnergyModel, PrescribedHydrologyModel},lm::LandHydrologyModel) where {FT, dm}

    function tendency_terms!(dY, Y, Ya, t)
        energy = model.energy_model
        dρe_int = dY.soil.ρe_int
        ρe_int = Y.soil.ρe_int
        ϑ_l = Ya.soil.ϑ_l
        θ_i = Ya.soil.θ_i

        # Parameters
        sp = model.soil_param_set
        param_set = lm.earth_param_set
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

        # boundary conditions
        cspace = axes(θ_i)
        faces = model.domain.boundary_tags
        bcs = getproperty.(Ref(model.boundary_conditions), faces)
        fluxes = (;
            zip(
                faces,
                boundary_fluxes.(
                    Ref(Y),
                    Ref(Ya),
                    bcs,
                    faces,
                    Ref(lm),
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
                Geometry.WVector(fluxes.top.fρe_int),
            ),
            bottom = Operators.SetValue(
                Geometry.WVector(fluxes.bottom.fρe_int),
            ),
        )
        @. dρe_int = -divf2c_heat(-interpc2f(κ) * gradc2f_heat(T))
    end
    return tendency_terms!
end

"""
    SubComponentModels.make_tendency_terms(model::SoilModel{FT, dm, SoilEnergyModel, SoilHydrologyModel{FT}}, ) where {FT, dm}

"""
function SubComponentModels.make_tendency_terms(model::SoilModel{FT, dm, SoilEnergyModel, SoilHydrologyModel{FT}},lm::LandHydrologyModel) where {FT, dm}
    function tendency_terms!(dY, Y, Ya, t)
        energy = model.energy_model
        hydrology = model.hydrology_model
        
        dϑ_l = dY.soil.ϑ_l
        dθ_i = dY.soil.θ_i
        dρe_int = dY.soil.ρe_int
        ϑ_l = Y.soil.ϑ_l
        θ_i = Y.soil.θ_i
        ρe_int = Y.soil.ρe_int
        zc = Ya.soil.zc

        # parameters
        sp = model.soil_param_set
        param_set = lm.earth_param_set
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

        # boundary conditions
        cspace = axes(θ_i)
        faces = model.domain.boundary_tags
        bcs = getproperty.(Ref(model.boundary_conditions), faces)
        fluxes = (;
            zip(
                faces,
                boundary_fluxes.(
                    Ref(Y),
                    Ref(Ya),
                    bcs,
                    faces,
                    Ref(lm),
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
                Geometry.WVector(fluxes.top.fρe_int),
            ),
            bottom = Operators.SetValue(
                Geometry.WVector(fluxes.bottom.fρe_int),
            ),
        )

        gradc2f_water = Operators.GradientC2F()
        divf2c_water = Operators.DivergenceF2C(
            top = Operators.SetValue(
                Geometry.WVector(fluxes.top.fϑ_l),
            ),
            bottom = Operators.SetValue(
                Geometry.WVector(fluxes.bottom.fϑ_l),
            ),
        )

        @. dϑ_l = -divf2c_water(-interpc2f(K) * gradc2f_water(h)) #Richards equation
        dθ_i .= zero_field(FT, cspace)

        @. dρe_int =
            -divf2c_heat(
                -interpc2f(κ) * gradc2f_heat(T) -
                interpc2f(ρe_int_l * K) * gradc2f_water(h),
            )
    end
    return tendency_terms!
end
