export make_rhs

function make_rhs(model::SoilModel)
    rhs_soil! = make_rhs!(model.energy_model, model.hydrology_model, model)
    function rhs!( dY, Y, p, t)
        rhs_soil!(dY,Y,p,t) 
        return dY
    end
    return rhs!
end

function make_rhs!(energy::PrescribedTemperatureModel, hydrology::PrescribedHydrologyModel, model)
    function rhs!(dY, Y, _, t)
        (Y_hydro, Y_energy) = Y.x
        @unpack ϑ_l, θ_i = Y_hydro
        @unpack ρe_int = Y_energy
        (dY_hydro, dY_energy) = dY.x
        dρe_int = dY_energy.ρe_int
        dϑ_l = dY_hydro.ϑ_l
        dθ_i = dY_hydro.θ_i

        #update the state Y with prescribed values at the current time
        ϑ_l = hydrology.ϑ_l_profile.(zc, t)
        θ_i = hydrology.θ_i_profile.(zc, t)
        θ_l = ϑ_l # eventually have a conversion to liquid water content
        T = energy.T_profile.(zc, t)
        ρc_s = volumetric_heat_capacity.(θ_l, θ_i, model.soil_param_set.ρc_ds, Ref(model.earth_param_set))
        ρe_int = volumetric_internal_energy.(θ_i, ρc_s, T, Ref(model.earth_param_set))
        
        # RHS is zero

        cs = axes(θ_i)
        zc = Fields.coordinate_field(cs)
        dϑ_l = Fields.zeros(typeof(θ_i),cs)
        dθ_i = Fields.zeros(typeof(θ_i),cs)
        dρe_int = Fields.zeros(typeof(ρe_int),cs)
        return dY
    end
    return rhs!
end

# Richards equation with temperature dependent hydraulic conductivity (optional)
function make_rhs!(energy::PrescribedTemperatureModel, hydrology::SoilHydrologyModel, model)
    function rhs!(dY, Y, _, t)
        (dY_hydro, dY_energy) = dY.x
        dϑ_l = dY_hydro.ϑ_l
        dθ_i = dY_hydro.θ_i
        dρe_int = dY_energy.ρe_int

        (Y_hydro, Y_energy) = Y.x
        @unpack ϑ_l, θ_i = Y_hydro
        @unpack ρe_int = Y_energy

        # boundary conditions and parameters
        bc = model.boundary_conditions
        top_water_flux, btm_water_flux = compute_vertical_flux(bc.top.hydrology), compute_vertical_flux(bc.bottom.hydrology)
        sp = model.soil_param_set
        param_set= model.earth_param_set
        @unpack ν, vgn, vgα, vgm, ksat, θr, S_s, ρc_ds = sp
        
        # Evaluate T at the current time, update ρe_int
        cs = axes(θ_i)
        zc = Fields.coordinate_field(cs)
        θ_l = ϑ_l
        T = energy.T_profile.(zc, t)
        ρc_s = volumetric_heat_capacity.(θ_l, θ_i, model.soil_param_set.ρc_ds, Ref(model.earth_param_set))
        ρe_int = volumetric_internal_energy.(θ_i, ρc_s, T, Ref(model.earth_param_set))
        dρe_int =Fields.zeros(eltype(ρe_int),cs)
        
        # Compute hydraulic head, conductivity
        S = effective_saturation.(θ_l; ν = ν, θr = θr)
        K = hydraulic_conductivity.(S; vgm = vgm, ksat = ksat)
        ψ = pressure_head.(
            S;
            vgn = vgn,
            vgα = vgα,
            vgm = vgm,
            ν = ν,
            θr = θr,
            S_s = S_s,
        )

        h = ψ .+ zc

        # right hand side operators
        interpc2f = Operators.InterpolateC2F()
        gradc2f_water = Operators.GradientC2F()
        divf2c_water = Operators.DivergenceF2C(
            top = Operators.SetValue(top_water_flux),
            bottom = Operators.SetValue(btm_water_flux),
        )
        
        @. dϑ_l = -(divf2c_water(-interpc2f(K) * gradc2f_water(h))) #Richards equation
        dθ_i = Fields.zeros(eltype(θ_i),cs)
        return dY
    end
    return rhs!
end

# Simple heat equation in soil. Water profile is prescribed. No heat flux due to diffusion of water.
function make_rhs!(energy::SoilEnergyModel, hydrology::PrescribedHydrologyModel, model::SoilModel)
    function rhs!(dY, Y, _, t)
        (dY_hydro, dY_energy) = dY.x
        dρe_int = dY_energy.ρe_int
        dϑ_l = dY_hydro.ϑ_l
        dθ_i = dY_hydro.θ_i
        (Y_hydro, Y_energy) = Y.x
        @unpack ϑ_l, θ_i = Y_hydro
        @unpack ρe_int = Y_energy

        # boundary conditions and parameters
        bc = model.boundary_conditions
        top_heat_flux, btm_heat_flux = compute_vertical_flux(bc.top.energy), compute_vertical_flux(bc.bottom.energy)
        sp = model.soil_param_set
        param_set= model.earth_param_set
        @unpack ν, ρc_ds, κ_sat_unfrozen, κ_sat_frozen = sp

        # update water content based on prescribed profiles, set RHS to zero.
        cs = axes(ρe_int)
        zc = Fields.coordinate_field(cs)
        ϑ_l = hydrology.ϑ_l_profile.(zc, t)
        θ_i = hydrology.θ_i_profile.(zc, t)
        dϑ_l = Fields.zeros(typeof(θ_i),cs)
        dθ_i = Fields.zeros(typeof(θ_i),cs)

        # Compute center values of everything
        ρc_s = volumetric_heat_capacity.(θ_l, θ_i, ρc_ds, Ref(param_set))
        T = temperature_from_ρe_int.(ρe_int, θ_i, ρc_s, Ref(param_set))
        κ_dry = k_dry(param_set, sp)
        S_r = relative_saturation.(θ_l, θ_i, ν)
        kersten = kersten_number.(θ_i, S_r, Ref(sp))
        κ_sat =
            saturated_thermal_conductivity.(θ_l, θ_i, κ_sat_unfrozen, κ_sat_frozen)
        κ = thermal_conductivity.(κ_dry, kersten, κ_sat)

        # rhs operators
        interpc2f = Operators.InterpolateC2F()
        gradc2f_heat = Operators.GradientC2F()
        divf2c_heat = Operators.DivergenceF2C(
            top = Operators.SetValue(top_heat_flux),
            bottom = Operators.SetValue(btm_heat_flux),
        )
        @. dρe_int =
            -divf2c_heat(
                -interpc2f(κ) * gradc2f_heat(T))
        return dY
    end
    return rhs!
end

# Full model without phase changes

function make_rhs!(energy::SoilEnergyModel, hydrology::SoilHydrologyModel, model::SoilModel)
    function rhs!(dY, Y, _, t)
        (dY_hydro, dY_energy) = dY.x
        dρe_int = dY_energy.ρe_int
        dϑ_l = dY_hydro.ϑ_l
        dθ_i = dY_hydro.θ_i
        (Y_hydro, Y_energy) = Y.x
        @unpack ϑ_l, θ_i = Y_hydro
        @unpack ρe_int = Y_energy

        # boundary conditions and parameters
        bc = model.boundary_conditions
        top_water_flux, btm_water_flux = compute_vertical_flux(bc.top.hydrology), compute_vertical_flux(bc.bottom.hydrology)
        top_heat_flux, btm_heat_flux = compute_vertical_flux(bc.top.energy), compute_vertical_flux(bc.bottom.energy)

        sp = model.soil_param_set
        param_set= model.earth_param_set
        @unpack ν, vgn, vgα, vgm, ksat, θr, S_s, ρc_ds, κ_sat_unfrozen, κ_sat_frozen = sp
        
        # Compute center values of everything
        θ_l = ϑ_l
        ρc_s = volumetric_heat_capacity.(θ_l, θ_i, ρc_ds, Ref(param_set))
        T = temperature_from_ρe_int.(ρe_int, θ_i, ρc_s, Ref(param_set))
        κ_dry = k_dry(param_set, sp)
        S_r = relative_saturation.(θ_l, θ_i, ν)
        kersten = kersten_number.(θ_i, S_r, Ref(sp))
        κ_sat =
            saturated_thermal_conductivity.(θ_l, θ_i, κ_sat_unfrozen, κ_sat_frozen)
        κ = thermal_conductivity.(κ_dry, kersten, κ_sat)
        ρe_int_l = volumetric_internal_energy_liq.(T, Ref(param_set))
        
        cs = axes(θ_i)
        zc = Fields.coordinate_field(cs)
        
        S = effective_saturation.(θ_l; ν = ν, θr = θr)
        K = hydraulic_conductivity.(S; vgm = vgm, ksat = ksat)
        ψ = pressure_head.(
            S;
            vgn = vgn,
            vgα = vgα,
            vgm = vgm,
            ν = ν,
            θr = θr,
            S_s = S_s,
        )
        h = ψ .+ zc

        # rhs operators
        interpc2f = Operators.InterpolateC2F()
        gradc2f_heat = Operators.GradientC2F()
        divf2c_heat = Operators.DivergenceF2C(
            top = Operators.SetValue(top_heat_flux),
            bottom = Operators.SetValue(btm_heat_flux),
        )
        
        gradc2f_water = Operators.GradientC2F()
        divf2c_water = Operators.DivergenceF2C(
            top = Operators.SetValue(top_water_flux),
            bottom = Operators.SetValue(btm_water_flux),
        )
        
        @. dϑ_l = -divf2c_water(-interpc2f(K) * gradc2f_water(h)) #Richards equation
        dθ_i = Fields.zeros(eltype(θ_i), cs) # zero unless we add phase change        

        @. dρe_int =
            -divf2c_heat(
                -interpc2f(κ) * gradc2f_heat(T) -
                interpc2f(ρe_int_l * K) * gradc2f_water(h),
            )
        return dY
    end
    return rhs!
end
