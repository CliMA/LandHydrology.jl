export make_rhs

function make_rhs(model::NoSoilModel)
    function rhs!( dY, Y, p, t)
    end
    return rhs!
end



function make_rhs(model::SoilModel)
    
    function rhs!( dY, Y, p, t)
        (Y_hydro, Y_energy) = Y.x
        (dY_hydro, dY_energy) = dY.x
        @unpack ϑ_l, θ_i = Y_hydro
        @unpack ρe_int = Y_energy
        
        dϑ_l = dY_hydro.ϑ_l
        dθ_i = dY_hydro.θ_i
        dρe_int = dY_energy.ρe_int


        heat_flux = model.boundary_conditions.energy_bc
        top_heat_flux, btm_heat_flux = heat_flux.top_flux, heat_flux.btm_flux
        water_flux = model.boundary_conditions.hydrology_bc
        top_water_flux, btm_water_flux = water_flux.top_flux, water_flux.btm_flux

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
        ψ = matric_potential.(S; vgn = vgn, vgα = vgα, vgm = vgm)
        h = ψ .+ zc
        
        interpc2f = Operators.InterpolateC2F()
        gradc2f_heat = Operators.GradientC2F()
        gradf2c_heat = Operators.GradientF2C(
            top = Operators.SetValue(top_heat_flux),
            bottom = Operators.SetValue(btm_heat_flux),
        )
        
        gradc2f_water = Operators.GradientC2F()
        gradf2c_water = Operators.GradientF2C(
            top = Operators.SetValue(top_water_flux),
            bottom = Operators.SetValue(btm_water_flux),
        )
        
        @. dϑ_l = -gradf2c_water(-interpc2f(K) * gradc2f_water(h)) #Richards equation
        @. dρe_int =
            -gradf2c_heat(
                -interpc2f(κ) * gradc2f_heat(T) -
                interpc2f(ρe_int_l * K) * gradc2f_water(h),
            )
        dθ_i = Fields.zeros(eltype(θ_i), cs)
        
        return dY
    end
    return rhs!
end

function make_rhs(model::RichardsEquation)
    
    function rhs!( dY, Y, p, t)
        (Y_hydro,) = Y.x
        (dY_hydro,) = dY.x
        @unpack ϑ_l, θ_i = Y_hydro
        
        dϑ_l = dY_hydro.ϑ_l
        dθ_i = dY_hydro.θ_i

        water_flux = model.boundary_conditions.hydrology_bc
        top_water_flux, btm_water_flux = water_flux.top_flux, water_flux.btm_flux


        sp = model.soil_param_set
        param_set= model.earth_param_set
        @unpack ν, vgn, vgα, vgm, ksat, θr, S_s = sp
        # Compute center values of everything
        θ_l = ϑ_l
        cs = axes(θ_l)
        zc = Fields.coordinate_field(cs)
        T = model.T_profile(zc, t)
      
        S = effective_saturation.(θ_l; ν = ν, θr = θr)
        K = hydraulic_conductivity.(S; vgm = vgm, ksat = ksat)
        ψ = pressure_head.(S; vgn = vgn, vgα = vgα, vgm = vgm, ν = ν, θr = θr, S_s = S_s)
        h = ψ .+ zc
        
        interpc2f = Operators.InterpolateC2F()
        gradc2f_water = Operators.GradientC2F()
        gradf2c_water = Operators.GradientF2C(
            top = Operators.SetValue(top_water_flux),
            bottom = Operators.SetValue(btm_water_flux),
        )
        
        @. dϑ_l = -gradf2c_water(-interpc2f(K) * gradc2f_water(h)) #Richards equation
   
        dθ_i = Fields.zeros(eltype(θ_i), cs)
        
        return dY
    end
    return rhs!
end


function make_rhs(model::SoilHeatEquation)
    
    function rhs!( dY, Y, p, t)
        (dY_energy,) = dY.x
        dρe_int = dY_energy.ρe_int
        heat_flux = model.boundary_conditions.energy_bc
        top_heat_flux, btm_heat_flux = heat_flux.top_flux, heat_flux.btm_flux
        
        sp = model.soil_param_set
        param_set= model.earth_param_set
        @unpack ν, ρc_ds, κ_sat_unfrozen, κ_sat_frozen = sp
        
        #the energy RHS needs water variables
        (Y_energy, ) = Y.x
        @unpack ρe_int = Y_energy

        sp = model.soil_param_set
        param_set= model.earth_param_set
        @unpack ν, ρc_ds, κ_sat_unfrozen, κ_sat_frozen = sp
        
        cs = axes(ρe_int)
        zc = Fields.coordinate_field(cs)
        θ_l = model.θ_l_profile(zc, t)
        θ_i = model.θ_i_profile(zc, t)
        
        # Compute center values of everything
        ρc_s = volumetric_heat_capacity.(θ_l, θ_i, ρc_ds, Ref(param_set))
        T = temperature_from_ρe_int.(ρe_int, θ_i, ρc_s, Ref(param_set))
        κ_dry = k_dry(param_set, sp)
        S_r = relative_saturation.(θ_l, θ_i, ν)
        kersten = kersten_number.(θ_i, S_r, Ref(sp))
        κ_sat =
            saturated_thermal_conductivity.(θ_l, θ_i, κ_sat_unfrozen, κ_sat_frozen)
        κ = thermal_conductivity.(κ_dry, kersten, κ_sat)
        
        interpc2f = Operators.InterpolateC2F()
        gradc2f_heat = Operators.GradientC2F()
        gradf2c_heat = Operators.GradientF2C(
            top = Operators.SetValue(top_heat_flux),
            bottom = Operators.SetValue(btm_heat_flux),
        )
        @. dρe_int =
            -gradf2c_heat(
                -interpc2f(κ) * gradc2f_heat(T))
        return dY
        
        return dY
    end
    return rhs!
end
