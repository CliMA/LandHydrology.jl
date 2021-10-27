#model.sources will be a tuple eventually?
export PhaseChange, add_source

function θ_star(T, θ_l, θ_i, ν, hm, Tfreeze, g, LH_f0, ρice,ρliq)
    θ_m = min(ρice* θ_i / ρliq + θ_l, ν)
    @unpack θr = hm
    S_m = (θ_m - θr)/(ν-θr)# νeff??
    ψ0 = matric_potential(hm, S_m)
    ψT = (LH_f0 / g / Tfreeze) * (T - Tfreeze)
    if T < Tfreeze
        θ_s = θr + (ν - θr) * inverse_matric_potential(hm, ψ0+ψT)
    else
        θ_s = θ_m
    end
    return θ_s
end
    

"""
    PhaseChange <: AbstractLandSource
The function which computes the freeze/thaw source term for Richard's equation.
"""
Base.@kwdef struct PhaseChange{FT <: AbstractFloat} <: AbstractLandSource
    "Typical resolution in the vertical"
    Δz::FT = FT(NaN)
    "O(1) factor multiplying the thermal time"
    γ_LTE::FT = FT(1)
end

                          
function make_source(model::SoilModel)
    source! = add_source(model.sources, model)
    function source_soil!(dY,Y,p,t)
        source!(dY,Y,p,t)
    end
    return source_soil!
end

function add_source(source::AbstractLandSource, model::SoilModel)
    function source!(dY,Y,p,t)
        nothing
    end
    return source!
end


function add_source(source::PhaseChange{FT}, model::SoilModel) where {FT}

    function source!(dY,Y,p,t)
        dϑ_l = dY.soil.ϑ_l
        dθ_i = dY.soil.θ_i
        ϑ_l = Y.soil.ϑ_l
        θ_i = Y.soil.θ_i
        ρe_int = Y.soil.ρe_int
        
        #parameters
        hydrology = model.hydrology_model
        param_set = model.earth_param_set
        sp = model.soil_param_set
        @unpack ν, ρc_ds,κ_sat_unfrozen,κ_sat_frozen  = sp
        hm = hydrology.hydraulic_model
        @unpack θr = hm
        _ρliq = FT(ρ_cloud_liq(param_set))
        _ρice = FT(ρ_cloud_ice(param_set))
        _Tfreeze = FT(T_freeze(param_set))
        _LH_f0 = FT(LH_f0(param_set))
        _g = FT(grav(param_set))
        
        
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
        θ_s = θ_star.(T, θ_l, θ_i, ν, Ref(hm), _Tfreeze, _g, _LH_f0, _ρice, _ρliq)
        
        Δz = source.Δz
        γ_LTE = source.γ_LTE
        τLTE = (ρc_s .* Δz^2 ./ κ) .*γ_LTE
        freeze_thaw =
            FT(1) ./ τLTE .* 
            _ρliq .*
            (θ_l .- θ_s) #*
        @. dϑ_l += - freeze_thaw / _ρliq
        @. dθ_i += freeze_thaw / _ρice
    end
    return source!
end
