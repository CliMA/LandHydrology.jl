#model.sources will be a tuple eventually?
export PhaseChange, NoSource

function θ_star(T, θ_l, θ_i, ν, hm, Tfreeze, g, LH_f0, ρice,ρliq)
    θ_m = min(ρice* θ_i / ρliq + θ_l, ν)
    @unpack θr = hm
    ψ0 = matric_potential(hm, θ_m)## isnt this wrong???
    ψT = (LH_f0 / g / Tfreeze) * (T - Tfreeze)
    if T < Tfreeze
        θ_s = θr + (ν - θr) * inverse_matric_potential(hm, ψ0+ψT)
    else
        θ_s = θ_l
    end
    return θ_s
end

struct NoSource <: AbstractLandSource{FT} end
    

"""
    PhaseChange <: AbstractLandSource
The function which computes the freeze/thaw source term for Richard's equation.
"""
Base.@kwdef struct PhaseChange{FT} <: AbstractLandSource{FT}
    "Typical resolution in the vertical"
    Δz::FT = FT(NaN)
    "O(1) factor multiplying the thermal time"
    γ_LTE::FT = FT(1)
end


function make_source(source::AbstractLandSource, energy::AbstractSoilComponentModel,
                      hydrology::AbstractSoilComponentModel, model::SoilModel)
    function source_soil!(dY,Y,p,t)
        nothing
    end
    return source_soil!
end

                          
function make_source(source::PhaseChange{FT},energy::SoilEnergyModel,
                      hydrology::SoilHydrologyModel{FT},
                      model::SoilModel,
                      ) where{FT}
        
    function source_soil!(dY,Y,p,t)
        ϑ_l = Y.ϑ_l
        θ_i = Y.θ_i
        ρe_int = Y.ρe_int
        
        #parameters
        
        param_set = model.earth_param_set
        sp = model.soil_param_set
        @unpack ν, ρc_ds  = sp
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
        τLTE = ρc_s .* Δz^2 ./ κ .*γ_LTE
        
        freeze_thaw =
            FT(1) ./ τLTE .* 
            _ρliq .*
            (θ_l .- θ_s) #*
        @. dϑ_l += -freeze_thaw ./ _ρliq
        @. dθ_i += freeze_thaw ./ _ρice
        return dY
    end
    return source_soil!
end
