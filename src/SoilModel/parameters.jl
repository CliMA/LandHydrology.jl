export SoilParams

abstract type ParameterStructure{FT <: AbstractFloat} end

"""
    struct SoilParams{FT} <: ParameterStructure{FT}
"""
struct SoilParams{FT} <: ParameterStructure{FT}
    ν::FT
    S_s::FT
    ν_ss_gravel::FT
    ν_ss_om::FT
    ν_ss_quartz::FT
    ρc_ds::FT
    κ_solid::FT
    ρp::FT
    κ_sat_unfrozen::FT
    κ_sat_frozen::FT
    a::FT
    b::FT
    κ_dry_parameter::FT
end
