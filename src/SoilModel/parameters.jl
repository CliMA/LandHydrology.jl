export SoilParams

abstract type ParameterStructure{FT <: AbstractFloat} end

"""
    SoilParams{FT} <: ParameterStructure{FT}
Defaults correspond to loam soil. 
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SoilParams{FT} <: ParameterStructure{FT}
    "Porosity"
    ν::FT = 0.43
    "Specific storage"
    S_s::FT = 1e-3
    "Volumetric fraction of soil solids in gravel"
    ν_ss_gravel::FT = 0.0
    "Volumetric fraction of soil solids in organic matter"
    ν_ss_om::FT = 0.0
    "Volumetric fraction of soil solids in quartz/sand"
    ν_ss_quartz::FT = 0.41
    "Volumetrcy heat capacicity of dry soil"
    ρc_ds::FT = 2700
    "Thermal conductivity of soil solids"
    κ_solid::FT = 3.97
    "Particle density"
    ρp::FT = 2700.0
    "Thermal conductivity of saturated unfrozen soil"
    κ_sat_unfrozen::FT = 1.72
    "Thermal conductivity of saturated frozen soil"
    κ_sat_frozen::FT = 3.13
    "Parameter needed in Balland & Arp model for Kersten number"
    a::FT = 0.24
    "Parameter needed in Balland & Arp model for Kersten number"
    b::FT = 18.1
    "Parameter needed in Balland & Arp model for κ_dry"
    κ_dry_parameter::FT = 0.053
end
