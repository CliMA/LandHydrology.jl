module LandHydrology
export SoilParams,
    SoilWaterParams,
    SoilHeatParams

abstract type ParameterStructure{FT <: AbstractFloat} end

struct SoilParams{FT} <: ParameterStructure{FT}
    ν::FT
    vgn::FT
    vgα::FT
    vgm::FT
    ksat::FT
    θr::FT
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

struct SoilWaterParams{FT} <: ParameterStructure{FT}
    ν::FT
    vgn::FT
    vgα::FT
    vgm::FT
    ksat::FT
    θr::FT
    S_s::FT
end


struct SoilHeatParams{FT} <: ParameterStructure{FT}
    ν::FT
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


include(joinpath("SoilModel", "SoilInterface.jl"))
include(joinpath("SoilModel", "SoilWaterParameterizations.jl"))
include(joinpath("SoilModel", "SoilHeatParameterizations.jl"))

end # module
