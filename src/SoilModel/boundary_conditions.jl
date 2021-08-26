export FluxBC, SoilBC, NoBC
abstract type AbstractBC end

struct NoBC <: AbstractBC
end
struct FluxBC{FT<:AbstractFloat} <: AbstractBC
    top_flux::FT
    btm_flux::FT
end

Base.@kwdef struct SoilBC{ebc<: AbstractBC, hbc<:AbstractBC}
    energy_bc::ebc = NoBC()
    hydrology_bc::hbc = NoBC()
end
