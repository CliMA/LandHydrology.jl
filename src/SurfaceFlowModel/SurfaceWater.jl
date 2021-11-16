module SurfaceWater
using ClimaCore: Fields, Operators
using LandHydrology.Domains: zero_field
import LandHydrology: Models, make_tendency_terms, default_initial_conditions, make_update_aux, initialize_states, AbstractModel, LandHydrologyModel, PrescribedAtmosState
using LandHydrology.SoilInterface: SoilModel
using LandHydrology.Domains: AbstractVerticalDomain, make_function_space, Column, coordinates, zero_field
export SurfaceWaterModel, get_surface_height

Base.@kwdef struct SurfaceWaterModel{FT<:AbstractFloat} <: AbstractModel
    domain::Column = Column(;zlim = (-0.05,0.05),nelements=1)
    name::Symbol = :sfc_water
end


### would like type sm to be SoilModel
function Models.make_tendency_terms(model::SurfaceWaterModel{FT},
                                    lm::LandHydrologyModel{FT, sm, SurfaceWaterModel{FT}, PrescribedAtmosState{FT, FT}},
                                    ) where {FT, sm}
    function tendency_terms!(dY, Y, Ya, t)
        infiltration = Ya.soil_infiltration
        cs = axes(Y.sfc_water.h)
        dY.sfc_water.h .= zero_field(FT, cs) .-(lm.atmos_state.precipitation - infiltration)
    end
    
    return tendency_terms!
end

function Models.default_initial_conditions(model::SurfaceWaterModel{FT}) where {FT}
    function default_ic(z::FT)
        return (; h = FT(0.0))
    end
    return Models.initialize_state(model, default_ic)
end

function Models.initialize_states(model::SurfaceWaterModel{FT} ,f::Function) where {FT}
    space_c, _ = make_function_space(model.domain)
    zc = coordinates(space_c)
    Y0 = Fields.FieldVector(; model.name => f.(zc, Ref(model)))
    Ya = Fields.FieldVector(; )
    return Y0, Ya
end

function get_surface_height(model::SurfaceWaterModel{FT}, Y) where {FT}
    return Operators.getidx(Y.sfc_water.h, Operators.Interior(),1)
end

function get_surface_height(model::AbstractModel, Y)
    return eltype(Y)(0.0)
end

        

                                                   
end
