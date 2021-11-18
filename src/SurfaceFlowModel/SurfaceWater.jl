module SurfaceWater
using ClimaCore: Fields, Operators
using LandHydrology.Domains: zero_field
import LandHydrology: SubComponentModels, make_tendency_terms, default_initial_conditions, make_update_aux, initialize_states, AbstractModel
using LandHydrology: LandHydrologyModel, PrescribedAtmosState
using LandHydrology.SoilInterface: SoilModel
using LandHydrology.Domains: AbstractVerticalDomain, make_function_space, Column, coordinates, zero_field
export SurfaceWaterModel, get_surface_height

"""
    SurfaceWaterModel{FT<:AbstractFloat} <: AbstractModel

A concrete type for a simple surface water model, which keeps track of the height of the surface
water as a prognostic variable. 

*****The `domain` fields should be reconsidered for non-column models. For column models,
this shouldn't need to be a Field.
"""
Base.@kwdef struct SurfaceWaterModel{FT<:AbstractFloat} <: AbstractModel
    domain::Column = Column(;zlim = (-0.05,0.05),nelements=1)
    name::Symbol = :sfc_water
end


"""
    SubComponentModels.make_tendency_terms(model::SurfaceWaterModel{FT},
                                lm::LandHydrologyModel{FT, sm, SurfaceWaterModel{FT},},
                                ) where {FT, sm}

Constructs the function which evaluates the tendency of the surface water height.

****This requires an atmosphere model, so we should also be using that in the parametric
type signature of the LandHydrologyModel.
"""
function SubComponentModels.make_tendency_terms(model::SurfaceWaterModel{FT},
                                    lm::LandHydrologyModel{FT, sm, SurfaceWaterModel{FT},},
                                    ) where {FT, sm}
    function tendency_terms!(dY, Y, Ya, t)
        infiltration = Ya.soil_infiltration
        cs = axes(Y.sfc_water.h)
        dY.sfc_water.h .= zero_field(FT, cs) .-(lm.atmos_state.precipitation - infiltration)
    end
    
    return tendency_terms!
end

"""
    SubComponentModels.default_initial_conditions(model::SurfaceWaterModel{FT},
                                                  lm::LandHydrologyModel{FT, sm, SurfaceWaterModel{FT}, }
                                                  ) where {FT, sm}

Sets up the default initial condition function for the `SurfaceWaterModel`.
"""
function SubComponentModels.default_initial_conditions(model::SurfaceWaterModel{FT},
                                                       lm::LandHydrologyModel{FT, sm, SurfaceWaterModel{FT}, }
                                                       ) where {FT, sm}
    function default_ic(z::FT, model::SurfaceWaterModel{FT})
        return (; h = FT(0.0))
    end
    return SubComponentModels.initialize_states(model, lm, (; :sfc_water => default_ic))
end

"""
    SubComponentModels.initialize_states(model::SurfaceWaterModel{FT},
                                         lm::LandHydrologyModel{FT, sm, SurfaceWaterModel{FT},},
                                         f::NamedTuple) where {FT, sm}

Given an initial condition function `land_f.sfc_water`, computes the initial prognostic state of the `SurfaceWaterModel`.

"""
function SubComponentModels.initialize_states(model::SurfaceWaterModel{FT},
                                              lm::LandHydrologyModel{FT, sm, SurfaceWaterModel{FT},},
                                              land_f::NamedTuple) where {FT, sm}
    f = getproperty(land_f, model.name)
    space_c, _ = make_function_space(model.domain)
    zc = coordinates(space_c)
    return f.(zc, Ref(model)), nothing
end

"""
    get_surface_height(model::SurfaceWaterModel{FT}, Y) where {FT}

Helper function which returns the height of the surface water.
"""
function get_surface_height(model::SurfaceWaterModel{FT}, Y) where {FT}
    return Operators.getidx(Y.sfc_water.h, Operators.Interior(),1)
end

"""
    get_surface_height(model::AbstractModel, Y) where {FT}

Helper function which returns the height of the surface water.
"""
function get_surface_height(model::AbstractModel, Y)
    return eltype(Y)(0.0)
end

        

                                                   
end
