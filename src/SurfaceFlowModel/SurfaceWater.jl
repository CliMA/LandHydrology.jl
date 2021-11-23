module SurfaceWater
import LandHydrology: SubComponentModels, make_tendency_terms, default_initial_conditions, make_update_aux, initialize_states, AbstractModel
using LandHydrology: LandHydrologyModel, PrescribedAtmosState
using LandHydrology.SoilInterface: SoilModel
export SurfaceWaterModel, get_surface_height

"""
    SurfaceWaterModel{FT<:AbstractFloat} <: AbstractModel

A concrete type for a simple surface water model, which keeps track of the height of the surface
water as a prognostic variable. 

*****Eventually, we will want a "domain" attribute for each
submodel, which specifies the coordinates on 
which the model variables exist. This may or may not need to be a field, depending
on the model.
"""
Base.@kwdef struct SurfaceWaterModel{FT<:AbstractFloat} <: AbstractModel
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
        dY.sfc_water.h = -(lm.atmos_state.precipitation - infiltration)
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
    function default_ic(model::SurfaceWaterModel{FT})
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
    # Get coordinates of the model, from model.domain
    # Evaluate initial conditions at those coordinates
    # For now, dont need to do this.
    return f(), nothing
end

"""
    get_surface_height(model::SurfaceWaterModel{FT}, Y) where {FT}

Helper function which returns the height of the surface water.
"""
function get_surface_height(model::SurfaceWaterModel{FT}, Y) where {FT}
    return Y.sfc_water.h
end

"""
    get_surface_height(model::AbstractModel, Y) where {FT}

Helper function which returns the height of the surface water.
"""
function get_surface_height(model::AbstractModel, Y)
    return eltype(Y)(0.0)
end

        

                                                   
end
