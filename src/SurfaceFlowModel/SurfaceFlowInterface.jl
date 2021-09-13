module SurfaceFlowInterface
using ClimaCore: Fields
using LandHydrology.Models: AbstractModel
using LandHydrology: LandHydrologyModel
using LandHydrology.SoilInterface: SoilHydrologyModel, SoilModel, return_fluxes
import LandHydrology:
    Models, LandHydrology, set_initial_state, get_initial_state
using UnPack
export SurfaceFlowModel, make_rhs
Base.@kwdef struct SurfaceFlowModel <: AbstractModel
    name::Symbol = :surface_flow
    variables::Tuple = (:h,)
end


function Models.make_rhs(model::SurfaceFlowModel)
    rhs_surface! = make_surface_rhs!(model)
    function rhs!(dY, Y, p, t)
        rhs_surface!(dY, Y, p, t)
        return dY
    end
    return rhs!
end

function make_surface_rhs!(model)

    function rhs!(dY, Y, _, t)
        dh = dY.surface_flow.h
        h = Y.surface_flow.h
        dh[1] = t
    end
    return rhs!
end

function LandHydrology.get_initial_state(
    model::SurfaceFlowModel,
    f::Function,
    t0::Real,
)
    @unpack h = f()
    return model.name => Fields.FieldVector(h = h)
end

function LandHydrology.set_initial_state(
    model::SurfaceFlowModel,
    f::Function,
    t0::Real,
)
    y = Dict(get_initial_state(model, f, t0))
    return Fields.FieldVector(; y...)
end

end
