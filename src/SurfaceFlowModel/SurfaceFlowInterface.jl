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


function Models.make_rhs(model::SurfaceFlowModel,land::LandHydrologyModel)
    rhs_surface! = make_surface_rhs!(model, land.soil.hydrology_model, land.soil)
    function rhs!(dY, Y, p, t)
        rhs_surface!(dY, Y, p, t)
        return dY
    end
    return rhs!
end

function make_surface_rhs!(model,hydrology::SoilHydrologyModel, soil::SoilModel)
    
    function rhs!(dY,Y,_,t)
        dh = dY.surface_flow.h
        h = Y.surface_flow.h
#        top_soil_bc = soil.boundary_conditions.top
#        precip = top_soil_bc.hydrology.precip(t)
#        cs = axes(Y.soil.ϑ_l)
#        fluxes = return_fluxes(Y.soil, top_soil_bc, :top, soil, cs,t)
#        infiltration = fluxes.fϑ_l
        dh[1] = t#precip - infiltration
    end
    return rhs!
end

function LandHydrology.get_initial_state(
    model::SurfaceFlowModel,
    f::Function,
    t0::Real,
)
    @unpack h = f()
    return model.name =>
        Fields.FieldVector(h = h,)
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
