module SurfaceFlow
import SoilInterface: compute_infiltration
import .Models: make_rhs, make_tendency_terms, default_initial_conditions, make_update_aux, initialize_states

Base.@kwdef struct SurfaceFlowModel{FT<:AbstractFloat} <: AbstractModel
    name::Symbol = :sfc_flow
end


#This does not need this function
#function Models.make_update_aux(model::SurfaceFlowModel)
#end


function Models.make_tendency_terms(model::SurfaceFlowModel{FT},
                                    lm::LandHydrologyModel{SoilModel,
                                                           SurfaceFlowModel{FT}
                                                           }
                                    ) where {FT, dm}
    function tendency_terms!(dY, Y, Ya, t)
        infiltration = SoilInterface.compute_infiltration(lm.soil.hydrology_model, Y, Ya)
        @. dY.sfc_flow.h = lm.atmos_state.precipitation(t) - infiltration
    end
    
    return tendency_terms!
end


# so you can run this model by itself
function Models.make_rhs(model::SurfaceFlowModel{FT})
    update_aux! = Models.make_update_aux(model)
    tendency_terms = Models.make_tendency_terms!(model)
    function rhs!(dY,Y,Ya,t)
        update_aux!(Ya,t)
        tendency_terms!(dY,Y,Ya,t)
    end
    return rhs!
end

function Models.default_initial_conditions(model::SurfaceFlowModel{FT}) where {FT}
    Y = Fields.FieldVector(; sfc_flow => [FT(0.0)])
end

function Models.initialize_states(model::SurfaceFlowModel{FT} ,f::Function, _) where {FT}
    return Fields.FieldVector(; sfc_flow => FT(f())), Fields.FieldVector(;)
end

end
