module SurfaceWater
import LandHydrology: Models, make_tendency_terms, default_initial_conditions, make_update_aux, initialize_states, AbstractModel, LandHydrologyModel
using LandHydrology.SoilInterface: SoilModel
using LandHydrology.Domains: AbstractVerticalDomain, make_function_space, Column
export SurfaceWaterModel

Base.@kwdef struct SurfaceWaterModel{FT<:AbstractFloat} <: AbstractModel
    domain::Column = Column(;zlim = [-0.05,-0.05],nelements=1)
    name::Symbol = :sfc_water
end


"""
    Models.make_update_aux(
        mode::AbstractSoilComponentModel
    )

Returns a function which does not update the auxilary state vector.
This is appropriate for models which do not add auxiliary state variables.
"""
function Models.make_update_aux(model::SurfaceWaterModel{FT},lm::LandHydrologyModel) where {FT}
    function update_aux!(Ya,Y, t)
        Ya.sfc_water.infiltration = FT(0.0) # could use Y.sfc_water.h in the one case, set h = 0 in the other case.
    end
    return update_aux!
end


# no need, can fall back on default
#function Models.make_update_aux(model::SurfaceWaterModel)
#end


function Models.make_tendency_terms(model::SurfaceWaterModel{FT},
                                    lm::LandHydrologyModel{SoilModel, SurfaceWaterModel{FT}, PrescribedAtmosState{FT}},
                                    ) where {FT}
    function tendency_terms!(dY, Y, Ya, t)
        infiltration = Ya.sfc_water.infiltration
        @. dY.sfc_water.h = lm.atmos_state.precipitation - infiltration
    end
    
    return tendency_terms!
end
#these wont work
function Models.default_initial_conditions(model::SurfaceWaterModel{FT}) where {FT}
    nothing
end

function Models.initialize_states(model::SurfaceWaterModel{FT} ,f::Function, _) where {FT}
    space_c, _ = make_function_space(model.domain)
    zc = coordinates(space_c)
    Y0 = Fields.FieldVector(; model.name => f.(zc, Ref(model)))
    Ya = Fields.FieldVector(; model.name => f_aux.(t0, 
    return Fields.FieldVector(; sfc_water => FT(f())), Fields.FieldVector(; sfc_water => 0.0)
end

end
