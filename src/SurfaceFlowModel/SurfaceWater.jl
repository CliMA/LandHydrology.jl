module SurfaceWater
using ClimaCore: Fields
import LandHydrology: Models, make_tendency_terms, default_initial_conditions, make_update_aux, initialize_states, AbstractModel, LandHydrologyModel, PrescribedAtmosState
using LandHydrology.SoilInterface: SoilModel, SoilHydrologyModel, AbstractSoilComponentModel
using LandHydrology.Domains: AbstractVerticalDomain, make_function_space, Column, coordinates, zero_field
export SurfaceWaterModel

Base.@kwdef struct SurfaceWaterModel{FT<:AbstractFloat} <: AbstractModel
    domain::Column = Column(;zlim = (-0.05,0.05),nelements=1)
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
        c_space= axes(Y.sfc_water.h)
        Ya.sfc_water.infiltration .= zero_field(FT, c_space)
    end
    return update_aux!
end



#function Models.make_update_aux(model::SurfaceWaterModel)
#end


function Models.make_tendency_terms(model::SurfaceWaterModel{FT},
                                    lm::LandHydrologyModel{FT, sm, SurfaceWaterModel{FT}, PrescribedAtmosState{FT, FT}},
                                    ) where {FT, sm}
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

function Models.initialize_states(model::SurfaceWaterModel{FT} ,f::Function) where {FT}
    space_c, _ = make_function_space(model.domain)
    zc = coordinates(space_c)
    Y0 = Fields.FieldVector(; model.name => f.(zc, Ref(model)))
    f_aux = initialize_aux_function(model)
    Ya = Fields.FieldVector(; model.name => f_aux.(zc))
    return Y0, Ya
end

function initialize_aux_function(model::SurfaceWaterModel{FT}) where {FT}
    function f_aux(z::FT)
        return (; infiltration = FT(0.0))
    end
    return f_aux
end

        

                                                   
end
