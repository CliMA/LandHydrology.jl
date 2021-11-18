module SubComponentModels
using LandHydrology: AbstractLandModel
export AbstractModel,
    default_initial_conditions,
    make_tendency_terms,
    make_update_aux,
    initialize_states,
    NotIncluded
"""
    AbstractModel

An abstract type for subcomponent models of the LandHydrologyModel
"""
abstract type AbstractModel end

"""
    NotIncluded <: AbstractModel

The concrete type for subcomponents when they are not included in the model (do not contribute prognostic or auxiliary variables).
"""
struct NotIncluded <: AbstractModel end


"""
    default_initial_conditions(model::AbstractModel, lm::AbstractLandModel)

Construct the initial conditions for subcomponent `model`, given any additional information 
(e.g. the shared parameter set) stored in the `lm` LandHydrologyModel.
"""
function default_initial_conditions(model::AbstractModel,lm::AbstractLandModel)
    if typeof(model) != NotIncluded
        println("Warning: no default initial conditions for ", model.name, " model.")
    end
              
    return nothing, nothing
end

"""
    initialize_states(model::AbstractModel, lm::AbstractLandModel)

Initialize the prognostic and auxiliary variables for the subcomponent `model`, given
any additional information required from the `lm` LandHydrologyModel. 

Currently, this should return a `Field` type object.
"""
function initialize_states(model::AbstractModel, lm::AbstractLandModel,land_f::NamedTuple)
    return nothing, nothing
end

"""
   make_tendency_terms(model::AbstractModel, lm::AbstractLandModel)

Constructs a function which, when evaluated, gives the value of a particular
(set of) tendency terms  appearing in the right hand side of a(n system of)
ordinary differential equation(s) for a given (set of) variable(s). 

The arguments of this function must be: `dY`, the right hand side of the prognostic
variables `Y`, required auxiliary parameters `Ya`, and the time `t`. This function should
update `dY` in place with the tendency under consideration.
"""
function make_tendency_terms(model::AbstractModel,  lm::AbstractLandModel)
    function tendency_terms!(dY, Y, Ya, t)
        nothing
    end
    return tendency_terms!
end


"""
   make_update_aux(model::AbstractModel)

Constructs a function which updates the auxilary parameters of the subcomponent model
`model`, given any additional information stored in the Land model `lm`. 

Auxiliary functions can be a function of space, time, or prognostic variables.
The returned function should update `Ya` in place.

"""
function make_update_aux(model::AbstractModel, lm::AbstractLandModel)
    function update_aux!(Ya, Y, t)
        nothing
    end
    return update_aux!
end




end
