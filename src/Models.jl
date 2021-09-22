module Models
export AbstractModel, default_initial_conditions
"""
    abstract type AbstractModel

An abstract type for models.

Eventually, the land model and all major subcomponents
will be of this type.
"""
abstract type AbstractModel end

"""
    default_initial_conditions(model)
Construct the initial conditions for `model`.
"""
function default_initial_conditions(model::AbstractModel) end

end
