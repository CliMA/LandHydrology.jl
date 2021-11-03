module Models
export AbstractModel, default_initial_conditions, make_rhs, make_update_aux, initialize_states
"""
    AbstractModel

An abstract type for models.

Eventually, the land model and all major subcomponents
will be of this type.
"""
abstract type AbstractModel end

"""
    default_initial_conditions(model::AbstractModel)
Construct the initial conditions for `model`.
"""
function default_initial_conditions(model::AbstractModel)
    return Fields.FieldVector(;)
end

"""
    initialize_states(model::AbstractModel)

Initialize the prognostic and auxiliary FieldVectors for `model`.
"""
function initialize_states(model::AbstractModel)
    return Fields.FieldVector(;)
end


"""
   make_rhs(model::AbstractModel)
Constructs the `rhs` (right hand side) function for the `model`, where 
the `model` defines a set of differential equations, prognostic
variables `Y`, and required parameters `Ya` for those equations, and the right hand side 
function updates in place `dY` with the continuous time derivative of each of the 
prognostic variables `Y`.

The original continuous differential equations can be a mix of partial and/or ordinary.
Given a model with e.g. a partial differential equation for `θ` and an
ordinary equation for `h`, the rhs function created by `make_rhs`
computes the *ode* rhs for `h` and the semi-discrete (ode) rhs for `θ`.
"""
function make_rhs(model::AbstractModel)
    function rhs!(dY, Y, Ya, t)
        nothing
    end
    return rhs!
end


"""
   make_update_aux(model::AbstractModel)

"""
function make_update_aux(model::AbstractModel)
    function update_aux!(Ya, t)
        nothing
    end
    return update_aux!
end




end
