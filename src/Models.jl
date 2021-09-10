module Models
export AbstractModel, make_rhs
"""
    abstract type AbstractModel

An abstract type for models.

Eventually, the land model and all major subcomponents
will be of this type.
"""
abstract type AbstractModel end

"""
    function make_rhs(model::AbstractModel)
"""
function make_rhs(model::AbstractModel)
    function rhs!(dY,Y,p,t)
        nothing
    end
    return rhs!
end

end
