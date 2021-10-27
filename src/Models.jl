module Models
export AbstractModel, AbstractLandSource, NoSource
"""
    AbstractModel

An abstract type for models.

Eventually, the land model and all major subcomponents
will be of this type.
"""
abstract type AbstractModel end

"""
    abstract type AbstractLandSource

An abstract type for source terms in the land model.

Source terms specifically refer to terms in PDEs that do not involve derivatives
of any kind. 
"""
abstract type AbstractLandSource end
struct NoSource <: AbstractLandSource end

end
