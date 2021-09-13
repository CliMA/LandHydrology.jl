module LandHydrology
export LandHydrologyModel, set_initial_state, get_initial_state, NotIncluded
using ClimaCore: Fields
include("Domains/Domains.jl")
include("Models.jl")
using .Models: AbstractModel
import .Models: make_rhs

struct NotIncluded <: AbstractModel end

Base.@kwdef struct LandHydrologyModel <: AbstractModel
    soil::AbstractModel = NotIncluded()
    surface_flow::AbstractModel = NotIncluded()
end

"""
    function Models.make_rhs(model::LandHydrologyModel)

Makes the rhs!(dY, Y, p t) ODE function for the entire land model.

Note that the default `make_rhs!(model::AbstractModel)` does not
alter any of the state variables. This is appropriate for subcomponents
 of type `NotIncluded`, for example, and may be appropriate for subcomponents
that are prescribed as well. The upshot is: if you do not define a 
method of `make_rhs!` for your specific model type, it will not adjust Y.
"""
function Models.make_rhs(model::LandHydrologyModel)
    subcomponents = propertynames(model)
    rhs_functions = []
    for sc_name in subcomponents
        sc = getproperty(model, sc_name)
        sc_rhs! = Models.make_rhs(sc)
        push!(rhs_functions, sc_rhs!)
    end

    function rhs!(dY, Y, p, t)
        for rhs_function! in rhs_functions
            rhs_function!(dY, Y, p, t)
        end
        return dY
    end
    return rhs!
end

"""
    function set_initial_state(model::LandHydrologyModel)

Compose initial state out of initial conditions for each subcomponent.

If your model subcomponent type is `NotIncluded`, no variables are added to the RHS.
"""
function set_initial_state(model::LandHydrologyModel, f::NamedTuple, t0::Real)
    subcomponents = propertynames(model)
    Ys = Dict()
    for sc_name in subcomponents
        sc = getproperty(model, sc_name)
        if typeof(sc) != NotIncluded
            f_sc = getproperty(f, sc_name)
            Y_sc = get_initial_state(sc, f_sc, t0)
            push!(Ys, Y_sc)
        end

    end
    Y = Fields.FieldVector(; Ys...)
    return Y

end

function get_initial_state(model::AbstractModel, f, t0::Real) end # not sure if we need this

include(joinpath("SoilModel", "SoilInterface.jl"))
include(joinpath("SurfaceFlowModel", "SurfaceFlowInterface.jl"))




end # module
