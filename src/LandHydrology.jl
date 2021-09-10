module LandHydrology
export LandHydrologyModel, set_initial_state
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
            rhs_function!(dY,Y,p,t)
        end
        return dY
    end
    return rhs!
end


function set_initial_state(model::LandHydrologyModel) end
include(joinpath("SoilModel", "SoilInterface.jl"))




end # module
