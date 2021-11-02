# Deal with default_initial_conditions
module LandHydrology
using CLIMAParameters
using DocStringExtensions

struct EarthParameterSet <: AbstractEarthParameterSet end

include("Domains/Domains.jl")
include("Models.jl")
using .Models: AbstractModel
import .Models: make_rhs, default_initial_conditions, make_update_aux, initialize_states


"""
    LandHydrologyModel{FT<: AbstractFloat, sm <: AbstractModel} <: AbstractModel

The structure which defines all of the physics and setup of the land 
hydrology model, including domains, parameters, prognostic variables, differential equations, and included interactions between subcomponents.
# Fields
$(DocStringExtensions.FIELDS)

"""
struct LandHydrologyModel{FT<: AbstractFloat, sm <: AbstractModel} <: AbstractModel
    "The soil model"
    soil::sm
end

"""
    function Models.make_rhs(model::LandHydrologyModel)

Makes the rhs!(dY, Y, Ya t) ODE function for the entire land model,
where dY is the ODE function evaluated at the current prognostic state `Y`,
time `t`, and parameters `Ya`. 

If your ODEs result from the semi-discrete form of a PDE,
the rhs function must give the rhs of the discretized system. For example, 
the equation for the single prognostic variable T

``
∂T/∂t = ∂²T/∂z²
``

 will become a set of ODEs

``
∂ T_i/∂t = f_i(T_i, T_i+1, T_i-1, Δz,...)
``

where the form of `f`, number of equations, etc, depends on your chosen 
spatial discretization. The `make_rhs` function returns a function that 
updates in place `dY = dT⃗ = (dT₁, dT₂,....)` with the values of `f⃗`.


If you do not define a 
method of `make_rhs!` for your specific model type, the default is not
to update the prognostic state Y.
"""
function Models.make_rhs(model::LandHydrologyModel)
    rhs_functions = map(propertynames(model)) do sc_name
        sc = getproperty(model, sc_name)
        Models.make_rhs(sc)
    end

    function rhs!(dY, Y, Ya, t)
        map(rhs_functions) do rhs_function
            rhs_function!(dY,Y, Ya, t)
        end
    end
    return rhs!
end

function make_update_aux(model::LandHydrologyModel)
    aux_functions = map(propertynames(model)) do sc_name
        sc = getproperty(model, sc_name)
        Models.make_update_aux(sc)
    end

    function update_aux!(Ya, t)
        map(aux_functions) do aux_function
            aux_function!(Ya, t)
        end
    end
    return update_aux!
end

function initialize_states(model::LandHydrologyModel)
end



include(joinpath("SoilModel", "SoilInterface.jl"))
include("Simulations/Simulations.jl")

end # module
