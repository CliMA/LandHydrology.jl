export initialize_prognostic, initialize_auxiliary, initialize_states

"""
    initialize_auxiliary(model::SoilModel, t0::Real, zc)

Initializes and returns the auxiliary state FieldVector.

This function creates an initialization function for each subcomponent of the model,
by calling `aux_vars` for that model, and then creates the final FieldVector of aux
states by evaluating those functions for each model on the coordinates of our domain
and at the initial time. The initial time is required as prescribed variables can be 
functions of space and time.
"""
function initialize_auxiliary(model::SoilModel, t0::Real, zc) # eventually to be called with LandModel type
    init_aux_soil = aux_vars(model) # eventually will be model.soil
    return Fields.FieldVector(; :zc => zc, model.name => init_aux_soil.(t0, zc))
end

"""
    aux_vars(model::SoilModel)

Returns a function of space and time which can be used to compute the initial
state of the auxiliary variables of the soil model at the initial time t0.

This implicitly defines which auxiliary variables are included for the model `model`.
"""
function aux_vars(model::SoilModel)
    init_aux_energy = aux_vars(model.energy_model)
    init_aux_hydrology = aux_vars(model.hydrology_model)
    function init_aux_soil(t, z)
        return (; init_aux_energy(t, z)..., init_aux_hydrology(t, z)...)
    end
    return init_aux_soil
end

"""
    aux_vars(m::PrescribedTemperatureModel)

Returns a function of space and time which can be used to compute the
initial state of the auxiliary variables of the PrescribedTemperatureModel
a time t0.

This implicitly defines which auxiliary variables are included for the model `m`.
"""
function aux_vars(m::PrescribedTemperatureModel)
    h = m.T_profile
    return (t, z) -> (; T = h(z, t))
end

"""
    aux_vars(m::PrescribedHydrologyModel)

Returns a function of space and time which can be used to compute the
initial state of the auxiliary variables of the PrescribedHydrologyModel
a time t0.

This implicitly defines which auxiliary variables are included for the model `m`.
"""
function aux_vars(m::PrescribedHydrologyModel)
    f = m.ϑ_l_profile
    g = m.θ_i_profile
    return (t, z) -> (; ϑ_l = f(z, t), θ_i = g(z, t))

end

"""
    aux_vars(m::AbstractSoilComponentModel)

Returns a function of space and time which can be used to compute the
initial state of the auxiliary variables of a default AbstractSoilComponentModel
a time t0.

This implicitly defines which auxiliary variables are included for the model `m`.
"""
function aux_vars(m::AbstractSoilComponentModel)
    return (t, z) -> (;)
end

"""
     initialize_prognostic(model::SoilModel, f::Function, zc)

This function evaluates the initial condition
function `f` at the domain points `zc` and creates the initial state vector `Y`.
"""
function initialize_prognostic(model::SoilModel, f::Function, zc)
    soil_ic_tuple = f.(zc, Ref(model))
    Y = Fields.FieldVector(; model.name => soil_ic_tuple)
    return Y
end

"""
     initialize_states(model::SoilModel, f::Function, t0::Real)

This function returns the initial prognostic and auxiliary states for the model,
given an initial condition function `f` for the prognostic variables, and
an initial time `t0`. 

In the future, this could be split into two functions, one for aux and one for prognostic variables,
if we can create the space twice or create the space elsewhere and pass in.
"""
function initialize_states(model::SoilModel, f::Function, t0::Real)
    space_c, _ = make_function_space(model.domain)
    zc = coordinates(space_c)
    Y0 = initialize_prognostic(model, f, zc)
    Ya0 = initialize_auxiliary(model, t0, zc)
    return Y0, Ya0
end
