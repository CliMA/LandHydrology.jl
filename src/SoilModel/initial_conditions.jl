"""
    Models.initialize_auxiliary(model::SoilModel, t0::Real, zc)

Initializes and returns the auxiliary state

The initial time is required as prescribed variables can be 
functions of space and time.
"""
function Models.initialize_auxiliary(model::SoilModel, t0::Real, zc)
    init_aux_coordinates = (z) -> (; zc = z)
    init_aux_energy = aux_vars(model.energy_model)
    init_aux_hydrology = aux_vars(model.hydrology_model)
    function init_aux_soil(t, z)
        return (; init_aux_coordinates(z)..., init_aux_energy(t, z)..., init_aux_hydrology(t, z)...)
    end
    return init_aux_soil.(t0, zc)
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
     Models.initialize_prognostic(model::SoilModel, f::Function, zc)

This function evaluates the initial condition
function `f` at the domain points `zc` and creates the initial state vector `Y`.
"""
function Models.initialize_prognostic(model::SoilModel, f::Function, zc)
    return f.(zc, Ref(model))
end

"""
     Models.initialize_states(model::SoilModel, f::Function, t0::Real)

This function returns the initial prognostic and auxiliary states for the model,
given an initial condition function `f` for the prognostic variables, and
an initial time `t0`. 

In the future, this could be split into two functions, one for aux and one for prognostic variables,
if we can create the space twice or create the space elsewhere and pass in.
"""
function Models.initialize_states(model::SoilModel, f::Function, t0::Real)
    space_c, _ = make_function_space(model.domain)
    zc = coordinates(space_c)
    Y0 = Fields.FieldVector(; model.name => Models.initialize_prognostic(model, f, zc)
    Ya0 = Fields.FieldVector(; model.name => Models.initialize_auxiliary(model, t0, zc)
    return Y0, Ya0
end
