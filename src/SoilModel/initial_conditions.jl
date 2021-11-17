export initialize_aux_function
"""
     create_aux_ic_function(model::SoilModel)

Creates and returns the auxilary state initialization function.

The auxiliary profiles (z,t) are determined via the fields
of the model; to change them from the defaults, they must
be adjusted at the soil component model level (model.energy_model,
model.hydrology_model).
"""
function initialize_aux_function(model::SoilModel)
    init_aux_coordinates = (z) -> (; zc = z)
    init_aux_energy = aux_vars(model.energy_model)
    init_aux_hydrology = aux_vars(model.hydrology_model)
    function f_aux(z)
        return (; init_aux_coordinates(z)..., init_aux_energy(z)..., init_aux_hydrology(z)...)
    end
    return f_aux
end

"""
    aux_vars(m::PrescribedTemperatureModel)

Returns a function of space and time which can be used to compute the
initial state of the auxiliary variables of the PrescribedTemperatureModel
a time t0.

This implicitly defines which auxiliary variables are included for the model `m`.
"""
function aux_vars(m::PrescribedTemperatureModel)
    return (z) -> (; T = typeof(z)(0.0))
end

"""
    aux_vars(m::PrescribedHydrologyModel)

Returns a function of space and time which can be used to compute the
initial state of the auxiliary variables of the PrescribedHydrologyModel
a time t0.

This implicitly defines which auxiliary variables are included for the model `m`.
"""
function aux_vars(m::PrescribedHydrologyModel)
    return (z) -> (; ϑ_l = typeof(z)(0.0), θ_i = typeof(z)(0.0))

end

"""
    aux_vars(m::AbstractSoilComponentModel)

Returns a function of space and time which can be used to compute the
initial state of the auxiliary variables of a default AbstractSoilComponentModel
a time t0.

This implicitly defines which auxiliary variables are included for the model `m`.
"""
function aux_vars(m::AbstractSoilComponentModel)
    return (z) -> (;)
end

"""
     SubComponentModels.initialize_states(model::SoilModel, f::Function, t0::Real)

This function returns the initial prognostic and auxiliary states for the model,
given an initial condition function `f` for the prognostic variables, and
an initial time `t0`. 

In the future, this could be split into two functions, one for aux and one for prognostic variables,
if we can create the space twice or create the space elsewhere and pass in.
"""
function SubComponentModels.initialize_states(model::SoilModel, lm::LandHydrologyModel, f::Function)
    space_c, _ = make_function_space(model.domain)
    zc = coordinates(space_c)
    f_aux = initialize_aux_function(model)
    return f.(zc, Ref(model), Ref(lm)), f_aux.(zc)
end
