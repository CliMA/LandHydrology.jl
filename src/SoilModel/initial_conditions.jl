export set_initial_state

"""
    set_initial_state(model::SoilModel, f::Function, t0::Real)

A function which sets the initial state for the soil model, given
a function of space (and time) `f`, as well as the initial time `t0`.
"""
function set_initial_state(model::SoilModel, f::Function, t0::Real)
    space_c, _ = make_function_space(model.domain)
    zc = coordinates(space_c)
    @unpack ϑ_l, θ_i, ρe_int = f.(zc, t0, Ref(model))
    Y = Fields.FieldVector(ϑ_l = ϑ_l, θ_i = θ_i, ρe_int = ρe_int)
    return Y
end
