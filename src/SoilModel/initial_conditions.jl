"""
    get_initial_state(model::SoilModel, f::Function, t0::Real)

A function which gets the initial state for the soil model, given
a function of space (and time) `f`, as well as the initial time `t0`, 
and returns it as a pair.

This pair can then be composed with other pairs to create a single
state vector Y for the Land model.
"""
function LandHydrology.get_initial_state(
    model::SoilModel,
    f::Function,
    t0::Real,
)
    space_c, _ = make_function_space(model.domain)
    zc = Fields.coordinate_field(space_c)
    @unpack ϑ_l, θ_i, ρe_int = f.(zc, t0, Ref(model))
    return model.name =>
        Fields.FieldVector(ϑ_l = ϑ_l, θ_i = θ_i, ρe_int = ρe_int)
end


"""
    set_initial_state(model::SoilModel, f::Function, t0::Real)

A function which gets the initial state for the soil model, given
a function of space (and time) `f`, as well as the initial time `t0`.

This is useful if you just want to run the soil model alone.
"""
function LandHydrology.set_initial_state(
    model::SoilModel,
    f::Function,
    t0::Real,
)
    y = Dict(get_initial_state(model, f, t0))
    return Fields.FieldVector(; y...)
end
