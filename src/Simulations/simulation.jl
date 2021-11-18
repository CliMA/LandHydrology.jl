"""
    struct Simulation <: AbstractSimulation

A simulation wraps an abstract LandHydrology `model` containing 
equation specifications and an instance of an `integrator` used for
time integration of the discretized model PDE.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Simulation{ML <: AbstractLandModel} <: AbstractSimulation
    "A LandHydrology model"
    model::ML
    "a DiffEqBase.jl integrator used for time"
    integrator::DiffEqBase.DEIntegrator
    "user defined callback operations"
    callbacks::Union{
        DiffEqBase.CallbackSet,
        DiffEqBase.DiscreteCallback,
        Nothing,
    }
end

"""
    Simulation(model::AbstractModel, method::AbstractTimestepper, 
               dt, tspan, init_state, callbacks, kwargs...)
Construct a `Simulation` for a `model` with a time stepping `method`,
initial conditions for the prognostic state `Y_init`, initial values of the
auxiliary state `p_init`, time step `Δt` for `tspan` time interval.

Note that if `Y_init` is supplied, `Ya_init` must also be supplied. If `Y_init` is
not supplied, a default state is created for both `Y_init` and `Ya_init`.
"""
function Simulation(
    model::AbstractLandModel,
    method;
    Y_init,
    dt,
    tspan,
    Ya_init,
    callbacks = nothing,
    kwargs...,
)

    # inital state is either default or set externally
    if Y_init isa Nothing
        println(
            "Creating default state for both prognostic and auxiliary states...",
        )
        Y, Ya = default_land_initial_conditions(model)
    else
        Y, Ya = Y_init, Ya_init
    end

    # contains all information about the 
    # pde systems jacobians and right-hand sides
    # to hook into the DiffEqBase.jl interface
    ode_function = make_rhs(model)

    # we use the DiffEqBase.jl interface
    # to set up and an ODE integrator that handles
    # integration in time and callbacks
    ode_problem = DiffEqBase.ODEProblem(ode_function, Y, tspan, Ya)
    integrator = DiffEqBase.init(
        ode_problem,
        method;
        dt = dt,
        callback = callbacks,
        kwargs...,
    )

    return Simulation(model, integrator, callbacks)
end

"""
    step!(simulation::AbstractSimulation, args...; kwargs...)
Step forward a `simulation` one time step.
"""
step!(simulation::AbstractSimulation, args...; kwargs...) =
    DiffEqBase.step!(simulation.integrator, args...; kwargs...)

"""
    run!(simulation::AbstractSimulation, args...; kwargs...)
Run a `simulation` to the end.
"""
run!(simulation::AbstractSimulation, args...; kwargs...) =
    DiffEqBase.solve!(simulation.integrator, args...; kwargs...)
