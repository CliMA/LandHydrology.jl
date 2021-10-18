"""
    struct Simulation <: AbstractSimulation
        # a LandHydrology model
        model::AbstractModel
        # a DiffEqBase.jl integrator used for time
        # stepping the simulation
        integrator::DiffEqBase.DEIntegrator
        # user defined callback operations 
        callbacks::Union{DiffEqBase.CallbackSet, DiffEqBase.DiscreteCallback, Nothing}
    end
A simulation wraps an abstract LandHydrology `model` containing 
equation specifications and an instance of an `integrator` used for
time integration of the discretized model PDE.
"""
struct Simulation{ML <: AbstractModel} <: AbstractSimulation
    model::ML
    integrator::DiffEqBase.DEIntegrator
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
initial conditions `Y_init`, time step `Δt` for `tspan` time interval.
"""
function Simulation(
    model::AbstractModel,
    method;
    Y_init,
    dt,
    tspan,
    p = nothing,
    callbacks = nothing,
    kwargs...,
)

    # inital state is either default or set externally 
    Y = Y_init isa Nothing ? default_initial_conditions(model) : Y_init

    # contains all information about the 
    # pde systems jacobians and right-hand sides
    # to hook into the DiffEqBase.jl interface
    ode_function = make_rhs(model)

    # we use the DiffEqBase.jl interface
    # to set up and an ODE integrator that handles
    # integration in time and callbacks
    ode_problem = DiffEqBase.ODEProblem(ode_function, Y_init, tspan, p)
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
