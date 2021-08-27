export SoilIC, create_initial_state, IC, NoIC
using RecursiveArrayTools: ArrayPartition

abstract type AbstractModelIC end
struct NoIC <: AbstractModelIC end

struct IC{F} <: AbstractModelIC
    f::F
end

Base.@kwdef struct SoilIC{h<: AbstractModelIC, e<: AbstractModelIC}
    hydrology::h = NoIC()
    energy::e = NoIC()
end

# could also dispatch on type of model, rather than type of IC?
function create_initial_state(model::AbstractSoilModel, ic::SoilIC, zc)
    Y_hydrology = create_initial_state(ic.hydrology, zc, model)
    Y_energy = create_initial_state(ic.energy, zc, model)
    f(x) = length(x) > 0
    prognostic_vars = filter(f, (Y_hydrology, Y_energy))
    # assert not empty?
    return ArrayPartition(prognostic_vars)
end

function create_initial_state(ic::NoIC, zc, model)
    return ()
end

function create_initial_state(ic::IC, zc, model)
    return ic.f.(zc, Ref(model))
end
