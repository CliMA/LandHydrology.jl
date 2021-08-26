export SoilIC, create_initial_state
using RecursiveArrayTools: ArrayPartition
abstract type AbstractModelIC end

struct SoilIC{H,E} <: AbstractModelIC
    hydrology::H
    energy::E
end

function SoilIC(;energy = energy, hydrology = hydrology)
    args = (hydrology, energy)
    return SoilIC{typeof.(args)...}(args...)
end

function create_initial_state(model::SoilModel, ic::SoilIC, zc)
    Y_hydrology = ic.hydrology.(zc, Ref(model))
    Y_energy = ic.energy.(zc, Ref(model))
    return ArrayPartition(Y_hydrology, Y_energy)
end
