using RecursiveArrayTools
export init_prognostic_vars
function init_prognostic_vars(
    soil::SoilModel,
    cell_center_space::Spaces.CenterFiniteDifferenceSpace,
)
    zc = Fields.coordinate_field(cell_center_space)
    hydrology = soil.hydrology_model
    energy = soil.energy_model
    msp = soil.soil_params
    param_set = soil.earth_params
    prognostic_hydrology = init_prognostic_vars(hydrology, msp, param_set, zc)
    prognostic_energy = init_prognostic_vars(energy, msp, param_set, zc)

    return ArrayPartition(prognostic_hydrology, prognostic_energy)
end

function init_prognostic_vars(
    hydrology::SoilHydrologyModel,
    soil_params,
    param_set,
    zc,
)
    return hydrology.initial_conditions.(zc)
end

function init_prognostic_vars(
    energy::PrescribedTemperatureModel,
    soil_params,
    param_set,
    zc,
)
    return ()
end

function init_prognostic_vars(
    energy::SoilEnergyModel,
    soil_params,
    param_set,
    zc,
)
    return energy.initial_conditions.(zc, Ref(soil_params), Ref(param_set))
end
