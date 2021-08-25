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
    f(x) = length(x) > 0
    prognostic_vars = filter(f, (prognostic_hydrology, prognostic_energy))
    return ArrayPartition(prognostic_vars)
end

#Alternatively, we might have initial_conditions = (soil_heat = init_energy, soil_water = init_hydrology...) and store in SoilModel instead? Then in the above function we'd still call init_prognostic_variables one by one.

#Another question - would it make more sense to make a single initialize function for all variables? the soil energy/hydrology divide is a little artificial, so it might make sense for that, but less so to have a single function for all land components. 
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
    energy::PrescribedHydrologyModel,
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
