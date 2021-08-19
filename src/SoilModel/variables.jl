function init_prognostic_vars(soil::SoilModel, cell_center_space::Spaces.CenterFiniteDifferenceSpace)
    zc = Fields.coordinate_field(cell_center_space)
    hydrology = soil.hydrology_model
    energy = soil.energy_model
    prognostic_hydrology = init_prognostic_vars(hydrology, zc)
    prognostic_energy = init_prognostic_vars(energy, zc)
    return ArrayPartition(merge(prognostic_hydrology, prognostic_energy))
    
end

function init_prognostic_vars(hydrology::SoilHydrologyModel, zc)
    ϑ_l, θ_i = hydrology.initial_conditions(zc)
    return (ϑ_l = ϑ_l, θ_i = θ_i)
end

function init_prognostic_vars(energy::PrescribedTemperatureModel, zc)
    return ()
end

    
