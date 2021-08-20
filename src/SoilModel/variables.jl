using RecursiveArrayTools

function init_prognostic_vars(soil::SoilModel, cell_center_space::Spaces.CenterFiniteDifferenceSpace)
    zc = Fields.coordinate_field(cell_center_space)
    hydrology = soil.hydrology_model
    energy = soil.energy_model
    prognostic_hydrology = init_prognostic_vars(hydrology, zc)
    prognostic_energy = init_prognostic_vars(energy, zc)
    #variables = [v for v in (prognostic_hydrology, prognostic_energy) if v !== nothing]
    return ArrayPartition(prognostic_hydrology, prognostic_energy)
end

function init_prognostic_vars(hydrology::SoilHydrologyModel, zc)
    return hydrology.initial_conditions.(zc)
end

function init_prognostic_vars(energy::PrescribedTemperatureModel, zc)
   return nothing 
end

#function init_prognostic_vars(energy::SoilEnergyModel, zc)
#    return energy.initial_conditions.(zc) 
#end

    
