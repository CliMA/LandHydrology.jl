module Interactions
import LandHydrology: Models, make_update_aux, AbstractModel
using LandHydrology: LandHydrologyModel, AbstractAtmosState, NoAtmosState
using LandHydrology.SoilInterface: SoilModel, SoilHydrologyModel, Dirichlet, SoilComponentBC, boundary_fluxes
using LandHydrology.SurfaceWater: SurfaceWaterModel, get_surface_height
using ClimaCore: Fields
using UnPack

function Models.make_update_aux(atmos_state::AbstractAtmosState, soil::SoilModel, sfc_water::AbstractModel, model::LandHydrologyModel)
    function update_aux!(Ya,Y,t)
        Ya.soil_infiltration = compute_infiltration(sfc_water, soil.hydrology_model, model, Y, Ya, t)
    end
end


function Models.make_update_aux(atmos_state::NoAtmosState, soil::SoilModel, sfc_water::AbstractModel, model::LandHydrologyModel)
    function update_aux!(Ya,Y,t)
        nothing
    end
end

include("runoff.jl")

end        
