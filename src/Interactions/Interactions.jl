module Interactions
import LandHydrology: SubComponentModels, make_update_aux, AbstractModel
using LandHydrology: LandHydrologyModel, AbstractAtmosState, NoAtmosState
using LandHydrology.SoilInterface: SoilModel, SoilHydrologyModel, Dirichlet, SoilComponentBC, boundary_fluxes
using LandHydrology.SurfaceWater: SurfaceWaterModel, get_surface_height
using ClimaCore: Fields
using UnPack

"""
     SubComponentModels.make_update_aux(atmos_state::PrescribedAtmosState,
                                            soil::SoilModel,
                                            sfc_water::AbstractModel,
                                            model::LandHydrologyModel
                                            )

Constructs the function `update_aux` which updates the auxiliary state
`Ya.soil_infiltration` for the case of surface water/ground water interactions,
assuming a prescribed precipitation value.

This variable is used by the soil model as a boundary condition, and by the
surface water model as a source term (runoff). 
"""
function SubComponentModels.make_update_aux(atmos_state::PrescribedAtmosState,
                                            soil::SoilModel,
                                            sfc_water::AbstractModel,
                                            model::LandHydrologyModel
                                            )
    function update_aux!(Ya,Y,t)
        Ya.soil_infiltration = compute_infiltration(sfc_water, soil.hydrology_model, model, Y, Ya, t)
    end
end

"""
     SubComponentModels.make_update_aux(atmos_state::AbstractAtmosState,
                                            soil::SoilModel,
                                            sfc_water::AbstractModel,
                                            model::LandHydrologyModel
                                            )
**** THink about this.
Constructs the function `update_aux` which updates the auxiliary state,
when no interaction is required.
"""
function SubComponentModels.make_update_aux(atmos_state::AbstractAtmosState, soil::SoilModel, sfc_water::AbstractModel, model::LandHydrologyModel)
    function update_aux!(Ya,Y,t)
        nothing
    end
end

include("runoff.jl")

end        
