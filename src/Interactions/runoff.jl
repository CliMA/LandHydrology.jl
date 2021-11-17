export compute_infiltration

"""
    compute_infiltration(model::AbstractModel,
                         soil_water::SoilHydrologyModel{FT},
                         lm::LandHydrologyModel,
                         Y::Fields.FieldVector,
                         Ya::Fields.FieldVector,
                         t::FT,
                         ) where {FT}

Computes the infiltration given the surface water state, the soil moisture
in the top layer, and precipitation rate.
"""
function compute_infiltration(model::AbstractModel,
                              soil_water::SoilHydrologyModel{FT},
                              lm::LandHydrologyModel,
                              Y::Fields.FieldVector,
                              Ya::Fields.FieldVector,
                              t::FT,
                              ) where {FT}
    @unpack S_s, ν  = lm.soil.soil_param_set
    precip = lm.atmos_state.precipitation
    surface_height = get_surface_height(model, Y)
    # Compute infiltration capacity
    ϑ_l⁺ = ν + S_s * surface_height
    @unpack fϑ_l = boundary_fluxes(Y,Ya,SoilComponentBC(hydrology = Dirichlet((t) -> ϑ_l⁺)), :top, lm, axes(Y.soil.ϑ_l), t)
    if precip < fϑ_l # more negative
        infiltration = fϑ_l
    else
        infiltration = precip
    end
    return infiltration
end


