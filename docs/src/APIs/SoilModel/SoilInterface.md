# Soil Interface

```@meta
CurrentModule = LandHydrology.SoilInterface
```

## Parameters
```@docs
SoilParams
```

## Models
```@docs
SoilModel
SoilHydrologyModel
SoilEnergyModel
PrescribedTemperatureModel
PrescribedHydrologyModel
```

## Right hand side
```@docs
make_rhs
```

## Boundary Conditions
```@docs
VerticalFlux
NoBC
SoilDomainBC
SoilComponentBC
compute_vertical_flux
```

## Initial Conditions
```@docs
set_initial_state
```