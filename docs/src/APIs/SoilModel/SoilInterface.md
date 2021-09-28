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
FreeDrainage
Dirichlet
SoilDomainBC
SoilComponentBC
```

## Initial Conditions
```@docs
set_initial_state
```