var documenterSearchIndex = {"docs":
[{"location":"APIs/SoilModel/SoilWaterParameterizations/#Soil-Water-Parameterizations","page":"Soil Water Parameterizations","title":"Soil Water Parameterizations","text":"","category":"section"},{"location":"APIs/SoilModel/SoilWaterParameterizations/","page":"Soil Water Parameterizations","title":"Soil Water Parameterizations","text":"CurrentModule = LandHydrology.SoilInterface.SoilWaterParameterizations","category":"page"},{"location":"APIs/SoilModel/SoilWaterParameterizations/","page":"Soil Water Parameterizations","title":"Soil Water Parameterizations","text":"vanGenuchten\nvolumetric_liquid_fraction\neffective_saturation\npressure_head\nhydraulic_conductivity\nhydrostatic_profile\nmatric_potential\ninverse_matric_potential\nviscosity_factor\nimpedance_factor\nAbstractConductivityFactor\nIceImpedance\nTemperatureDependentViscosity\nNoEffect","category":"page"},{"location":"APIs/SoilModel/SoilWaterParameterizations/#LandHydrology.SoilInterface.SoilWaterParameterizations.vanGenuchten","page":"Soil Water Parameterizations","title":"LandHydrology.SoilInterface.SoilWaterParameterizations.vanGenuchten","text":"vanGenuchten{FT} <: AbstractHydraulicsModel{FT}\n\nThe necessary parameters for the van Genuchten hydraulic model;  defaults are for Yolo light clay.\n\nFields\n\nn\nPore size distribution index, unitless\nα\nInverse of the air entry potential (1/m)\nm\nExponent parameter - determined as 1-1/n\nθr\nResidual water content\nKsat\nSaturated hydraulic conductivity, m/s\n\n\n\n\n\n","category":"type"},{"location":"APIs/SoilModel/SoilWaterParameterizations/#LandHydrology.SoilInterface.SoilWaterParameterizations.volumetric_liquid_fraction","page":"Soil Water Parameterizations","title":"LandHydrology.SoilInterface.SoilWaterParameterizations.volumetric_liquid_fraction","text":"volumetric_liquid_fraction(\n    ϑ_l::FT,\n    ν_eff::FT,\n) where {FT}\n\nComputes the volumetric liquid fraction from the effective porosity and the augmented liquid fraction.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilWaterParameterizations/#LandHydrology.SoilInterface.SoilWaterParameterizations.effective_saturation","page":"Soil Water Parameterizations","title":"LandHydrology.SoilInterface.SoilWaterParameterizations.effective_saturation","text":"effective_saturation(porosity::FT, ϑ_l::FT, θr::FT)\n\nComputes the effective saturation of the soil, given the augmented liquid fraction ϑ_l, the residual water fraction θr, and a porosity.\n\nThis porosity can be effective (ν-θ_i) or true (ν) depending on the application. Note that ϑl is defined to be larger than θr. If `ϑl-θr` is negative,  hydraulic functions that take it as an argument will return  imaginary numbers, resulting in domain errors.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilWaterParameterizations/#LandHydrology.SoilInterface.SoilWaterParameterizations.pressure_head","page":"Soil Water Parameterizations","title":"LandHydrology.SoilInterface.SoilWaterParameterizations.pressure_head","text":"pressure_head(hm::vanGenuchten{FT}, ϑ_l::FT, ν_eff::FT, S_s::FT)\n\nComputes the pressure head according to the hydraulics model hm.\n\nFor unsaturated soil (ϑ_l < ν_eff, or S_l_eff < 1), the matric potential is used. It is evaluated at an effective saturation computed using the effective porosity. This ensures that the pressure head is continuous across the saturated and unsaturated  border. For saturated soil, the pressure head is positive.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilWaterParameterizations/#LandHydrology.SoilInterface.SoilWaterParameterizations.hydraulic_conductivity","page":"Soil Water Parameterizations","title":"LandHydrology.SoilInterface.SoilWaterParameterizations.hydraulic_conductivity","text":"hydraulic_conductivity(hm::vanGenuchten{FT}, S::FT, viscosity_f::FT, impedance_f::FT) where {FT}\n\nComputes the hydraulic conductivity given a hydraulics model hm,  the effective saturation, and the additional  viscosity and ice impedance factors, which must be computed prior  to evaluating this function.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilWaterParameterizations/#LandHydrology.SoilInterface.SoilWaterParameterizations.hydrostatic_profile","page":"Soil Water Parameterizations","title":"LandHydrology.SoilInterface.SoilWaterParameterizations.hydrostatic_profile","text":"hydrostatic_profile(hm, z, zmin, ν)\n\nReturns the augmented liquid fraction corresponding to a hydrostatic profile in variably saturated soil, with the boundary between the saturated and unsaturated zones at z_∇.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilWaterParameterizations/#LandHydrology.SoilInterface.SoilWaterParameterizations.matric_potential","page":"Soil Water Parameterizations","title":"LandHydrology.SoilInterface.SoilWaterParameterizations.matric_potential","text":"matric_potential(hm::vanGenuchten{FT}, S)\n\nComputes the matric potential according the the hydraulics model hm, given an effective saturation S.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilWaterParameterizations/#LandHydrology.SoilInterface.SoilWaterParameterizations.inverse_matric_potential","page":"Soil Water Parameterizations","title":"LandHydrology.SoilInterface.SoilWaterParameterizations.inverse_matric_potential","text":"inverse_matric_potential(\n    hm::vanGenuchten{FT},\n    ψ::FT\n) where {FT}\n\nCompute the effective saturation given the matric potential ψ, using the van Genuchten formulation hm.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilWaterParameterizations/#LandHydrology.SoilInterface.SoilWaterParameterizations.viscosity_factor","page":"Soil Water Parameterizations","title":"LandHydrology.SoilInterface.SoilWaterParameterizations.viscosity_factor","text":"viscosity_factor(\n    vm::NoEffect{FT},\n    _...,\n) where {FT}\n\nReturns a multiplicative factor of 1.0, which does not account for the temperature dependence of the conductivity.\n\n\n\n\n\nviscosity_factor(\n    vm::TemperatureDependentViscosity{FT},\n    T::FT,\n) where {FT}\n\nReturns the multiplicative factor which accounts for the temperature dependence of the conductivity.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilWaterParameterizations/#LandHydrology.SoilInterface.SoilWaterParameterizations.impedance_factor","page":"Soil Water Parameterizations","title":"LandHydrology.SoilInterface.SoilWaterParameterizations.impedance_factor","text":"impedance_factor(\n    imp::NoEffect{FT},\n    _...,\n) where {FT}\n\nReturns a multiplicative factor of 1.0, which does not account for the ice on conductivity.\n\n\n\n\n\nimpedance_factor(\n    imp::IceImpedance{FT},\n    f_i::FT,\n) where {FT}\n\nReturns the multiplicative factor when an effect due to the fraction of  ice f_i is desired. \n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilWaterParameterizations/#LandHydrology.SoilInterface.SoilWaterParameterizations.AbstractConductivityFactor","page":"Soil Water Parameterizations","title":"LandHydrology.SoilInterface.SoilWaterParameterizations.AbstractConductivityFactor","text":"AbstractConductivityFactor{FT <: AbstractFloat}\n\nAn abstract type for multiplicative factors that affect the hydraulic conductivity, increasing or decreasing it relative to the formula that depends only on liquid water content.\n\nExamples include a multiplicative factor due to ice impedance or one due to the  temperature dependence of the viscosity of water.\n\n\n\n\n\n","category":"type"},{"location":"APIs/SoilModel/SoilWaterParameterizations/#LandHydrology.SoilInterface.SoilWaterParameterizations.IceImpedance","page":"Soil Water Parameterizations","title":"LandHydrology.SoilInterface.SoilWaterParameterizations.IceImpedance","text":"IceImpedance{FT} <: AbstractConductivityFactor{FT}\n\nThe factor type which alters the conductivity to account for the effects of ice in the pore space.\n\nThe empirical factor comes from Lundin (1990).\n\n\n\n\n\n","category":"type"},{"location":"APIs/SoilModel/SoilWaterParameterizations/#LandHydrology.SoilInterface.SoilWaterParameterizations.TemperatureDependentViscosity","page":"Soil Water Parameterizations","title":"LandHydrology.SoilInterface.SoilWaterParameterizations.TemperatureDependentViscosity","text":"TemperatureDependentViscosity{FT} <: AbstractConductivityFactor{FT}\n\nThe factor type which alters the conductivity to account for the effects of temperature on the viscosity of water.\n\n\n\n\n\n","category":"type"},{"location":"APIs/SoilModel/SoilWaterParameterizations/#LandHydrology.SoilInterface.SoilWaterParameterizations.NoEffect","page":"Soil Water Parameterizations","title":"LandHydrology.SoilInterface.SoilWaterParameterizations.NoEffect","text":"NoEffect{FT} <: AbstractConductivityFactor{FT}\n\nThe default type of conductivity factor, used to set a particular multiplicative factor to unity.\n\n\n\n\n\n","category":"type"},{"location":"APIs/SoilModel/SoilHeatParameterizations/#Soil-Heat-Parameterizations","page":"Soil Heat Parameterizations","title":"Soil Heat Parameterizations","text":"","category":"section"},{"location":"APIs/SoilModel/SoilHeatParameterizations/","page":"Soil Heat Parameterizations","title":"Soil Heat Parameterizations","text":"CurrentModule = LandHydrology.SoilInterface.SoilHeatParameterizations","category":"page"},{"location":"APIs/SoilModel/SoilHeatParameterizations/#Heat-functions","page":"Soil Heat Parameterizations","title":"Heat functions","text":"","category":"section"},{"location":"APIs/SoilModel/SoilHeatParameterizations/","page":"Soil Heat Parameterizations","title":"Soil Heat Parameterizations","text":"volumetric_heat_capacity\nvolumetric_internal_energy\nsaturated_thermal_conductivity\nthermal_conductivity\nrelative_saturation\nkersten_number\nk_solid\nk_dry\nksat_unfrozen\nksat_frozen\nvolumetric_internal_energy_liq\ntemperature_from_ρe_int","category":"page"},{"location":"APIs/SoilModel/SoilHeatParameterizations/#LandHydrology.SoilInterface.SoilHeatParameterizations.volumetric_heat_capacity","page":"Soil Heat Parameterizations","title":"LandHydrology.SoilInterface.SoilHeatParameterizations.volumetric_heat_capacity","text":"volumetric_heat_capacity(\n    θ_l::FT,\n    θ_i::FT,\n    ρc_ds::FT,\n    param_set::AbstractParameterSet\n) where {FT}\n\nCompute the expression for volumetric heat capacity.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilHeatParameterizations/#LandHydrology.SoilInterface.SoilHeatParameterizations.volumetric_internal_energy","page":"Soil Heat Parameterizations","title":"LandHydrology.SoilInterface.SoilHeatParameterizations.volumetric_internal_energy","text":"volumetric_internal_energy(\n    θ_i::FT,\n    ρc_s::FT,\n    T::FT,\n    param_set::AbstractParameterSet\n) where {FT}\n\nCompute the expression for volumetric internal energy.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilHeatParameterizations/#LandHydrology.SoilInterface.SoilHeatParameterizations.saturated_thermal_conductivity","page":"Soil Heat Parameterizations","title":"LandHydrology.SoilInterface.SoilHeatParameterizations.saturated_thermal_conductivity","text":"saturated_thermal_conductivity(\n    θ_l::FT,\n    θ_i::FT,\n    κ_sat_unfrozen::FT,\n    κ_sat_frozen::FT\n) where {FT}\n\nCompute the expression for saturated thermal conductivity of soil matrix.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilHeatParameterizations/#LandHydrology.SoilInterface.SoilHeatParameterizations.thermal_conductivity","page":"Soil Heat Parameterizations","title":"LandHydrology.SoilInterface.SoilHeatParameterizations.thermal_conductivity","text":"thermal_conductivity(\n    κ_dry::FT,\n    K_e::FT,\n    κ_sat::FT\n) where {FT}\n\nCompute the expression for thermal conductivity of soil matrix.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilHeatParameterizations/#LandHydrology.SoilInterface.SoilHeatParameterizations.relative_saturation","page":"Soil Heat Parameterizations","title":"LandHydrology.SoilInterface.SoilHeatParameterizations.relative_saturation","text":"relative_saturation(\n        θ_l::FT,\n        θ_i::FT,\n        porosity::FT\n) where {FT}\n\nCompute the expression for relative saturation.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilHeatParameterizations/#LandHydrology.SoilInterface.SoilHeatParameterizations.kersten_number","page":"Soil Heat Parameterizations","title":"LandHydrology.SoilInterface.SoilHeatParameterizations.kersten_number","text":"kersten_number(\n    θ_i::FT,\n    S_r::FT,\n    soil_param_functions::PS\n) where {FT, PS}\n\nCompute the expression for the Kersten number.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilHeatParameterizations/#LandHydrology.SoilInterface.SoilHeatParameterizations.k_solid","page":"Soil Heat Parameterizations","title":"LandHydrology.SoilInterface.SoilHeatParameterizations.k_solid","text":"function k_solid(\n    ν_ss_om::FT,\n    ν_ss_quartz::FT,\n    κ_quartz::FT,\n    κ_minerals::FT,\n    κ_om::FT,\n) where {FT}\n\nComputes the thermal conductivity of the solid material in soil. The _ss_ subscript denotes that the volumetric fractions of the soil components are referred to the soil solid components, not including the pore space.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilHeatParameterizations/#LandHydrology.SoilInterface.SoilHeatParameterizations.k_dry","page":"Soil Heat Parameterizations","title":"LandHydrology.SoilInterface.SoilHeatParameterizations.k_dry","text":"function k_dry(\n    param_set::AbstractParameterSet\n    soil_param_functions::PS,\n) where {PS}\n\nComputes the thermal conductivity of dry soil.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilHeatParameterizations/#LandHydrology.SoilInterface.SoilHeatParameterizations.ksat_unfrozen","page":"Soil Heat Parameterizations","title":"LandHydrology.SoilInterface.SoilHeatParameterizations.ksat_unfrozen","text":"function ksat_unfrozen(\n    κ_solid::FT,\n    porosity::FT,\n    κ_l::FT\n) where {FT}\n\nComputes the thermal conductivity for saturated unfrozen soil.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilHeatParameterizations/#LandHydrology.SoilInterface.SoilHeatParameterizations.ksat_frozen","page":"Soil Heat Parameterizations","title":"LandHydrology.SoilInterface.SoilHeatParameterizations.ksat_frozen","text":"function ksat_frozen(\n    κ_solid::FT,\n    porosity::FT,\n    κ_ice::FT\n) where {FT}\n\nComputes the thermal conductivity for saturated frozen soil.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilHeatParameterizations/#LandHydrology.SoilInterface.SoilHeatParameterizations.volumetric_internal_energy_liq","page":"Soil Heat Parameterizations","title":"LandHydrology.SoilInterface.SoilHeatParameterizations.volumetric_internal_energy_liq","text":"volumetric_internal_energy_liq(\n    T::FT,\n    T_ref::FT,\n) where {FT}\n\nCompute the expression for the volumetric internal energy of liquid water.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilHeatParameterizations/#LandHydrology.SoilInterface.SoilHeatParameterizations.temperature_from_ρe_int","page":"Soil Heat Parameterizations","title":"LandHydrology.SoilInterface.SoilHeatParameterizations.temperature_from_ρe_int","text":"function temperature_from_ρe_int(\n    ρe_int::FT,\n    θ_i::FT,\n    ρc_s::FT,\n    param_set::AbstractParameterSet\n) where {FT}\n\nComputes the temperature of soil given θ_i and volumetric internal energy ρe_int.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilInterface/#Soil-Interface","page":"SoilInterface","title":"Soil Interface","text":"","category":"section"},{"location":"APIs/SoilModel/SoilInterface/","page":"SoilInterface","title":"SoilInterface","text":"CurrentModule = LandHydrology.SoilInterface","category":"page"},{"location":"APIs/SoilModel/SoilInterface/#Parameters","page":"SoilInterface","title":"Parameters","text":"","category":"section"},{"location":"APIs/SoilModel/SoilInterface/","page":"SoilInterface","title":"SoilInterface","text":"SoilParams","category":"page"},{"location":"APIs/SoilModel/SoilInterface/#LandHydrology.SoilInterface.SoilParams","page":"SoilInterface","title":"LandHydrology.SoilInterface.SoilParams","text":"struct SoilParams{FT} <: ParameterStructure{FT}\n\n\n\n\n\n","category":"type"},{"location":"APIs/SoilModel/SoilInterface/#Models","page":"SoilInterface","title":"Models","text":"","category":"section"},{"location":"APIs/SoilModel/SoilInterface/","page":"SoilInterface","title":"SoilInterface","text":"SoilModel\nSoilHydrologyModel\nSoilEnergyModel\nPrescribedTemperatureModel\nPrescribedHydrologyModel","category":"page"},{"location":"APIs/SoilModel/SoilInterface/#LandHydrology.SoilInterface.SoilModel","page":"SoilInterface","title":"LandHydrology.SoilInterface.SoilModel","text":"SoilModel{FT, domain, em <: AbstractSoilModel, hm <: AbstractSoilModel, bc, A,B}\n\nThe model type for the soil model.\n\nFields\n\ndomain\nenergy_model\nSoil energy model - prescribed or dynamics\nhydrology_model\nSoil hydrology model - prescribed or dynamic\nboundary_conditions\nBoundary conditions tuple\nsoil_param_set\nSoil parameters\nearth_param_set\nEarth parameter set\nname\nname\nvariables\nvariables\n\n\n\n\n\n","category":"type"},{"location":"APIs/SoilModel/SoilInterface/#LandHydrology.SoilInterface.SoilHydrologyModel","page":"SoilInterface","title":"LandHydrology.SoilInterface.SoilHydrologyModel","text":"SoilHydrologyModel\n\nThe model type to be used when the user wants to simulate the flow of water in soil by solving Richards equation.\n\nFields\n\nhydraulic_model\nviscosity_factor\nimpedance_factor\n\n\n\n\n\n","category":"type"},{"location":"APIs/SoilModel/SoilInterface/#LandHydrology.SoilInterface.SoilEnergyModel","page":"SoilInterface","title":"LandHydrology.SoilInterface.SoilEnergyModel","text":"SoilEnergyModel\n\nThe model type to be used when the user wants to simulate heat transfer in soil by solving the heat partial differential equation.\n\nFields\n\n\n\n\n\n","category":"type"},{"location":"APIs/SoilModel/SoilInterface/#LandHydrology.SoilInterface.PrescribedTemperatureModel","page":"SoilInterface","title":"LandHydrology.SoilInterface.PrescribedTemperatureModel","text":"Base.@kwdef struct PrescribedTemperatureModel <: AbstractSoilComponentModel\n    \"Profile of (z,t) for temperature\"\n     T_profile::Function = (z,t) -> eltype(z)(288)\nend\n\nThe model type to be used when the user does not wish to solve the heat partial differential equation, but instead wishes to prescibe a temperature profile in the soil.\n\nThis is useful for situations where Richards Equation alone is sufficient. Because the hydraulic conductivity can be a function of temperature, a temperature profile can be supplied in order to simulate that. The default is 288K across the domain, the reference temperature for the viscosity effect.\n\nFields\n\nT_profile\nProfile of (z,t) for temperature\n\n\n\n\n\n","category":"type"},{"location":"APIs/SoilModel/SoilInterface/#LandHydrology.SoilInterface.PrescribedHydrologyModel","page":"SoilInterface","title":"LandHydrology.SoilInterface.PrescribedHydrologyModel","text":"Base.@kwdef struct PrescribedHydrologyModel <: AbstractSoilComponentModel\n    \"Profile of (z,t) for ϑ_l\"\n    ϑ_l_profile::Function = (z,t) -> eltype(z)(0.0)\n    \"Profile of (z,t) for θ_i\"\n    θ_i_profile::::Function = (z,t) -> eltype(z)(0.0)\nend\n\nThe model type to be used when the user does not wish to solve the Richards equation, but instead wishes to prescibe a water content profile in the soil.\n\nThis is useful for situations where only the heat equation is to be solved. Because the thermal conductivity and heat capacities depend on water content, a water profile must be defined to solve the heat equation. The default for both ice and liquid water content is zero, which applies for totally dry soil. \n\nFields\n\nϑ_l_profile\nProfile of (z,t) for ϑ_l\nθ_i_profile\nProfile of (z,t) for θ_i\n\n\n\n\n\n","category":"type"},{"location":"APIs/SoilModel/SoilInterface/#Right-hand-side","page":"SoilInterface","title":"Right hand side","text":"","category":"section"},{"location":"APIs/SoilModel/SoilInterface/","page":"SoilInterface","title":"SoilInterface","text":"make_rhs","category":"page"},{"location":"APIs/SoilModel/SoilInterface/#LandHydrology.SoilInterface.make_rhs","page":"SoilInterface","title":"LandHydrology.SoilInterface.make_rhs","text":"make_rhs(model::SoilModel)\n\nA function which takes a model::AbstractModel as argument,  and returns function which computes the rhs of a set of ordinary differential equations corresponding to that model.\n\nCurrently, the arguments of the returned rhs function are configured for use with OrdinaryDiffEq.jl. For the soil model, the rhs function depends on the type of the  components of the model (the energy and hydrology models), as well as  whether additional sources are included.\n\n\n\n\n\n","category":"function"},{"location":"APIs/SoilModel/SoilInterface/#Boundary-Conditions","page":"SoilInterface","title":"Boundary Conditions","text":"","category":"section"},{"location":"APIs/SoilModel/SoilInterface/","page":"SoilInterface","title":"SoilInterface","text":"VerticalFlux\nNoBC\nFreeDrainage\nDirichlet\nSoilDomainBC\nSoilComponentBC","category":"page"},{"location":"APIs/SoilModel/SoilInterface/#LandHydrology.SoilInterface.VerticalFlux","page":"SoilInterface","title":"LandHydrology.SoilInterface.VerticalFlux","text":"VerticalFlux{f <: AbstractFloat} <: AbstractBC\n\nThe BC type to be used for prescribed vertical boundary fluxes. The flux is assumed to be of the form\n\nF = f z\n\nwhere f is the value supplied by the user (currently a constant).\n\nFields\n\nflux\nScalar flux; positive = aligned with ẑ\n\n\n\n\n\n","category":"type"},{"location":"APIs/SoilModel/SoilInterface/#LandHydrology.SoilInterface.NoBC","page":"SoilInterface","title":"LandHydrology.SoilInterface.NoBC","text":"NoBC <: AbstractBC\n\nThe BC type to be used when equations do not require boundary conditions, e.g. for Prescribed Models.\n\n\n\n\n\n","category":"type"},{"location":"APIs/SoilModel/SoilInterface/#LandHydrology.SoilInterface.FreeDrainage","page":"SoilInterface","title":"LandHydrology.SoilInterface.FreeDrainage","text":"FreeDrainage <: AbstractBC\n\nA BC type for use with the SoilHydrologyModel (Richards Equation), at the bottom of the domain, setting \n\nh = 1\n\n\n\n\n\n","category":"type"},{"location":"APIs/SoilModel/SoilInterface/#LandHydrology.SoilInterface.Dirichlet","page":"SoilInterface","title":"LandHydrology.SoilInterface.Dirichlet","text":"struct Dirichlet{f} <: AbstractBC\n\nA BC type setting a boundary value of the state ϑ_l or temperature T. This can be a function of time.\n\nNote that we apply boundary conditions to temperature even though volumetric internal energy is the prognostic variable.\n\nFields\n\nstate_value\nState Value; (t) -> f(t)\n\n\n\n\n\n","category":"type"},{"location":"APIs/SoilModel/SoilInterface/#LandHydrology.SoilInterface.SoilDomainBC","page":"SoilInterface","title":"LandHydrology.SoilInterface.SoilDomainBC","text":"struct SoilDomainBC{D, TBC, BBC}\n\nA container holding the SoilComponentBC for each boundary face.\n\nEach field value should be of type SoilComponentBC. This doesn't do  what we want. Ideally the fields would change depending on the domain -  e.g. for a Column, they are top and bottom, but for a 3D domain, they might be top, bottom, xleft, xright, yleft, yright.\n\nFields\n\ntop\nSoilComponentBC for the top of the domain\nbottom\nSoilComponentBC for the bottom of the domain\n\n\n\n\n\n","category":"type"},{"location":"APIs/SoilModel/SoilInterface/#LandHydrology.SoilInterface.SoilComponentBC","page":"SoilInterface","title":"LandHydrology.SoilInterface.SoilComponentBC","text":"struct SoilComponentBC{ebc <: AbstractBC, hbc <: AbstractBC}\n\nA container for holding the boundary conditions for the components of the soil model.\n\nThe values must be of type AbstractBC; the two components are energy and hydrology. Each boundary will have a SoilComponentBC object associated with it.\n\nFields\n\nenergy\nBC for the heat equation\nhydrology\nBC for ϑ_l\n\n\n\n\n\n","category":"type"},{"location":"APIs/SoilModel/SoilInterface/#Initial-Conditions","page":"SoilInterface","title":"Initial Conditions","text":"","category":"section"},{"location":"APIs/SoilModel/SoilInterface/","page":"SoilInterface","title":"SoilInterface","text":"set_initial_state","category":"page"},{"location":"APIs/SoilModel/SoilInterface/#LandHydrology.SoilInterface.set_initial_state","page":"SoilInterface","title":"LandHydrology.SoilInterface.set_initial_state","text":"set_initial_state(model::SoilModel, f::Function, t0::Real)\n\nA function which sets the initial state for the soil model, given a function of space (and time) f, as well as the initial time t0.\n\n\n\n\n\n","category":"function"},{"location":"#LandHydrology.jl","page":"Home","title":"LandHydrology.jl","text":"","category":"section"},{"location":"Contributing/#Contributing","page":"Contribution guide","title":"Contributing","text":"","category":"section"},{"location":"Contributing/","page":"Contribution guide","title":"Contribution guide","text":"Thank you for contributing to LandHydrology! We encourage Pull Requests (PRs). Please do not hesitate to ask questions.","category":"page"},{"location":"Contributing/#Some-useful-tips","page":"Contribution guide","title":"Some useful tips","text":"","category":"section"},{"location":"Contributing/","page":"Contribution guide","title":"Contribution guide","text":"When you start working on a new feature branch, make sure you start from master by running: git checkout master.\nMake sure you add tests for your code in test/ and appropriate documentation in the code and/or in docs/. All exported functions and structs must be documented.\nWhen your PR is ready for review, clean up your commit history by squashing and make sure your code is current with ClimateMachine master by rebasing.","category":"page"},{"location":"Contributing/#Continuous-integration","page":"Contribution guide","title":"Continuous integration","text":"","category":"section"},{"location":"Contributing/","page":"Contribution guide","title":"Contribution guide","text":"After rebasing your branch, you can ask for review. Fill out the template and provide a clear summary of what your PR does. When a PR is created or updated, a set of automated tests are run on the PR in our continuous integration (CI) system.","category":"page"},{"location":"Contributing/#Automated-testing","page":"Contribution guide","title":"Automated testing","text":"","category":"section"},{"location":"Contributing/","page":"Contribution guide","title":"Contribution guide","text":"Currently a number of checks are run per commit for a given PR.","category":"page"},{"location":"Contributing/","page":"Contribution guide","title":"Contribution guide","text":"JuliaFormatter checks if the PR is formatted with .dev/climaformat.jl.\nDocumentation rebuilds the documentation for the PR and checks if the docs are consistent and generate valid output.\nTests runs the file test/runtests.jl,  using Pkg.test(). These are a mix of unit tests and fast integration tests.","category":"page"},{"location":"Contributing/","page":"Contribution guide","title":"Contribution guide","text":"We use bors to manage merging PR's in the the LandHydrology repo. If you're a collaborator and have the necessary permissions, you can type bors try in a comment on a PR to have integration test suite run on that PR, or bors r+ to try and merge the code.  Bors ensures that all integration tests for a given PR always pass before merging into master.","category":"page"},{"location":"APIs/LandHydrology/#Land-Hydrology","page":"Models","title":"Land Hydrology","text":"","category":"section"},{"location":"APIs/LandHydrology/","page":"Models","title":"Models","text":"CurrentModule = LandHydrology","category":"page"},{"location":"APIs/LandHydrology/#Domains","page":"Models","title":"Domains","text":"","category":"section"},{"location":"APIs/LandHydrology/","page":"Models","title":"Models","text":"LandHydrology.Domains.Column\nLandHydrology.Domains.make_function_space\nLandHydrology.Domains.AbstractVerticalDomain","category":"page"},{"location":"APIs/LandHydrology/#LandHydrology.Domains.Column","page":"Models","title":"LandHydrology.Domains.Column","text":"Column{FT} <: AbstractVerticalDomain\n\nA struct holding the necessary information  to construct a domain, a mesh, a center and face space, etc. For use when a finite difference in 1D is suitable.\n\nFields\n\nzlim\nDomain interval limits, (zmin, zmax)\nnelements\nNumber of elements used to discretize the interval\nboundary_tags\nBoundary face identifiers\n\n\n\n\n\n","category":"type"},{"location":"APIs/LandHydrology/#LandHydrology.Domains.make_function_space","page":"Models","title":"LandHydrology.Domains.make_function_space","text":"make_function_space(domain::Column)\n\nReturns the center and face space z values of the  column domain.\n\n\n\n\n\n","category":"function"},{"location":"APIs/LandHydrology/#LandHydrology.Domains.AbstractVerticalDomain","page":"Models","title":"LandHydrology.Domains.AbstractVerticalDomain","text":"AbstractVerticalDomain{FT}\n\nAn abstract type for vertical domains, using the floating point type FT. \n\n\n\n\n\n","category":"type"},{"location":"APIs/LandHydrology/#Models","page":"Models","title":"Models","text":"","category":"section"},{"location":"APIs/LandHydrology/","page":"Models","title":"Models","text":"LandHydrology.Models.AbstractModel","category":"page"},{"location":"APIs/LandHydrology/#LandHydrology.Models.AbstractModel","page":"Models","title":"LandHydrology.Models.AbstractModel","text":"abstract type AbstractModel\n\nAn abstract type for models.\n\nEventually, the land model and all major subcomponents will be of this type.\n\n\n\n\n\n","category":"type"}]
}
