include("SoilModel/richards_equation.jl")
#Currently these tests dont use the new interface b/c dirichlet BC dont work yet.
include("SoilModel/richards_vg.jl")
include("SoilModel/heat_analytic_unit_test.jl")
#include("SoilModel/dirichlet_bc_as_flux.jl")
