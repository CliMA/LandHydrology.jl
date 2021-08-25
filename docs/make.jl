using Documenter
using LandHydrology

makedocs(
    sitename = "LandHydrology.jl",
    authors = "CliMA Land Model Team"
    format = Documenter.HTML(collapselevel = 1, mathengine = MathJax3()),
    pages = [
    "Home" => "index.md",
    ],
    modules = [LandHydrology]
)

deploydocs(repo = "github.com/CliMA/LandHydrology.jl.git", devbranch = "main")
