push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using Documenter
using LandHydrology
include("list_of_apis.jl")

pages = Any[
    "Home" => "index.md",
    "APIs" => apis,
    "Contribution guide" => "Contributing.md",
]

mathengine = MathJax(
    Dict(
        :TeX => Dict(
            :equationNumbers => Dict(:autoNumber => "AMS"),
            :Macros => Dict(),
        ),
    ),
)


format = Documenter.HTML(
    prettyurls = !isempty(get(ENV, "CI", "")),
    collapselevel = 1,
    mathengine = mathengine,
)

makedocs(
    sitename = "LandHydrology.jl",
    authors = "Clima Land Model Team",
    format = format,
    pages = pages,
    checkdocs = :exports,
    doctest = true,
    strict = false,
    clean = true,
    modules = [LandHydrology],
)

deploydocs(
    repo = "github.com/CliMA/LandHydrology.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)
