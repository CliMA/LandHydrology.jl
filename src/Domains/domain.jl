"""
    Column{FT} <: AbstractVerticalDomain

A struct holding the necessary information 
to construct a domain, a mesh, a center and face
space, etc. For use when a finite difference in
1D is suitable.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Column{FT} <: AbstractVerticalDomain{FT}
    "Domain interval limits, (zmin, zmax)"
    zlim::Tuple{FT, FT}
    "Number of elements used to discretize the interval"
    nelements::Int32
    "Boundary face identifiers"
    boundary_tags::Tuple{Symbol, Symbol}
end

"""
    function Column(FT::DataType = Float64; zlim, nelements)

Outer constructor for the `Column` type.

The `boundary_tags` field values are used to label the boundary faces 
at the top and bottom of the domain.
"""
function Column(FT::DataType = Float64; zlim, nelements)
    @assert zlim[1] < zlim[2]
    boundary_tags = (:bottom, :top)
    return Column{FT}(zlim, nelements, boundary_tags)
end

Base.ndims(::Column) = 1

Base.length(domain::Column) = domain.zlim[2] - domain.zlim[1]

Base.size(domain::Column) = length(domain)

function Base.show(io::IO, domain::Column)
    min = domain.zlim[1]
    max = domain.zlim[2]
    printstyled(io, "[", color = 226)
    astring = @sprintf("%0.1f", min)
    bstring = @sprintf("%0.1f", max)
    printstyled(astring, ", ", bstring, color = 7)
    printstyled(io, "]", color = 226)
end


"""
    make_function_space(domain::Column)

Returns the center and face space z values of the 
column domain.
"""
function make_function_space(domain::Column{FT}) where {FT}
    column = ClimaCore.Domains.IntervalDomain(
        ClimaCore.Geometry.ZPoint{FT}(domain.zlim[1]),
        ClimaCore.Geometry.ZPoint{FT}(domain.zlim[2]);
        boundary_tags = domain.boundary_tags,
    )
    mesh = Meshes.IntervalMesh(column; nelems = domain.nelements)
    center_space = Spaces.CenterFiniteDifferenceSpace(mesh)
    face_space = Spaces.FaceFiniteDifferenceSpace(center_space)

    return center_space, face_space
end
