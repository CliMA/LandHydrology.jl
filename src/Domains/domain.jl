"""
    Column{FT} <: AbstractVerticalDomain
"""
struct Column{FT} <: AbstractVerticalDomain{FT}
    zlim::Tuple{FT, FT}
    nelements::Int32
end

function Column(FT::DataType = Float64; zlim, nelements)
    @assert zlim[1] < zlim[2]
    return Column{FT}(zlim, nelements)
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
"""
function make_function_space(domain::Column{FT}) where {FT}
    column = ClimaCore.Domains.IntervalDomain(
        domain.zlim[1],
        domain.zlim[2];
        x3boundary = (:bottom, :top),
    )
    mesh = Meshes.IntervalMesh(column; nelems = domain.nelements)
    center_space = Spaces.CenterFiniteDifferenceSpace(mesh)
    face_space = Spaces.FaceFiniteDifferenceSpace(center_space)

    return center_space, face_space
end
