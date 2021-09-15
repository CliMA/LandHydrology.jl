function instantiate_column(FT)
    domain = Column(FT, zlim = (0.0, 1.0), nelements = 2)
    check1 = domain.zlim == (0.0, 1.0)
    check2 = domain.nelements == 2

    return check1 && check2
end
float_types = (Float32, Float64)
@testset "Domains" begin
    @info "Testing LandHydrology.Domains..."

    @testset "Domains" begin
        for FT in float_types
            @test instantiate_column(FT)

            # Test ndims
            I = Column(FT, zlim = (0.0, 1.0), nelements = 2)
            @test ndims(I) == 1

            # Test length
            I = Column(FT, zlim = (1.0, 2.0), nelements = 2)
            @test length(I) == 1.0

            # Test size
            I = Column(FT, zlim = (1.0, 4.0), nelements = 2)
            @test size(I) == 3.0

            # Test show functions
            I = Column(FT, zlim = (0.0, 1.0), nelements = 2)
            show(I)
            println()
            @test I isa Column{FT}
        end
    end
end
