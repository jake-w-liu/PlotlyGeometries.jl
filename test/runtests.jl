using PlotlyGeometries
using PlotlySupply
using Test

@testset "PlotlyGeometries.jl" begin
    @testset "Base Geometry + Plot Mutation" begin
        c = cuboids([0, 0, 0], [1, 2, 3], "red"; opc=0.5)
        @test c[:type] == "mesh3d"

        gtrans!(c, [1, 2, 3])
        @test isapprox(c.x[1], 0.5; atol=1e-9)
        @test isapprox(c.y[1], 1.0; atol=1e-9)
        @test isapprox(c.z[1], 1.5; atol=1e-9)

        fig = Plot([c], blank_layout())
        add_ref_axes!(fig)
        @test length(fig.data) == 7

        add_arrows!(fig, [0, 0, 0], [1, 0, 0], 1.0, "black")
        add_text!(fig, [0, 0, 0], "origin", "black")
        @test length(fig.data) == 10

        set_view!(fig, 45, 30)
        @test haskey(fig.layout, :scene)
        @test haskey(fig.layout.scene, :camera)
        @test haskey(fig.layout.scene[:camera], :eye)
    end

    @testset "New Primitives" begin
        cyl = cylinders([0, 0, 0], 1.0, 2.0, "z", "royalblue"; tres=24)
        fru = cones([0, 0, 0], 1.0, 0.5, 2.0, "z", "orange"; tres=24)
        con = cones([0, 0, 0], 1.0, 2.0, "z", "orange"; tres=24)
        tor = tori([0, 0, 0], 3.0, 1.0, "z", "teal"; ures=24, vres=12)
        pln = planes([0, 0, 0], [2.0, 3.0], "y", "gray")
        dsk = disks([0, 0, 0], 1.5, "x", "purple"; tres=24)

        for geo in (cyl, fru, con, tor, pln, dsk)
            @test geo[:type] == "mesh3d"
            @test length(geo.i) > 0
            @test length(geo.j) == length(geo.i)
            @test length(geo.k) == length(geo.i)
        end
    end

    @testset "Concave Polygon Triangulation" begin
        pts = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [1.0, 0.3, 0.0], # concave notch
            [0.0, 1.0, 0.0],
        ]
        poly = polygons(pts, "cyan")
        @test poly[:type] == "mesh3d"
        @test length(poly.i) == length(pts) - 2
    end

    @testset "Group Transforms + Edge Guards" begin
        l1 = lines([0, 0, 0], [1, 0, 0], "red")
        l2 = lines([0, 1, 0], [1, 1, 0], "red")
        geos = [l1, l2]

        gtrans!(geos, [1, 0, 0])
        @test l1.x == [1, 2]
        @test l2.x == [1, 2]

        grot!(geos, 90, [0, 0, 1], [0, 0, 0])
        @test isapprox(l1.x[1], 0.0; atol=1e-9)
        @test isapprox(l1.y[1], 1.0; atol=1e-9)

        @test_throws AssertionError add_arrows!(Plot([cuboids([0, 0, 0], [1, 1, 1])], blank_layout()), [0, 0, 0], [0, 0, 0], 1.0, "black")
        @test_throws AssertionError grot!(cuboids([0, 0, 0], [1, 1, 1]), 30, [0, 0, 0], [0, 0, 0])
        @test_throws AssertionError polygons([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], 3, "blue")
    end
end
