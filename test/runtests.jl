using PlotlyGeometries
using PlotlyJS
using Test

@testset "PlotlyGeometries.jl" begin
    # Write your tests here.
    b1 = boxes([0, 0, 0], [1, 2, 3], "pink", 0.2)
    fig = plot(b1, blank_layout())
    display(fig)

    str = readline()

end
