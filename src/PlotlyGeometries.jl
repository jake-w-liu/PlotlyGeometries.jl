module PlotlyGeometries

using PlotlyJS
using Combinatorics
using LinearAlgebra
using BatchAssign

export cuboids, cubes, squares, ellipsoids, spheres, lines, polygons, 
    grot!, gtrans!, sort_pts, sort_pts!, 
    add_ref_axes!, add_arrows!, add_text!, blank_layout, set_view!

include("api.jl")

end
