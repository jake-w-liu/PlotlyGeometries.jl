module PlotlyGeometries

using PlotlyJS
using Combinatorics
using LinearAlgebra

export cubes, squares, ellipsoids, spheres, lines, polygons, 
    rot!, trans!, sort_pts, sort_pts!, add_ref_axes, add_arrows, add_text,
    blank_layout

include("apis.jl")

end
