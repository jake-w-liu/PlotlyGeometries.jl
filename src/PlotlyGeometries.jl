module PlotlyGeometries

using PlotlyJS
using Combinatorics
using LinearAlgebra

export cubes, squares, polygons, ellipsoids, spheres, lines, create_mesh, 
    rot!, trans!, sort_pts!, add_ref_axes, add_arrows, add_text,
    blank_layout

include("apis.jl")

end
