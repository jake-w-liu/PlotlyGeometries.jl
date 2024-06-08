module PlotlyGeometries

using PlotlyJS
using Combinatorics
using LinearAlgebra

export boxes, squares, polygons, 
    ellipsoids, ellipsoids, lines, arrows,
    create_mesh, rotate!, translate!, sort_pts!, add_ref_axes, no_grids, blank_layout

include("apis.jl")

end
