# Manual for PlotlyGeometries.jl

## Overview

This manual contains the APIs provided `PlotlyGeometries.jl`. In order to use this module, `PlotlyJS.jl` should be installed first. 

## APIs

### Geometry Creation

#### cubes

```julia
cubes(origin::Vector{<:Real}, dimension::Vector{<:Real}, color::String, opc::Real=1)
```

Creates a 3D box mesh centered at the given origin with specified dimensions and color.


- `origin::Vector{<:Real}`: A vector of three Reals specifying the center of the box.
- `dimension::Vector{<:Real}`: A vector of three Reals specifying the dimensions (width, height, depth) of the box.
- `color::String`: A string specifying the color of the box.
- `opc::Real`: (optional) A Real specifying the opacity of the box. Default is 1.

___

#### squares

```julia
squares(origin::Vector{<:Real}, side::Real, color::String, mode::String="z", opc::Real=1)
```

Creates a 2D square mesh centered at the given origin with the specified side length and color.


- `origin::Vector{<:Real}`: A vector of three Reals specifying the center of the square.
- `side::Real`: A Real specifying the side length of the square.
- `color::String`: A string specifying the color of the square.
- `mode::String`: (optional) A string specifying the orientation of the square ("x", "y", or "z"). Default is "z".
- `opc::Real`: (optional) A Real specifying the opacity of the square. Default is 1.

___

#### ellipsoids

```julia
ellipsoids(origin::Vector{<:Real}, par::Vector{<:Real}, color::String, opc::Real=1, tres=60, pres=30)
```

Creates a 3D ellipsoid mesh.


- `origin::Vector{<:Real}`: The center of the ellipsoid.
- `par::Vector{<:Real}`: Parameters of the ellipsoid (a, b, c).
- `color::String`: The color of the ellipsoid.
- `opc::Real`: The opacity of the ellipsoid. Default is 1.
- `tres`: The resolution of the mesh grid (theta). Default is 60.
- `pres`: The resolution of the mesh grid (phi). Default is 30.

___

#### spheres

```julia
spheres(origin::Vector{<:Real}, r::Real, color::String, opc::Real=1, tres=60, pres=30)
```

Creates a 3D sphere mesh.


- `origin::Vector{<:Real}`: The center of the sphere.
- `r::Real`: Radius of the sphere.
- `color::String`: The color of the sphere.
- `opc::Real`: The opacity of the sphere. Default is 1.
- `tres`: The resolution of the mesh grid (theta). Default is 60.
- `pres`: The resolution of the mesh grid (phi). Default is 30.

___

#### lines

```julia
lines(pt1::Vector{<:Real}, pt2::Vector{<:Real}, color::String, opc::Real=1, style="")
```

Creates a 3D line between two points.


- `pt1::Vector{<:Real}`: Starting point of the line.
- `pt2::Vector{<:Real}`: Ending point of the line.
- `color::String`: The color of the line.
- `opc::Real`: The opacity of the line. Default is 1.
- `style`: The line style (e.g., "solid", "dash"). Default is "".

___

#### polygons

```julia
polygons(pts::Vector, color::String, opc::Real=1)
```

Creates a polygon mesh from a set of points.


- `pts::Vector`: List of points defining the polygon.
- `color::String`: The color of the polygon.
- `opc::Real`: The opacity of the polygon. Default is 1.

```julia
polygons(::Vector, ng::Int, color::String, opc::Real=1)
```

Creates a group of polygons from a set of points and a specified number of vertices per polygon.


- `pts::Vector`: List of points defining the mesh.
- `ng::Int`: Number of vertices per polygon.
- `color::String`: The color of the mesh.
- `opc::Real`: The opacity of the mesh. Default is 1.

___

### Geometry Manipulation

#### trans!

```julia
trans!(geo::GenericTrace, dis::Vector{<:Real})
```

Translates a 3D geometry by a specified displacement vector.


- `geo::GenericTrace`: The geometry to translate.
- `dis::Vector{<:Real}`: A vector of three Reals specifying the translation distances for the x, y, and z axes.

___

#### rot!

```julia
rot!(geo::GenericTrace, rotang::Vector{<:Real}, center::Vector{<:Real}=[0])
```

Rotates a 3D geometry around a specified center point.


- `geo::GenericTrace`: The 3D geometry to be rotated, which must have `x`, `y`, and `z` coordinates.
- `rotang::Vector{<:Real}`: A vector of three Tait–Bryan rotation angles in degrees for rotations around the x, y, and z axes respectively.
- `center::Vector{<:Real}`: The center point of rotation. Default is `[0]`, which means the rotation center will be set at the geometric center of the object.

___

#### sort_pts

```julia
sort_pts(pts::Vector)
```

Sorts points based on their angular position relative to the centroid.


- `pts::Vector`: List of points to be sorted.

___

#### sort_pts!

```julia
sort_pts!(pts::Vector)
```

Sorts points in place based on their angular position relative to the centroid.


- `pts::Vector`: List of points to be sorted.

___

### Additional Features

#### add_ref_axes

```julia
add_ref_axes(plt::PlotlyJS.SyncPlot, origin::Vector{<:Real}=[0, 0, 0], r::Real=1)
```

Adds reference axes (x, y, z) to a plot.


- `plt::PlotlyJS.SyncPlot`: The plot to which the axes will be added.
- `origin::Vector{<:Real}`: The origin point of the axes.
- `r::Real`: The length of the reference axes.
___

