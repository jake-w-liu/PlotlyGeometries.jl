# Manual for PlotlyGeometries.jl

## Overview

This manual contains the API provided by `PlotlyGeometries.jl`. Outputs of the functionalities are mostly traces to be added on Plotly plots. In order to use this module, `PlotlyJS.jl` should be installed first. 

## API

### Geometry Creation

#### cuboids

```julia
cuboids(origin::Vector{<:Real}, dimension::Vector{<:Real}, color::String=""; opc::Real=1)
```
Creates a 3D box mesh centered at the given origin with specified dimensions and color.

##### Arguments
- `origin::Vector{<:Real}`: A vector of three Reals specifying the center of the box.
- `dimension::Vector{<:Real}`: A vector of three Reals specifying the dimensions (width, height, depth) of the box.
- `color::String`: A string specifying the color of the box.

##### Keywords
- `opc`: (optional) A Real specifying the opacity of the box. Default is 1.   

#### cubes

```julia
cubes(origin::Vector{<:Real}, side::Real, color::String=""; opc::Real=1)
```

Creates a 3D cube mesh centered at the given origin with specified dimensions and color.

##### Arguments
- `origin::Vector{<:Real}`: A vector of three Reals specifying the center of the cube.
- `side::Real`: Side length of the cube.
- `color::String`: A string specifying the color of the cube.

##### Keywords
- `opc`: (optional) A Real specifying the opacity of the cube. Default is 1. 

___

#### squares

```julia
squares(origin::Vector{<:Real}, side::Real, color::String="", mode::String="z"; opc::Real=1)
```

Creates a 2D square mesh centered at the given origin with the specified side length and color.

##### Arguments
- `origin::Vector{<:Real}`: A vector of three Reals specifying the center of the square.
- `side::Real`: A Real specifying the side length of the square.
- `color::String`: A string specifying the color of the square.
- `mode`::String: (optional) A string specifying the orientation of the square ("x", "y", or "z"). Default is "z".

##### Keywords
- `opc`: (optional) A Real specifying the opacity of the square. Default is 1.

___

#### ellipsoids

```julia
ellipsoids(origin::Vector{<:Real}, par::Vector{<:Real}, color::String=""; opc::Real=1, tres=61, pres=31, ah::Real=0)
```

Creates a 3D ellipsoid mesh.

##### Arguments
- `origin::Vector{<:Real}`: The center of the ellipsoid.
- `par::Vector{<:Real}`: Parameters of the ellipsoid (a, b, c).
- `color::String`: The color of the ellipsoid.

##### Keywords
- `opc`: The opacity of the ellipsoid. Default is 1.
- `tres`: The resolution of the mesh grid (theta). Default is 61.
- `pres`: The resolution of the mesh grid (phi). Default is 31.
- `ah`: alphahull value. Default is 0.
___

#### spheres

```julia
spheres(origin::Vector{<:Real}, r::Real, color::String=""; opc::Real=1, tres=60, pres=30, ah::Real=0)

```

Creates a 3D sphere mesh.

##### Arguments
- `origin::Vector{<:Real}`: The center of the ellipsoid.
- `r::Real`: Radius of the sphere.
- `color::String`: The color of the ellipsoid.

##### Keywords
- `opc`: The opacity of the ellipsoid. Default is 1.
- `tres`: The resolution of the mesh grid (theta). Default is 60.
- `pres`: The resolution of the mesh grid (phi). Default is 30.
- `ah`: alphahull value. Default is 0.

___

#### lines

```julia
lines(pt1::Vector{<:Real}, pt2::Vector{<:Real}, color::String; opc::Real=1, style="")

```

Creates a 3D line between two points.

##### Arguments
- `pt1::Vector{<:Real}`: Starting point of the line.
- `pt2::Vector{<:Real}`: Ending point of the line.
- `color::String`: The color of the line.

##### Keywords
- `opc`: The opacity of the line. Default is 1.
- `style`: The line style (e.g., "solid", "dash"). Default is "".

___

#### polygons

```julia
polygons(pts::Vector, color::String; opc::Real=1, ah::Real=0)
```

Creates a polygon mesh from a set of points (form around the mid point of the set of points).

##### Arguments
- `pts::::Vector`: List of points defining the polygon.
- `color::String`: The color of the polygon.

##### Keywords
- `opc`: The opacity of the polygon. Default is 1.
- `ah`: alphahull value. Default is 0.

```julia
polygons(pts::Vector, ng::Int, color::String; opc::Real=1, ah::Real=0)
```

Creates a group of polygons from a set of points and a specified number of vertices per polygon.

##### Arguments
- `pts::Vector`: List of points defining the mesh.
- `ng::Int`: Number of vertices per polygon.
- `color::String`: The color of the mesh.

##### Keywords
- `opc`: The opacity of the mesh. Default is 1.
- `ah`: alphahull value.

___

### Geometry Manipulation

#### translation

```julia
gtrans!(geo::GenericTrace, dis::Vector{<:Real})
```

Translates a 3D geometry by a specified displacement vector.

##### Arguments
- `geo::GenericTrace`: The geometry to translate.
- `dis::Vector{<:Real}`: A vector of three Reals specifying the translation distances for the x, y, and z axes.

___

#### rotation (according to Tait–Bryan angle)

```julia
grot!(geo::GenericTrace, rotang::Vector{<:Real}, center::Vector{<:Real}=[0])
```

Rotates a 3D geometry around a specified center point. (Tait–Bryan rotation)

##### Arguments
- `geo::GenericTrace`: The 3D geometry to be rotated, which must have `x`, `y`, and `z` coordinates.
- `rotang::Vector{<:Real}`: A vector of three Tait–Bryan rotation angles in degrees  for rotations around the x, y, and z axes respectively.
- `center::Vector{<:Real}`: The center point of rotation. Default is `[0]`, which means the rotation center will be set at the geometric center of the object.

___

#### rotation (according to axis)

```julia
grot!(geo::GenericTrace, ang::Real, axis::Vector{<:Real}, origin::Vector{<:Real}=[0])
```

Rotates the geometry by the specified angle `ang` around the axis `axis` and origin `origin`.

##### Arguments
- `geo::GenericTrace`: The geometry to be rotated.
- `ang::Real`: The rotation angle.
- `axis::Vector{<:Real}`: The rotation axis.
- `origin::Vector{<:Real}=[0]`: The rotation origin. Defaults to the center of the geometry if not specified.

___

### Additional Features

#### sort points

```julia
sort_pts(pts::Vector)
```

Sorts points based on their angular position relative to the centroid.

##### Arguments
- `pts::Vector`: List of points to be sorted.

___

#### inplace sort points

```julia
sort_pts!(pts::Vector)
```

Sorts points in place based on their angular position relative to the centroid.

##### Arguments
- `pts::Vector`: List of points to be sorted.

___

#### add reference axis

```julia
add_ref_axes!(plt::PlotlyJS.SyncPlot, origin::Vector{<:Real}, r::Real)
```

Adds reference axes (x, y, z) to a plot.

##### Arguments
- `plt::PlotlyJS.SyncPlot`: The plot to which the axes will be added.
- `origin::Vector{<:Real}`: The origin point of the axes.
- `r::Real`: The length of the reference axes.
___

```julia
add_ref_axes!(plt::PlotlyJS.SyncPlot, origin::Vector{<:Real}, r::Vector{<:Real}=[1, 1, 1])
```

Adds reference axes (x, y, z) to a plot.

##### Arguments
- `plt::PlotlyJS.SyncPlot`: The plot to which the axes will be added.
- `origin::Vector{<:Real}`: The origin point of the axes.
- `r::Vector{<:Real}`: The lengths of the reference axes.
___

#### add arrows

```julia
add_arrows!(plt::PlotlyJS.SyncPlot, origin::Vector{<:Real}, dir::Vector{<:Real}, len::Float64=1.0, color::String=""; opc::Real=1)

```

Creates a 3D arrow starting from a point and pointing in a given direction.

##### Arguments
- `plt::PlotlyJS.SyncPlot`: The plot to which the axes will be added.
- `origin::Vector{<:Real}`: The starting point of the arrow.
- `dir::Vector{<:Real}`: The direction vector of the arrow.
- `len::Real`: length of the arrow
- `color::String`: The color of the arrow.

##### Keywords
- `opc`: The opacity of the arrow. Default is 1.

___

#### add texts

```julia
add_text!(plt::PlotlyJS.SyncPlot, origin::Vector{<:Real}, text::String, color::String="")
```

Add text to plot. 

##### Arguments
- plt::PlotlyJS.SyncPlot: Plot to add text.
- origin::Vector{<:Real}: origin of the text.
- text::String:: text to be added.
- color::String=: color of the text.

___

#### blank layout

```julia
blank_layout()
```

Return blank layout for easy use.
___
