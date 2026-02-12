# Manual for PlotlyGeometries.jl

## Overview

`PlotlyGeometries.jl` builds and manipulates 3D geometry traces for Plotly figures, using `PlotlySupply.jl`.

Most constructors return a `mesh3d` or `scatter3d` trace (a `GenericTrace`), and mutating utilities operate in-place on those traces.

For plotting and display, use `PlotlySupply`:

```julia
using PlotlySupply
using PlotlyGeometries

geo = cuboids([0, 0, 0], [1, 2, 3], "steelblue"; opc=0.4)
fig = plot(geo, blank_layout())
add_ref_axes!(fig)
set_view!(fig, 35, 20)
display(fig)
```

## API Summary

### Geometry Constructors

- `cuboids`
- `cubes`
- `squares`
- `ellipsoids`
- `spheres`
- `cylinders`
- `cones` (cone + frustum overloads)
- `disks`
- `planes`
- `tori`
- `lines`
- `polygons` (single polygon + grouped polygons)

### Geometry Transforms

- `gtrans!` (single trace or vector of traces)
- `grot!` (Tait-Bryan angles; or axis-angle; both support single trace or vector of traces)

### Utilities

- `sort_pts`
- `sort_pts!`

### Plot Helpers

- `add_ref_axes!` (scalar or per-axis lengths)
- `add_arrows!`
- `add_text!`
- `blank_layout`
- `set_view!`

## Conventions

- Coordinates are 3D vectors `[x, y, z]`.
- `axis` options are `"x"`, `"y"`, `"z"`.
- If `color == ""`, a random RGB color is generated.
- Angles are in degrees for rotation/view APIs.
- Plot helper functions accept `Union{Plot, SyncPlot}` from `PlotlySupply`.

---

## Geometry Constructors

### `cuboids`

```julia
cuboids(origin::Vector{<:Real}, dimension::Vector{<:Real}, color::String=""; opc::Real=1)
```

Creates a box mesh centered at `origin`.

Arguments:
- `origin`: `[x, y, z]`
- `dimension`: `[dx, dy, dz]`
- `color`

Keywords:
- `opc`: opacity

Returns:
- `mesh3d` trace

Validation:
- `length(origin) == 3`
- `length(dimension) == 3`

### `cubes`

```julia
cubes(origin::Vector{<:Real}, side::Real, color::String=""; opc::Real=1)
```

Convenience wrapper for `cuboids` with equal side lengths.

Validation:
- `length(origin) == 3`
- `side > 0`

### `squares`

```julia
squares(origin::Vector{<:Real}, side::Real, mode::String="z", color::String=""; opc::Real=1)
```

Creates a square mesh oriented by `mode`:
- `"z"`: square lies in XY plane
- `"x"`: square lies in YZ plane
- `"y"`: square lies in XZ plane

Validation:
- `length(origin) == 3`

### `ellipsoids`

```julia
ellipsoids(origin::Vector{<:Real}, par::Vector{<:Real}, color::String=""; opc::Real=1, tres=61, pres=31, ah::Real=0)
```

Creates an ellipsoid parameterized by semi-axes `par = [a, b, c]`.

Keywords:
- `opc`: opacity
- `tres`: theta resolution
- `pres`: phi resolution
- `ah`: `alphahull` passed to Plotly mesh

Validation:
- `length(origin) == 3`
- `length(par) == 3`
- `all(par .> 0)`
- `tres > 0`
- `pres > 0`

### `spheres`

```julia
spheres(origin::Vector{<:Real}, r::Real, color::String=""; opc::Real=1, tres=60, pres=30, ah::Real=0)
```

Sphere wrapper around `ellipsoids(origin, [r, r, r], ...)`.

Validation:
- `length(origin) == 3`
- `r > 0`

### `cylinders`

```julia
cylinders(origin::Vector{<:Real}, r::Real, h::Real, axis::String="z", color::String=""; opc::Real=1, tres::Int=60, caps::Bool=true)
```

Creates a cylinder centered at `origin`.

Arguments:
- `r`: radius
- `h`: height
- `axis`: axis of the cylinder (`"x"`, `"y"`, `"z"`)

Keywords:
- `opc`
- `tres`: circumferential resolution
- `caps`: whether to close top and bottom

Validation:
- `length(origin) == 3`
- `r > 0`
- `h > 0`
- `tres >= 3`
- `axis in ("x", "y", "z")`

### `cones` (cone overload)

```julia
cones(origin::Vector{<:Real}, r::Real, h::Real, axis::String="z", color::String=""; opc::Real=1, tres::Int=60, caps::Bool=true)
```

Creates a cone (`r` at base, `0` at top).

### `cones` (frustum overload)

```julia
cones(origin::Vector{<:Real}, r1::Real, r2::Real, h::Real, axis::String="z", color::String=""; opc::Real=1, tres::Int=60, caps::Bool=true)
```

Creates a frustum/cone with base radius `r1` and top radius `r2`.

Validation:
- `length(origin) == 3`
- `r1 >= 0`
- `r2 >= 0`
- `r1 + r2 > 0`
- `h > 0`
- `tres >= 3`
- `axis in ("x", "y", "z")`

### `disks`

```julia
disks(origin::Vector{<:Real}, r::Real, axis::String="z", color::String=""; opc::Real=1, tres::Int=60)
```

Creates a circular disk.

Validation:
- `length(origin) == 3`
- `r > 0`
- `tres >= 3`
- `axis in ("x", "y", "z")`

### `planes`

```julia
planes(origin::Vector{<:Real}, dims::Vector{<:Real}, axis::String="z", color::String=""; opc::Real=1)
```

Creates a rectangular plane centered at `origin`, with dimensions `dims = [d1, d2]`, normal to `axis`.

Validation:
- `length(origin) == 3`
- `length(dims) == 2`
- `all(dims .> 0)`
- `axis in ("x", "y", "z")`

### `tori`

```julia
tori(origin::Vector{<:Real}, R::Real, r::Real, axis::String="z", color::String=""; opc::Real=1, ures::Int=61, vres::Int=31)
```

Creates a torus with:
- `R`: major radius
- `r`: minor radius

Keywords:
- `ures`: major-circle resolution
- `vres`: minor-circle resolution

Validation:
- `length(origin) == 3`
- `R > 0`
- `r > 0`
- `ures >= 3`
- `vres >= 3`
- `axis in ("x", "y", "z")`

### `lines`

```julia
lines(pt1::Vector{<:Real}, pt2::Vector{<:Real}, color::String=""; opc::Real=1, style="")
```

Creates a `scatter3d` line trace.

Keywords:
- `opc`
- `style`: line dash style

Validation:
- `length(pt1) == 3`
- `length(pt2) == 3`

### `polygons` (single polygon)

```julia
polygons(pts::Vector, color::String=""; opc::Real=1, ah::Real=0)
```

Creates a polygon mesh from coplanar points.

Implementation detail:
- Uses ear-clipping triangulation; supports concave polygons.
- Expects a simple, non-self-intersecting boundary.

Validation:
- each point has length 3 and real entries
- `length(pts) >= 3`

### `polygons` (multiple polygons with fixed group size)

```julia
polygons(pts::Vector, ng::Int, color::String=""; opc::Real=1, ah::Real=0)
```

Interprets `pts` as consecutive groups of `ng` points and triangulates each group independently.

Validation:
- each point has length 3 and real entries
- `ng > 0`
- `length(pts) >= ng`
- `length(pts) % ng == 0`

---

## Geometry Transforms

### `gtrans!` (single trace)

```julia
gtrans!(geo::GenericTrace, dis::Vector{<:Real})
```

Translates `geo` in-place by displacement `dis`.

Validation:
- `length(dis) == 3`

### `gtrans!` (vector of traces)

```julia
gtrans!(geos::AbstractVector{<:GenericTrace}, dis::Vector{<:Real})
```

Applies translation to each trace.

### `grot!` (Tait-Bryan, single trace)

```julia
grot!(geo::GenericTrace, rotang::Vector{<:Real}, center::Vector{<:Real}=[0])
```

In-place rotation with Tait-Bryan angles `[alpha, beta, gamma]` in degrees.

Behavior:
- If `center == [0]`, uses the geometry centroid.

Validation:
- `length(rotang) == 3`
- if center specified, `length(center) == 3`

### `grot!` (Tait-Bryan, vector of traces)

```julia
grot!(geos::AbstractVector{<:GenericTrace}, rotang::Vector{<:Real}, center::Vector{<:Real}=[0])
```

Applies same rotation to each trace.

### `grot!` (axis-angle, single trace)

```julia
grot!(geo::GenericTrace, ang::Real, axis::Vector{<:Real}, origin::Vector{<:Real}=[0])
```

In-place rotation around `axis` by angle `ang` (degrees), using Rodrigues' formula.

Behavior:
- If `origin == [0]`, uses geometry centroid.

Validation:
- `length(axis) == 3`
- `norm(axis) > 0`
- if origin specified, `length(origin) == 3`

### `grot!` (axis-angle, vector of traces)

```julia
grot!(geos::AbstractVector{<:GenericTrace}, ang::Real, axis::Vector{<:Real}, origin::Vector{<:Real}=[0])
```

Applies same axis-angle rotation to each trace.

---

## Utilities

### `sort_pts`

```julia
sort_pts(pts::Vector)
```

Returns a sorted copy of points by angular order around centroid.

Validation:
- all points are 3D

### `sort_pts!`

```julia
sort_pts!(pts::Vector)
```

Sorts points in-place by angular order around centroid.

Validation:
- all points are 3D

---

## Plot Helpers

All helpers below operate on:

```julia
plt::Union{Plot, SyncPlot}
```

### `add_ref_axes!` (uniform length)

```julia
add_ref_axes!(plt::Union{Plot, SyncPlot}, origin::Vector{<:Real}=[0, 0, 0], r::Real=1)
```

Adds x/y/z colored axes with cone tips and annotation labels.

Validation:
- `length(origin) == 3`
- `r > 0`

### `add_ref_axes!` (per-axis lengths)

```julia
add_ref_axes!(plt::Union{Plot, SyncPlot}, origin::Vector{<:Real}, r::Vector{<:Real})
```

Same as above but with per-axis lengths `[rx, ry, rz]`.

Validation:
- `length(origin) == 3`
- `length(r) == 3`
- `all(r .> 0)`

### `add_arrows!`

```julia
add_arrows!(plt::Union{Plot, SyncPlot}, origin::Vector{<:Real}, dir::Vector{<:Real}, len::Real=1.0, color::String=""; opc::Real=1, endpoint::Bool=true, asize::Real=len)
```

Adds a cone+line arrow.

Behavior:
- `endpoint=true`: arrow starts at `origin` and points toward endpoint.
- `endpoint=false`: `origin` is arrow center.
- line width scales with `asize/len`.

Validation:
- `length(origin) == 3`
- `length(dir) == 3`
- `norm(dir) > 0`
- `len > 0`
- `asize > 0`

### `add_text!`

```julia
add_text!(plt::Union{Plot, SyncPlot}, origin::Vector{<:Real}, text::String, color::String="")
```

Adds a 3D text annotation as a `scatter3d` text trace.

### `blank_layout`

```julia
blank_layout()
```

Returns a layout with hidden axes/grid and `scene.aspectmode = "data"`.

### `set_view!`

```julia
set_view!(plt::Union{Plot, SyncPlot}, az::Real, el::Real, r::Real=1.25 * sqrt(3))
```

Sets camera eye from azimuth/elevation/radius.

Behavior:
- `el == ±90` is nudged slightly to avoid singularity.
- Uses `relayout!(..., scene_camera=...)`.

---

## Notes on Assertions and Errors

This package primarily uses `@assert` for argument validation. Invalid inputs therefore raise `AssertionError`.

Common invalid-input cases:
- wrong vector sizes (not 3D where required)
- non-positive radii/heights/resolutions
- zero direction vector in `add_arrows!`
- zero axis vector in axis-angle `grot!`
- `polygons(pts, ng)` where `length(pts)` is not a multiple of `ng`

---

## Examples

See `examples/`:
- `examples/ex_basics.jl`: comprehensive geometry showcase
- `examples/ex_polygons.jl`: polygon and grouped polygon usage
