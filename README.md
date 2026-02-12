# PlotlyGeometries.jl

[![Build Status](https://github.com/akjake616/PlotlyGeometries.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/akjake616/PlotlyGeometries.jl/actions/workflows/CI.yml)

`PlotlyGeometries.jl` is a Julia package designed for creating and manipulating 3D geometrical shapes and visualizations using [`PlotlySupply.jl`](https://github.com/jake-w-liu/PlotlySupply.jl). This package provides a variety of functions to easily generate and customize 3D shapes such as boxes, spheres, ellipsoids, cylinders, cones/frustums, tori, planes, disks, lines, and arrows, as well as utility functions for transformations and visual enhancements.I hope this package will be useful for those trying to create better illustrations for their academic research using Julia :wink:


<p align="center">
  <img alt="PlotlyGeometries.jl" src="./media/illus.gif" width="50%" height="auto" />
</p>

## Installation

To install `PlotlyGeometries.jl`, use the following command in the Julia REPL:

```julia
using Pkg
Pkg.add("PlotlyGeometries")
```

## Learn by Examples

Please refer to the examples to get familiar with some basic usages:

- `ex_basics.jl`: Full showcase of classic + newly added geometries.
- `ex_polygons.jl`: Use polygons to build geometries.

## Highlights

- Concave polygon triangulation using ear clipping.
- Batch transforms via `gtrans!`/`grot!` on vectors of traces.
- Plot helpers for reference axes, arrows, text annotations, and camera view.

## Usage

Please refer to the [user manual](./docs/MANUAL.md) in the docs folder.
