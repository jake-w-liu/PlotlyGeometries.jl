push!(LOAD_PATH, "../src")
using PlotlyJS
using PlotlyGeometries

b1 = boxes([0, 0, 0], [1, 2, 3], "pink", 0.2)

rotate!(b1, [10, 20, 30])

fig = plot(b1, blank_layout())
add_ref_axes(fig, [0, 0, 0], 1)

display(fig)

# e1 = ellipsoids([0, 5, 0], [3, 2, 1], "purple", [0, 0, 0], 0.3)
# addtraces!(plt, e1)
