push!(LOAD_PATH, "../src")
using PlotlyJS
using PlotlyGeometries

# create a simple cube
c1 = cubes([0, 0, 0], [1, 2, 3], "pink", 0.2)

# translate the cube
trans!(c1, [2, 4, 6])

# create fig with PlotlyJS (use blank_layout() to easily create a blank layout)
fig = plot(c1, blank_layout())

# add reference axes
add_ref_axes(fig, [0, 0, 0], 1)

# show figure
display(fig)

# create ellipsoid
e1 = ellipsoids([0, 5, 0], [3, 1, 1], "purple", 0.2)

# rotate ellipsoid (Tait–Bryan angles)
rot!(e1, [10, 20, 30])

# add elipsoid to figure
add_trace!(fig, e1)

# add arrow and text 
add_arrows(fig, [2, 4, 8], [0, 0, -1], "black", 1/4)
add_text(fig, [2, 4, 8.5], "this is a cube", "black")

add_arrows(fig, [0, 6, 2], [0, -1, -1], "black", 1/4)
add_text(fig, [0, 6, 2.5], "this is an ellipsoid", "red")






