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
e1 = ellipsoids([0, 5, 0], [3, 1, 1], "purple", 0.1)

# rotate ellipsoid (Tait–Bryan angles)
rot!(e1, [10, 20, 30])

# add elipsoid to figure
addtraces!(fig, e1)

# create sphere
s1 = spheres([3, -2, 0], 2, "coral", 0.1)
addtraces!(fig, s1)

# add line
l1 = lines([2, 4, 6], [0, 5, 0], "orange", 1, "dash")
addtraces!(fig, l1)

# add arrow and text 
add_arrows(fig, [2, 4, 8], [0, 0, -1], "black", 1/4)
add_text(fig, [2, 4, 8.5], "this is a cube", "black")

add_arrows(fig, [0, 6, 2], [0, -1, -1], "black", 1/4)
add_text(fig, [0, 6, 2.5], "this is an ellipsoid", "red")






