using PlotlySupply
using PlotlyGeometries

# Classic geometries
c1 = cuboids([0, 0, 0], [1, 2, 3], "deepskyblue"; opc=0.25)
e1 = ellipsoids([0, 5, 0], [3, 1, 1], "purple"; opc=0.18)
s1 = spheres([3, -2, 0], 2, "coral"; opc=0.18)
l1 = lines([2, 4, 6], [0, 5, 0], "orange"; style="dash")

# Transform classic geometries
gtrans!(c1, [2, 4, 6])
grot!(c1, 45, [0, 0, 1])
grot!(e1, [10, 20, 30])

# New primitives
cy1 = cylinders([-4.0, 0.0, 0.0], 0.8, 2.5, "z", "steelblue"; opc=0.45, tres=50)
co1 = cones([-1.5, 0.0, 0.0], 0.9, 2.0, "z", "tomato"; opc=0.45, tres=50)
fr1 = cones([1.5, 0.0, 0.0], 1.0, 0.45, 2.0, "z", "goldenrod"; opc=0.45, tres=50)
t1 = tori([4.0, 0.0, 0.0], 1.3, 0.45, "z", "teal"; opc=0.35, ures=70, vres=40)
pl1 = planes([0.0, -3.2, 0.0], [2.2, 2.2], "y", "lightgray"; opc=0.40)
dk1 = disks([0.0, 3.0, 0.0], 1.2, "x", "purple"; opc=0.50, tres=60)

# Group transform support on multiple traces
grot!([cy1, co1, fr1, t1, pl1, dk1], [20, 0, 25])

# Create figure with all geometries
fig = plot([c1, e1, s1, l1, cy1, co1, fr1, t1, pl1, dk1], blank_layout())
add_ref_axes!(fig, [0, 0, 0], [1.5, 1.5, 1.5])

# Labels and arrows
add_arrows!(fig, [2, 4, 8], [0, 0, -1], 0.5, "black")
add_text!(fig, [2, 4, 8.5], "this is a cube")

add_arrows!(fig, [0, 6, 2], [0, -1, -1], 0.5, "black")
add_text!(fig, [0, 6, 2.5], "this is an ellipsoid")

add_arrows!(fig, [3, -2, 2], [0, 0, -1], 0.5, "black")
add_text!(fig, [3, -2, 2.5], "this is a sphere")

add_text!(fig, [-4.1, 0.0, 1.8], "cylinder", "steelblue")
add_text!(fig, [-1.6, 0.0, 1.7], "cone", "tomato")
add_text!(fig, [1.5, 0.0, 1.7], "frustum", "goldenrod")
add_text!(fig, [4.0, 0.0, 1.6], "torus", "teal")
add_text!(fig, [0.0, -3.0, 1.4], "plane", "black")
add_text!(fig, [0.0, 3.0, 1.4], "disk", "purple")

set_view!(fig, 35, 20)
display(fig)
