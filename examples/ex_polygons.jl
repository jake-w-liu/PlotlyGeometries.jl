using PlotlyJS
using PlotlyGeometries

# create two triangles 
pts = []
push!(pts, [0, 0, 0])
push!(pts, [0, 1, 0])
push!(pts, [1, 0, 0])

p1 = polygons(pts, "aqua"; opc=0.5)

for p in pts
    p[3] += 1
end

p2 = polygons(pts, "aqua"; opc=0.5)

fig = plot([p1, p2], blank_layout())
display(fig)

add_ref_axes!(fig)

# create a set of rectangles
pts = []

# first side
push!(pts, [0, 0, 0])
push!(pts, [1, 0, 0])
push!(pts, [1, 0, 1])
push!(pts, [0, 0, 1])

# second side
push!(pts, [0, 0, 0])
push!(pts, [0, 1, 0])
push!(pts, [0, 1, 1])
push!(pts, [0, 0, 1])

# third side
push!(pts, [0, 1, 0])
push!(pts, [1, 0, 0])
push!(pts, [0, 1, 1])
push!(pts, [1, 0, 1])

p3 = polygons(pts, 4, "yellow"; opc=0.5)
addtraces!(fig, p3)

for n = 0:360
    set_view!(fig, n, 90)
    sleep(0.1)
end

