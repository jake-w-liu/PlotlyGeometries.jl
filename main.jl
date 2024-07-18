using Pkg, Revise
Pkg.activate(".")
using PlotlyJS, PlotlyGeometries

fig = plot()
fig.plot.layout = blank_layout()
display(fig)

add_ref_axes!(fig, [0, 0, 0], [2, 1, 1])
add_arrows!(fig, [0, 0, 0], [1, 0, 0], 1)