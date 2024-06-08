# using PlotlyJS
# using Combinatorics
# using LinearAlgebra

"""
    boxes(origin::Vector{<:Number}, dimension::Vector{<:Number}, color::String, opc=1)

Creates a 3D box mesh centered at the given origin with specified dimensions and color.

# Arguments
- `origin::Vector{<:Number}`: A vector of three numbers specifying the center of the box.
- `dimension::Vector{<:Number}`: A vector of three numbers specifying the dimensions (width, height, depth) of the box.
- `color::String`: A string specifying the color of the box.
- `opc`: (optional) A number specifying the opacity of the box. Default is 1.                                   
"""
function boxes(origin::Vector{<:Number}, dimension::Vector{<:Number}, color::String, opc=1)
    @assert length(origin) == 3
    @assert length(dimension) == 3

    x1 = origin[1] - dimension[1] / 2
    x2 = origin[1] + dimension[1] / 2
    y1 = origin[2] - dimension[2] / 2
    y2 = origin[2] + dimension[2] / 2
    z1 = origin[3] - dimension[3] / 2
    z2 = origin[3] + dimension[3] / 2

    x = [x1, x1, x2, x2, x1, x1, x2, x2]
    y = [y1, y2, y2, y1, y1, y2, y2, y1]
    z = [z1, z1, z1, z1, z2, z2, z2, z2]
    i = [7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2]
    j = [3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3]
    k = [0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6]

    return mesh3d(x=x, y=y, z=z,
        i=i, j=j, k=k,
        alphahull=1,
        flatshading=true,
        color=color,
        opacity=opc,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0,
        ),
    )
end

"""
squares(origin::Vector{<:Number}, side::Number, color::String, mode="z", opc=1)

Creates a 2D square mesh centered at the given origin with the specified side length and color.

# Arguments
- `origin::Vector{<:Number}`: A vector of three numbers specifying the center of the square.
- `side::Number`: A number specifying the side length of the square.
- `color::String`: A string specifying the color of the square.
- `mode`::String: (optional) A string specifying the orientation of the square ("x", "y", or "z"). Default is "z".
- `opc`: (optional) A number specifying the opacity of the square. Default is 1.
"""
function squares(origin::Vector{<:Number}, side::Number, color::String, mode="z", opc=1)
    @assert length(origin) == 3

    if mode == "x"
        x = [origin[1], origin[1], origin[1], origin[1]]
        y = [origin[2] - side / 2, origin[2] + side / 2, origin[2] + side / 2, origin[2] - side / 2]
        z = [origin[3] - side / 2, origin[3] - side / 2, origin[3] + side / 2, origin[3] + side / 2]
        i = [0, 2]
        j = [1, 3]
        k = [2, 0]
    elseif mode == "y"
        x = [origin[1] - side / 2, origin[1] - side / 2, origin[1] + side / 2, origin[1] + side / 2]
        y = [origin[2], origin[2], origin[2], origin[2]]
        z = [origin[3] - side / 2, origin[3] + side / 2, origin[3] + side / 2, origin[3] - side / 2]
        i = [0, 2]
        j = [1, 3]
        k = [2, 0]
    else
        x = [origin[1] - side / 2, origin[1] - side / 2, origin[1] + side / 2, origin[1] + side / 2]
        y = [origin[2] - side / 2, origin[2] + side / 2, origin[2] + side / 2, origin[2] - side / 2]
        z = [origin[3], origin[3], origin[3], origin[3]]
        i = [0, 2]
        j = [1, 3]
        k = [2, 0]
    end

    return mesh3d(x=x, y=y, z=z, i=i, j=j, k=k,
        color=color,
        opacity=opc,
        alphahull=1,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0
        ),
    )
end

"""
    polygons(pts::Vector{Vector{<:Number}}, color::String, opc=1)

Creates a polygon mesh from a set of points.

# Arguments
- `pts::Vector{Vector{<:Number}}`: List of points defining the polygon.
- `color::String`: The color of the polygon.
- `opc`: The opacity of the polygon. Default is 1.
"""
function polygons(pts::Vector{Vector{<:Number}}, color::String, opc=1)
    @assert all(length.(pts) .== 3)

    N = length(pts)
    sort_pts!(pts)
    x = []
    y = []
    z = []
    for i in eachindex(pts)
        push!(x, pts[i][1])
        push!(y, pts[i][2])
        push!(z, pts[i][3])
    end

    mid = zeros(3)
    for n in eachindex(pts)
        for m in 1:3
            mid[m] += pts[n][m]
        end
    end
    mid = mid ./ N
    push!(x, mid[1])
    push!(y, mid[2])
    push!(z, mid[3])

    a = 0:1:length(pts)
    i = []
    j = []
    k = []

    for n = 0:length(a)-2
        push!(i, a[mod(n + 0, length(a) - 1)+1])
        push!(j, a[mod(n + 1, length(a) - 1)+1])
        push!(k, a[end])
    end

    return mesh3d(x=x, y=y, z=z,
        i=i, j=j, k=k,
        color=color,
        opacity=opc,
        alphahull=1,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0
        ),
    )
end

"""
    ellipsoids(origin::Vector{<:Number}, par::Vector{<:Number}, color::String, rotang::Vector{<:Number}=[0, 0, 0], opc=1, res=25)

Creates a 3D ellipsoid mesh.

# Arguments
- `origin::Vector{<:Number}`: The center of the ellipsoid.
- `par::Vector{<:Number}`: Parameters of the ellipsoid (a, b, c).
- `color::String`: The color of the ellipsoid.
- `rotang::Vector{<:Number}`: Rotation angles in degrees (alpha, beta, gamma). Default is [0, 0, 0].
- `opc`: The opacity of the ellipsoid. Default is 1.
- `res`: The resolution of the mesh grid. Default is 25.
"""
function ellipsoids(origin::Vector{<:Number}, par::Vector{<:Number}, color::String, opc=1, res=25)
    @assert length(origin) == 3
    @assert length(par) == 3

    a = par[1]
    b = par[2]
    c = par[3]



    P, T = meshgrid(
        LinRange(0, 2 * pi, res),
        LinRange(0, pi, res),
    )

    x = sin.(T) .* cos.(P) .* a
    y = sin.(T) .* sin.(P) .* b
    z = cos.(T) .* c
    x = x[:]
    y = y[:]
    z = z[:]



    return mesh3d(x=x, y=y, z=z,
        alphahull=0,
        flatshading=true,
        color=color,
        opacity=opc,
        lighting=attr(
            diffuse=0.1,
            specular=2.0,
            roughness=0.5
        ),
    )
end

"""
    lines(pt1::Vector{<:Number}, pt2::Vector{<:Number}, color::String, opc=1, style="")

Creates a 3D line between two points.

# Arguments
- `pt1::Vector{<:Number}`: Starting point of the line.
- `pt2::Vector{<:Number}`: Ending point of the line.
- `color::String`: The color of the line.
- `opc`: The opacity of the line. Default is 1.
- `style`: The line style (e.g., "solid", "dash"). Default is "".
"""
function lines(pt1::Vector{<:Number}, pt2::Vector{<:Number}, color::String, opc=1, style="")

    x = [pt1[1], pt2[1]]
    y = [pt1[2], pt2[2]]
    z = [pt1[3], pt2[3]]

    return scatter3d(x=x, y=y, z=z,
        mode="lines",
        line=attr(
            color=color,
            dash=style,
        ),
        showlegend=false,
        alphahull=1,
        flatshading=true,
        opacity=opc,
    )
end

"""
arrows(origin::Vector{<:Number}, direction::Vector{<:Number}, color::String, opc=1)

Creates a 3D arrow starting from a point and pointing in a given direction.

# Arguments
- `origin::Vector{<:Number}`: The starting point of the arrow.
- `direction::Vector{<:Number}`: The direction vector of the arrow.
- `color::String`: The color of the arrow.
- `opc`: The opacity of the arrow. Default is 1.
"""
function arrows(origin::Vector{<:Number}, direction::Vector{<:Number}, color::String, opc=1)
    @assert length(origin) == 3
    @assert length(direction) == 3

    c = cone(x=[origin[1] + direction[1] / 2], y=[origin[2] + direction[2] / 2], z=[origin[3] + direction[3] / 2], u=[direction[1]], v=[direction[2]], w=[direction[3]],
        colorscale=[[0, color], [1, color]],
        opacity=opc,
        showscale=false)

    l = scatter3d(x=[origin[1] - direction[1] / 2, origin[1] + direction[1] / 2], y=[origin[2] - direction[2] / 2, origin[2] + direction[2] / 2], z=[origin[3] - direction[3] / 2, origin[3] + direction[3] / 2],
        line=attr(color=color, width=4),
        mode="lines",
        opacity=opc,
        showlegend=false)

    return [c, l]
end


"""
    create_mesh(pts::Vector{Vector{<:Number}}, ng::Int, color::String, opc=1)

Creates a mesh from a set of points and a specified number of vertices per polygon.

# Arguments
- `pts::Vector{Vector{<:Number}}`: List of points defining the mesh.
- `ng::Int`: Number of vertices per polygon.
- `color::String`: The color of the mesh.
- `opc`: The opacity of the mesh. Default is 1.
"""
function create_mesh(pts::Vector{Vector{<:Number}}, ng::Int, color::String, opc=1)
    @assert all(length.(pts) .== 3)

    N = length(pts)

    x = []
    y = []
    z = []
    i = []
    j = []
    k = []
    Ng = floor(Int, N / ng)
    for p = 1:Ng
        ptsg = []
        for m = 1:ng
            push!(ptsg, pts[(p-1)*ng+m])
        end
        sort_pts!(ptsg)
        for i in eachindex(ptsg)
            push!(x, ptsg[i][1])
            push!(y, ptsg[i][2])
            push!(z, ptsg[i][3])
        end

        mid = zeros(3)
        for n in eachindex(ptsg)
            for m in 1:3
                mid[m] += ptsg[n][m]
            end
        end
        mid = mid ./ ng
        push!(x, mid[1])
        push!(y, mid[2])
        push!(z, mid[3])

        a = 0:1:length(ptsg)
        for n = 0:length(a)-2
            push!(i, a[mod(n + 0, length(a) - 1)+1] + (p - 1) * (ng + 1))
            push!(j, a[mod(n + 1, length(a) - 1)+1] + (p - 1) * (ng + 1))
            push!(k, a[end] + (p - 1) * (ng + 1))
        end
    end

    return mesh3d(x=x, y=y, z=z,
        i=i, j=j, k=k,
        color=color,
        opacity=opc,
        alphahull=1,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0
        ),
    )
end

"""
    rotate!(geo::GenericTrace, rotang::Vector{<:Number}, center::Vector{<:Number}=[0])

Rotates a 3D geometry around a specified center point.

# Arguments
- `geo::GenericTrace`: The 3D geometry to be rotated, which must have `x`, `y`, and `z` coordinates.
- `rotang::Vector{<:Number}`: A vector of three Tait–Bryan rotation angles in degrees `[alpha, beta, gamma]` for rotations around the x, y, and z axes respectively.
- `center::Vector{<:Number}`: The center point of rotation. Default is `[0]`, which means the rotation center will be set at the geometric center of the object.

"""
function rotate!(geo::GenericTrace, rotang::Vector{<:Number}, center::Vector{<:Number}=[0])
    @assert length(rotang) == 3

    alpha = rotang[1]
    beta = rotang[2]
    gama = rotang[3]

    Rx = [1 0 0;
        0 cosd(alpha) -sind(alpha);
        0 sind(alpha) cosd(alpha)]
    Ry = [cosd(beta) 0 sind(beta);
        0 1 0;
        -sind(beta) 0 cosd(beta)]
    Rz = [cosd(gama) -sind(gama) 0;
        sind(gama) cosd(gama) 0;
        0 0 1]

    R = Rz * Ry * Rx

    pos = []
    for n in eachindex(geo.x)
        push!(pos, [geo.x[n], geo.y[n], geo.z[n]])
    end

    if center == [0] # rotation center set at the geometry center
        center = sum(pos) ./ length(pos)
    else
        @assert length(center) == 3
    end

    @inbounds for n in eachindex(geo.x)
        vec = R * (pos[n] - center) 
        geo.x[n] = vec[1] + center[1]
        geo.y[n] = vec[2] + center[2]
        geo.z[n] = vec[3] + center[3]
    end
end

"""
    sort_pts!(pts::Vector{Vector{<:Number}})

Sorts points in place based on their angular position relative to the centroid.

# Arguments
- `pts::Vector{Vector{<:Number}}`: List of points to be sorted.
"""
function sort_pts!(pts::Vector{Vector{<:Number}})
    @assert all(length.(pts) .== 3)

    N = length(pts)
    ang = zeros(length(pts))
    mid = zeros(3)
    for n in eachindex(pts)
        for m in 1:3
            mid[m] += pts[n][m]
        end
    end
    mid = mid ./ N

    c = collect(combinations(1:N, 3))
    thtr = 0
    phir = 0
    for n in eachindex(c)
        vec = cross(pts[c[n][2]] .- pts[c[n][1]], pts[c[n][3]] .- pts[c[n][1]])
        if norm(vec) == 0
            continue
        else
            vec = vec ./ norm(vec)
            thtr = acosd(vec[3])
            phir = atand(vec[2], vec[1])
            break
        end
    end

    Ry = [
        cosd(thtr) 0 -sind(thtr);
        0 1 0;
        sind(thtr) 0 cosd(thtr)
    ]
    Rz = [
        cosd(phir) sind(phir) 0;
        -sind(phir) cosd(phir) 0;
        0 0 1
    ]
    R = Ry * Rz
    pts_rot = similar(pts)
    mid_rot = R * mid
    for n in eachindex(pts)
        pts_rot[n] = R * pts[n]
        ang[n] = atan(pts_rot[n][2] - mid_rot[2], pts_rot[n][1] - mid_rot[1])
    end
    pts .= pts[sortperm(ang)]
end



"""
    add_ref_axes(plt, origin::Vector{<:Number}, r::Number)

Adds reference axes (x, y, z) to a plot.

# Arguments
- `plt`: The plot to which the axes will be added.
- `origin::Vector{<:Number}`: The origin point of the axes.
- `r::Number`: The length of the reference axes.
"""
function add_ref_axes(plt, origin, r)
    cx = cone(x=[r + origin[1]], y=[origin[2]], z=[origin[3]], u=[r / 10], v=[0], w=[0],
        colorscale=[[0, "red"], [1, "red"]],
        showscale=false)
    cy = cone(x=[origin[1]], y=[r + origin[2]], z=[origin[3]], u=[0], v=[r / 10], w=[0],
        colorscale=[[0, "green"], [1, "green"]],
        showscale=false)
    cz = cone(x=[origin[1]], y=[origin[2]], z=[r + origin[3]], u=[0], v=[0], w=[r / 10],
        colorscale=[[0, "blue"], [1, "blue"]],
        showscale=false)
    lx = scatter3d(x=[origin[1], r + origin[1]], y=[origin[2], origin[2]], z=[origin[3], origin[3]],
        line=attr(color="red"),
        mode="lines",
        showlegend=false)
    ly = scatter3d(x=[origin[1], origin[1]], y=[origin[2], r + origin[2]], z=[origin[3], origin[3]],
        line=attr(color="green"),
        mode="lines",
        showlegend=false)
    lz = scatter3d(x=[origin[1], origin[1]], y=[origin[2], origin[2]], z=[origin[3], r + origin[3]],
        line=attr(color="blue"),
        mode="lines",
        showlegend=false)
    addtraces!(plt, cx)
    addtraces!(plt, cy)
    addtraces!(plt, cz)
    addtraces!(plt, lx)
    addtraces!(plt, ly)
    addtraces!(plt, lz)
    relayout!(plt, scene=attr(
        annotations=[
            attr(
                showarrow=false,
                x=origin[1] + 1.1 * r,
                y=origin[2],
                z=origin[3],
                text="x",
                font=attr(color="red")
            ),
            attr(
                showarrow=false,
                x=origin[1],
                y=origin[2] + 1.1 * r,
                z=origin[3],
                text="y",
                font=attr(color="green")
            ),
            attr(
                showarrow=false,
                x=origin[1],
                y=origin[2],
                z=origin[3] + 1.1 * r,
                text="z",
                font=attr(color="blue")
            ),
        ]
    ))
end

"""
    no_grids(plt)

Disable grid for a plot.
"""
function no_grids(plt)
    relayout!(plt,
        scene=attr(
            xaxis=attr(
                visible=false,
                showgrid=false
            ),
            yaxis=attr(
                visible=false,
                showgrid=false
            ),
            zaxis=attr(
                visible=false,
                showgrid=false
            ),
        ))
end

"""
    blank_layout()

Return blank layout.
"""
function blank_layout()
    layout = Layout(
        scene=attr(
            xaxis=attr(
                visible=false,
                showgrid=false
            ),
            yaxis=attr(
                visible=false,
                showgrid=false
            ),
            zaxis=attr(
                visible=false,
                showgrid=false
            ),
            aspectmode="data",
        ),
    )
    return layout
end