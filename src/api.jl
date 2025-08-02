"""
    cuboids(origin::Vector{<:Real}, dimension::Vector{<:Real}, color::String=""; opc::Real=1)

    Creates a 3D box mesh centered at the given origin with specified dimensions and color.

    # Arguments
    - `origin::Vector{<:Real}`: A vector of three Reals specifying the center of the box.
    - `dimension::Vector{<:Real}`: A vector of three Reals specifying the dimensions (width, height, depth) of the box.
    - `color::String`: A string specifying the color of the box.
    
    # Keywords
    - `opc`: (optional) A Real specifying the opacity of the box. Default is 1.                                   
"""
function cuboids(origin::Vector{<:Real}, dimension::Vector{<:Real}, color::String=""; opc::Real=1)
    @assert length(origin) == 3
    @assert length(dimension) == 3

    if color == ""
        @all r g b = round(Int, rand() * 255)
        color = "rgb($r, $g, $b)"
    end

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
    cubes(origin::Vector{<:Real}, side::Real, color::String=""; opc::Real=1)

    Creates a 3D cube mesh centered at the given origin with specified dimensions and color.

    # Arguments
    - `origin::Vector{<:Real}`: A vector of three Reals specifying the center of the cube.
    - `side::Real`: Side length of the cube.
    - `color::String`: A string specifying the color of the cube.
    
    # Keywords
    - `opc`: (optional) A Real specifying the opacity of the cube. Default is 1.                                   
"""
function cubes(origin::Vector{<:Real}, side::Real, color::String=""; opc::Real=1)
    @assert length(origin) == 3
    @assert side > 0

    return cuboids(origin, [side, side, side], color; opc=opc)
end

"""
    squares(origin::Vector{<:Real}, side::Real, mode::String="z", color::String=""; opc::Real=1)
    
    Creates a 2D square mesh centered at the given origin with the specified side length and color.

    # Arguments
    - `origin::Vector{<:Real}`: A vector of three Reals specifying the center of the square.
    - `side::Real`: A Real specifying the side length of the square.
    - `mode`::String: (optional) A string specifying the orientation of the square ("x", "y", or "z"). Default is "z".
    - `color::String`: A string specifying the color of the square.

    
    # Keywords
    - `opc`: (optional) A Real specifying the opacity of the square. Default is 1.
"""
function squares(origin::Vector{<:Real}, side::Real, mode::String="z", color::String=""; opc::Real=1)
    @assert length(origin) == 3

    if color == ""
        @all r g b = round(Int, rand() * 255)
        color = "rgb($r, $g, $b)"
    end

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
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0
        ),
    )
end

"""
    ellipsoids(origin::Vector{<:Real}, par::Vector{<:Real}, color::String=""; opc::Real=1, tres=61, pres=31, ah::Real=0)

    Creates a 3D ellipsoid mesh.

    # Arguments
    - `origin::Vector{<:Real}`: The center of the ellipsoid.
    - `par::Vector{<:Real}`: Parameters of the ellipsoid (a, b, c).
    - `color::String`: The color of the ellipsoid.
    
    # Keywords
    - `opc`: The opacity of the ellipsoid. Default is 1.
    - `tres`: The resolution of the mesh grid (theta). Default is 61.
    - `pres`: The resolution of the mesh grid (phi). Default is 31.
    - `ah`: alphahull value. Default is 0.
"""
function ellipsoids(origin::Vector{<:Real}, par::Vector{<:Real}, color::String=""; opc::Real=1, tres=61, pres=31, ah::Real=0)
    @assert length(origin) == 3
    @assert length(par) == 3
    @assert all(par .> 0)
    @assert tres > 0
    @assert pres > 0

    if color == ""
        @all r g b = round(Int, rand() * 255)
        color = "rgb($r, $g, $b)"
    end

    phi = LinRange(0, 360, pres)
    tht = LinRange(0, 180, tres)

    x = sind.(tht) .* cosd.(phi') .* par[1] .+ origin[1]
    y = sind.(tht) .* sind.(phi') .* par[2] .+ origin[2]
    z = cosd.(tht * ones(length(phi))') .* par[3] .+ origin[3]
    x = x[:]
    y = y[:]
    z = z[:]

    return mesh3d(x=x, y=y, z=z,
        alphahull=ah,
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
    spheres(origin::Vector{<:Real}, r::Real, color::String=""; opc::Real=1, tres=60, pres=30, ah::Real=0)

    Creates a 3D sphere mesh.

    # Arguments
    - `origin::Vector{<:Real}`: The center of the ellipsoid.
    - `r::Real`: Radius of the sphere.
    - `color::String`: The color of the ellipsoid.

    # Keywords
    - `opc`: The opacity of the ellipsoid. Default is 1.
    - `tres`: The resolution of the mesh grid (theta). Default is 60.
    - `pres`: The resolution of the mesh grid (phi). Default is 30.
    - `ah`: alphahull value. Default is 0.
"""
function spheres(origin::Vector{<:Real}, r::Real, color::String=""; opc::Real=1, tres=60, pres=30, ah::Real=0)
    @assert length(origin) == 3
    @assert r > 0

    if color == ""
        @all r g b = round(Int, rand() * 255)
        color = "rgb($r, $g, $b)"
    end

    return ellipsoids(origin, [r, r, r], color; opc=opc, tres=tres, pres=pres, ah=ah)
end

"""
    lines(pt1::Vector{<:Real}, pt2::Vector{<:Real}, color::String; opc::Real=1, style="")

    Creates a 3D line between two points.

    # Arguments
    - `pt1::Vector{<:Real}`: Starting point of the line.
    - `pt2::Vector{<:Real}`: Ending point of the line.
    - `color::String`: The color of the line.

    # Keywords
    - `opc`: The opacity of the line. Default is 1.
    - `style`: The line style (e.g., "solid", "dash"). Default is "".
"""
function lines(pt1::Vector{<:Real}, pt2::Vector{<:Real}, color::String=""; opc::Real=1, style="")
    @assert length(pt1) == 3
    @assert length(pt2) == 3

    if color == ""
        @all r g b = round(Int, rand() * 255)
        color = "rgb($r, $g, $b)"
    end

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
        flatshading=true,
        opacity=opc,
    )
end

"""
    polygons(pts::Vector, color::String; opc::Real=1, ah::Real=0)

    Creates a polygon mesh from a set of points (form around the mid point of the set of points).

    # Arguments
    - `pts::::Vector`: List of points defining the polygon.
    - `color::String`: The color of the polygon.

    # Keywords
    - `opc`: The opacity of the polygon. Default is 1.
    - `ah`: alphahull value. Default is 0.
"""
function polygons(pts::Vector, color::String=""; opc::Real=1, ah::Real=0)
    @assert all(length.(pts) .== 3)

    if color == ""
        @all r g b = round(Int, rand() * 255)
        color = "rgb($r, $g, $b)"
    end

    for vec in pts
        for num in vec
            @assert isreal(num)
        end
    end

    N = length(pts)
    pts_copy = sort_pts(pts)
    @all x y z = []
    for i in eachindex(pts_copy)
        push!(x, pts_copy[i][1])
        push!(y, pts_copy[i][2])
        push!(z, pts_copy[i][3])
    end

    mid = zeros(3)
    for n in eachindex(pts_copy)
        for m in 1:3
            mid[m] += pts_copy[n][m]
        end
    end
    mid .= mid ./ N
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
        alphahull=ah,
        color=color,
        opacity=opc,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0
        ),
    )
end

"""
    polygons(pts::Vector, ng::Int, color::String; opc::Real=1, ah::Real=0)

    Creates a group of polygons from a set of points and a specified number of vertices per polygon.

    # Arguments
    - `pts::Vector`: List of points defining the mesh.
    - `ng::Int`: Number of vertices per polygon.
    - `color::String`: The color of the mesh.

    # Keywords
    - `opc`: The opacity of the mesh. Default is 1.
    - `ah`: alphahull value.
"""
function polygons(pts::Vector, ng::Int, color::String=""; opc::Real=1, ah::Real=0)
    @assert all(length.(pts) .== 3)

    if color == ""
        @all r g b = round(Int, rand() * 255)
        color = "rgb($r, $g, $b)"
    end

    for vec in pts
        for num in vec
            @assert isreal(num)
        end
    end
    @assert ng > 0

    N = length(pts)

    @all x y z i j k = []
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
        mid .= mid ./ ng
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
        alphahull=ah,
        color=color,
        opacity=opc,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0
        ),
    )
end

"""
    gtrans!(geo::GenericTrace, dis::Vector{<:Real})

    Translates a 3D geometry by a specified displacement vector.

    # Arguments
    - `geo::GenericTrace`: The geometry to translate.
    - `dis::Vector{<:Real}`: A vector of three Reals specifying the translation distances for the x, y, and z axes.
"""
function gtrans!(geo::GenericTrace, dis::Vector{<:Real})
    @inbounds for n in eachindex(geo.x)
        geo.x[n] += dis[1]
        geo.y[n] += dis[2]
        geo.z[n] += dis[3]
    end
end

"""
    grot!(geo::GenericTrace, rotang::Vector{<:Real}, center::Vector{<:Real}=[0])

    Rotates a 3D geometry around a specified center point. (Tait–Bryan rotation)

    # Arguments
    - `geo::GenericTrace`: The 3D geometry to be rotated, which must have `x`, `y`, and `z` coordinates.
    - `rotang::Vector{<:Real}`: A vector of three Tait–Bryan rotation angles in degrees  for rotations around the x, y, and z axes respectively.
    - `center::Vector{<:Real}`: The center point of rotation. Default is `[0]`, which means the rotation center will be set at the geometric center of the object.
"""
function grot!(geo::GenericTrace, rotang::Vector{<:Real}, center::Vector{<:Real}=[0])
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
    grot!(geo::GenericTrace, ang::Real, axis::Vector{<:Real}, origin::Vector{<:Real}=[0])

    Rotates the geometry by the specified angle `ang` around the axis `axis` and origin `origin`.
    
    # Arguments
    - `geo::GenericTrace`: The geometry to be rotated.
    - `ang::Real`: The rotation angle.
    - `axis::Vector{<:Real}`: The rotation axis.
    - `origin::Vector{<:Real}=[0]`: The rotation origin. Defaults to the center of the geometry if not specified.
"""
function grot!(geo::GenericTrace, ang::Real, axis::Vector{<:Real}, origin::Vector{<:Real}=[0])
    @assert length(axis) == 3

    pos = []
    for n in eachindex(geo.x)
        push!(pos, [geo.x[n], geo.y[n], geo.z[n]])
    end

    axis_norm = axis ./ norm(axis)
    vrot = similar(pos)

    if origin == [0] # rotation center set at the geometry center
        origin = sum(pos) ./ length(pos)
    else
        @assert length(origin) == 3
    end

    v = similar(origin)
    for n in eachindex(vrot)
        v .= (pos[n] .- origin)
        vrot[n] = cosd(ang) * v + sind(ang) * cross(axis_norm, v) + (1 - cosd(ang)) * dot(axis_norm, v) * axis_norm
        pos[n] = vrot[n] .+ origin
    end

    geo.x = getindex.(pos, 1)
    geo.y = getindex.(pos, 2)
    geo.z = getindex.(pos, 3)
end

"""
    sort_pts(pts::Vector)

    Sorts points in place based on their angular position relative to the centroid.

    # Arguments
    - `pts::Vector`: List of points to be sorted.
"""
function sort_pts(pts::Vector)
    @assert all(length.(pts) .== 3)

    N = length(pts)
    ang = zeros(length(pts))
    mid = zeros(3)
    for n in eachindex(pts)
        for m in 1:3
            mid[m] += pts[n][m]
        end
    end
    mid .= mid ./ N

    c = collect(combinations(1:N, 3))
    thtr = 0
    phir = 0

    vec = similar(mid)
    for n in eachindex(c)
        vec .= cross(pts[c[n][2]] .- pts[c[n][1]], pts[c[n][3]] .- pts[c[n][1]])
        if norm(vec) == 0
            continue
        else
            vec .= vec ./ norm(vec)
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
    pts_rot .= pts[sortperm(ang)]

    return pts_rot
end

"""
    sort_pts!(pts::Vector)

    Sorts points in place based on their angular position relative to the centroid.

    # Arguments
    - `pts::Vector`: List of points to be sorted.
"""
function sort_pts!(pts::Vector)
    @assert all(length.(pts) .== 3)

    N = length(pts)
    ang = zeros(length(pts))
    mid = zeros(3)
    for n in eachindex(pts)
        for m in 1:3
            mid[m] += pts[n][m]
        end
    end
    mid .= mid ./ N

    c = collect(combinations(1:N, 3))
    thtr = 0
    phir = 0
    vec = similar(mid)
    for n in eachindex(c)
        vec .= cross(pts[c[n][2]] .- pts[c[n][1]], pts[c[n][3]] .- pts[c[n][1]])
        if norm(vec) == 0
            continue
        else
            vec .= vec ./ norm(vec)
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
    return nothing
end

"""
    add_ref_axes!(plt::PlotlyJS.SyncPlot, origin::Vector{<:Real}=[0, 0, 0], r::Real=1)

    Adds reference axes (x, y, z) to a plot.

    # Arguments
    - `plt::PlotlyJS.SyncPlot`: The plot to which the axes will be added.
    - `origin::Vector{<:Real}`: The origin point of the axes.
    - `r::Real`: The length of the reference axes.
"""
function add_ref_axes!(plt::PlotlyJS.SyncPlot, origin::Vector{<:Real}=[0, 0, 0], r::Real=1)
    @assert length(origin) == 3
    @assert r > 0

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
    return nothing
end

"""
    add_ref_axes!(plt::PlotlyJS.SyncPlot, origin::Vector{<:Real}, r::Vector{<:Real})

    Adds reference axes (x, y, z) to a plot.

    # Arguments
    - `plt::PlotlyJS.SyncPlot`: The plot to which the axes will be added.
    - `origin::Vector{<:Real}`: The origin point of the axes.
    - `r::Vector{<:Real}`: The lengths of the reference axes.
"""
function add_ref_axes!(plt::PlotlyJS.SyncPlot, origin::Vector{<:Real}, r::Vector{<:Real})
    @assert length(origin) == 3
    @assert length(r) == 3
    @assert all(r .> 0)

    csize = minimum([r[1] / 10, r[2] / 10, r[3] / 10])
    cx = cone(x=[r[1] + origin[1]], y=[origin[2]], z=[origin[3]], u=[csize], v=[0], w=[0],
        colorscale=[[0, "red"], [1, "red"]],
        showscale=false)
    cy = cone(x=[origin[1]], y=[r[2] + origin[2]], z=[origin[3]], u=[0], v=[csize], w=[0],
        colorscale=[[0, "green"], [1, "green"]],
        showscale=false)
    cz = cone(x=[origin[1]], y=[origin[2]], z=[r[3] + origin[3]], u=[0], v=[0], w=[csize],
        colorscale=[[0, "blue"], [1, "blue"]],
        showscale=false)
    lx = scatter3d(x=[origin[1], r[1] + origin[1]], y=[origin[2], origin[2]], z=[origin[3], origin[3]],
        line=attr(color="red"),
        mode="lines",
        showlegend=false)
    ly = scatter3d(x=[origin[1], origin[1]], y=[origin[2], r[2] + origin[2]], z=[origin[3], origin[3]],
        line=attr(color="green"),
        mode="lines",
        showlegend=false)
    lz = scatter3d(x=[origin[1], origin[1]], y=[origin[2], origin[2]], z=[origin[3], r[3] + origin[3]],
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
                x=origin[1] + 1.1 * r[1],
                y=origin[2],
                z=origin[3],
                text="x",
                font=attr(color="red")
            ),
            attr(
                showarrow=false,
                x=origin[1],
                y=origin[2] + 1.1 * r[2],
                z=origin[3],
                text="y",
                font=attr(color="green")
            ),
            attr(
                showarrow=false,
                x=origin[1],
                y=origin[2],
                z=origin[3] + 1.1 * r[3],
                text="z",
                font=attr(color="blue")
            ),
        ]
    ))
    return nothing
end

"""
    add_arrows!(plt::PlotlyJS.SyncPlot, origin::Vector{<:Real}, dir::Vector{<:Real}, len::Real=1.0, color::String=""; opc::Real=1, endpoint::Bool=true, asize::Real=len)
    
    Creates a 3D arrow starting from a point and pointing in a given direction.

    # Arguments
    - `plt::PlotlyJS.SyncPlot`: The plot to which the axes will be added.
    - `origin::Vector{<:Real}`: The starting point of the arrow.
    - `dir::Vector{<:Real}`: The direction vector of the arrow.
    - `len::Real`: length of the arrow
    - `color::String`: The color of the arrow.

    # Keywords
    - `opc`: The opacity of the arrow. Default is 1.
    - `asize` Size factor of the arrow cone. Default is `len`.
"""
function add_arrows!(plt::PlotlyJS.SyncPlot, origin::Vector{<:Real}, dir::Vector{<:Real}, len::Real=1.0, color::String=""; opc::Real=1, endpoint::Bool=true, asize::Real=len)
    @assert length(origin) == 3
    @assert length(dir) == 3
    @assert len > 0
    @assert asize > 0

    if color == ""
        @all r g b = round(Int, rand() * 255)
        color = "rgb($r, $g, $b)"
    end

    dir_norm = dir / norm(dir) * len

    if endpoint == true
        c = cone(x=[origin[1] + dir_norm[1]], y=[origin[2] + dir_norm[2]], z=[origin[3] + dir_norm[3]], u=[dir_norm[1]/len*asize], v=[dir_norm[2]]/len*asize, w=[dir_norm[3]/len*asize],
            colorscale=[[0, color], [1, color]],
            opacity=opc,
            showscale=false)
        l = scatter3d(x=[origin[1], origin[1] + dir_norm[1]], y=[origin[2], origin[2] + dir_norm[2]], z=[origin[3], origin[3] + dir_norm[3]],
            line=attr(color=color, width=2 *dir_norm),
            mode="lines",
            opacity=opc,
            showlegend=false)
    else # origin at center
        c = cone(x=[origin[1] + dir_norm[1] / 2], y=[origin[2] + dir_norm[2] / 2], z=[origin[3] + dir_norm[3] / 2], u=[dir_norm[1]/len*asize], v=[dir_norm[2]/len*asize], w=[dir_norm[3]/len*asize],
            colorscale=[[0, color], [1, color]],
            opacity=opc,
            showscale=false)
        l = scatter3d(x=[origin[1] - dir_norm[1] / 2, origin[1] + dir_norm[1] / 2], y=[origin[2] - dir_norm[2] / 2, origin[2] + dir_norm[2] / 2], z=[origin[3] - dir_norm[3] / 2, origin[3] + dir_norm[3] / 2],
            line=attr(color=color, width=2 * dir_norm),
            mode="lines",
            opacity=opc,
            showlegend=false)
    end

    addtraces!(plt, c)
    addtraces!(plt, l)

    return nothing
end

"""
    add_text!(plt::PlotlyJS.SyncPlot, origin::Vector{<:Real}, text::String, color::String="")

    Add text to plot. 

    # Arguments
    - plt::PlotlyJS.SyncPlot: Plot to add text.
    - origin::Vector{<:Real}: origin of the text.
    - text::String:: text to be added.
    - color::String=: color of the text.
"""
function add_text!(plt::PlotlyJS.SyncPlot, origin::Vector{<:Real}, text::String, color::String="")

    if color == ""
        @all r g b = round(Int, rand() * 255)
        color = "rgb($r, $g, $b)"
    end

    addtraces!(plt, scatter3d(x=[origin[1]], y=[origin[2]], z=[origin[3]],
        mode="text", text=[text], textposition="middle center",
        textfont=attr(
            color=color,
        ),
        showlegend=false,
    ))
    return nothing
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

"""
    set_view!(plt::PlotlyJS.SyncPlot, az::Real, el::Real, r::Real=1)

    Set az/el (deg) view of the plot.

    # Arguments
    - plt::PlotlyJS.SyncPlot: plot to be modified.
    - az::Real: az value.
    - el::Real: el value.
    - r::Real=1: distance.
"""
function set_view!(plt::PlotlyJS.SyncPlot, az::Real, el::Real, r::Real=1.25 * sqrt(3))
    if el == 90
        el -= 0.01
    elseif el == -90
        el += 0.01
    end
    x = r * cosd(el) * cosd(az)
    y = r * cosd(el) * sind(az)
    z = r * sind(el)
    if !haskey(plt.plot.layout.scene_camera, :eye)
        plt.plot.layout.scene_camera = attr(eye=attr(x=1.25, y=1.25, z=1.25))
    end
    plt.plot.layout.scene_camera[:eye][:x] = x
    plt.plot.layout.scene_camera[:eye][:y] = y
    plt.plot.layout.scene_camera[:eye][:z] = z
    react!(plt, plt.plot.data, plt.plot.layout)
    return nothing
end


