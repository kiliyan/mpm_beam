# helper.jl
using StaticArrays

const Vec2 = SVector{2,Float64}
const Mat22 = SMatrix{2,2,Float64,4}
const Vec6 = SVector{6,Float64}
const Mat66 = SMatrix{6,6,Float64,36}
const Vec = Vec2


smoothstep(s) = 3s^2 - 2s^3
function ramp_factor(t, T_ramp; use_smoothstep::Bool=true)
    if T_ramp <= 0
        return 1.0
    end
    s = t / T_ramp
    if s <= 0
        return 0.0
    elseif s >= 1
        return 1.0
    else
        return use_smoothstep ? smoothstep(s) : s
    end
end

# -----------------------------------------
# 1D & 2D hat kernels (shape functions)
# -----------------------------------------
@inline function hat1d(x_node::Float64, x_point::Float64, h::Float64)
    ξ = (x_point - x_node) / h
    r = abs(ξ)
    if r <= 1.0
        return 1.0 - r
    else
        return 0.0
    end
end

@inline function hat2d(x_node::Vec2, x_point::Vec2, h::Float64)
    wx = hat1d(x_node[1], x_point[1], h)
    wy = hat1d(x_node[2], x_point[2], h)
    return wx * wy
end


# Grid node for 2D grid
mutable struct GridNode
    x::Vec2
    v::Vec2
    v_old::Vec2
    ang_v::Float64
    ang_v_old::Float64
    m::Float64
    I::Float64
    f::Vec2
    M::Float64
end

function make_grid(h::Float64, L::Float64)
    # 2D grid around the beam
    x_coords = -2h:h:L+3h            #######   decide the range of the grid in x
    y_coords = -3h:h:h               #######   decide the range of the grid in y
    nx = length(x_coords)
    ny = length(y_coords)

    grid = Array{GridNode}(undef, nx, ny)
    for i in 1:nx
        for j in 1:ny
            x = Vec2(x_coords[i], y_coords[j])
            grid[i, j] = GridNode(
                x,
                Vec2(0.0, 0.0),
                Vec2(0.0, 0.0),
                0.0, 0.0,
                0.0, 0.0,
                Vec2(0.0, 0.0),
                0.0
            )
        end
    end
    return grid, x_coords, y_coords
end

#Beam particle struct
mutable struct BeamParticle
    x::Vec2      # center
    m::Float64
    I::Float64
    v::Vec2
    v_pic::Vec2
    v_flip::Vec2
    ang_v::Float64
    ang_v_pic::Float64
    ang_v_flip::Float64

    x_initial::Vec2

    l::Float64
    l_initial::Float64
    c1_x::Vec2
    c2_x::Vec2
    c1_theta::Float64
    c2_theta::Float64

    c1_v::Vec2
    c2_v::Vec2
    c1_ang_v::Float64
    c2_ang_v::Float64

    c1_M::Float64
    c2_M::Float64
    c1_f::Vec2
    c2_f::Vec2

    Δc1_M::Float64
    Δc2_M::Float64
    Δc1_f::Vec2
    Δc2_f::Vec2

    T::Mat22
    T_new::Mat22
    Δf_g::Vec6
    Δq_g::Vec6
    K_g::Mat66

    center_flag::Bool
end

function make_beam_particles(L::Float64,
    A::Float64,
    h_rect::Float64,
    ρ0::Float64,
    n_particles::Int, clearance::Vec2=Vec2(0.0, 0.0))

    @assert n_particles ≥ 1 "n_particles must be at least 1"

    ℓ = L / n_particles
    particles = Vector{BeamParticle}(undef, n_particles)

    for k in 1:n_particles
        # element k spans [(k-1)*ℓ, k*ℓ] along x
        x1 = clearance[1] + (k - 1) * ℓ
        x2 = clearance[1] + k * ℓ

        c1_x = Vec2(x1, clearance[2])
        c2_x = Vec2(x2, clearance[2])
        x_center = Vec2((x1 + x2) / 2, 0.0)

        l = ℓ
        l_initial = ℓ

        m = ρ0 * l * A
        # approximate out-of-plane inertia for this segment
        I = m * ((l^2 + h_rect^2) / 12.0)

        # Local->global triad and element stiffness (placeholder E=1, Iz=1)
        T, L0 = transform_local_to_global_2d_new(c1_x, c2_x)
        K_local = beam_local_stiffness_new(1.0, A, 1.0, L0)
        K_g = frame_K_global_2d_new(T, 1.0, A, 1.0, L0, K_local)  # overwrite later in main

        particles[k] = BeamParticle(
            x_center, m, I,
            Vec2(0.0, 0.0),          # v
            Vec2(0.0, 0.0),          # v_pic
            Vec2(0.0, 0.0),          # v_flip
            0.0, 0.0, 0.0,           # ang_v, ang_v_pic, ang_v_flip
            x_center,                # x_initial
            l, l_initial,
            c1_x, c2_x,
            0.0, 0.0,                # c1_theta, c2_theta
            Vec2(0.0, 0.0), Vec2(0.0, 0.0),  # c1_v, c2_v
            0.0, 0.0,                # c1_ang_v, c2_ang_v
            0.0, 0.0,                # c1_M, c2_M
            Vec2(0.0, 0.0), Vec2(0.0, 0.0),  # c1_f, c2_f
            0.0, 0.0,                # Δc1_M, Δc2_M
            Vec2(0.0, 0.0), Vec2(0.0, 0.0),  # Δc1_f, Δc2_f
            T, T,
            zeros(Vec6), zeros(Vec6),
            K_g,
            false                     # center_flag
        )
    end

    return particles
end

# ----------------- VTK output helpers -----------------
function write_grid_vtk(step::Int,
    grid::Array{GridNode,2},
    x_coords::AbstractVector,
    y_coords::AbstractVector,
    outdir::AbstractString)
    nx, ny = size(grid)
    z_coords = [0.0]  # 2D grid in z=0 plane

    fname = joinpath(outdir, @sprintf("grid_%06d", step))
    vtk = vtk_grid(fname, x_coords, y_coords, z_coords)

    vx = Array{Float64}(undef, nx, ny, 1)
    vy = similar(vx)
    vmag = similar(vx)
    m_arr = similar(vx)

    for i in 1:nx, j in 1:ny
        node = grid[i, j]
        vx[i, j, 1] = node.v[1]
        vy[i, j, 1] = node.v[2]
        vmag[i, j, 1] = sqrt(node.v[1]^2 + node.v[2]^2)
        m_arr[i, j, 1] = node.m
    end

    vtk["vx", VTKPointData()] = vx
    vtk["vy", VTKPointData()] = vy
    vtk["vmag", VTKPointData()] = vmag
    vtk["m", VTKPointData()] = m_arr

    vtk_save(vtk)
end

function write_particles_vtk(step::Int,
    particles::Vector{BeamParticle},
    outdir::AbstractString)

    N = length(particles)
    npts = 3 * N                     # center, c1, c2 for each particle

    # points matrix: size (3, npts), columns = [x, y, z]ᵀ
    points = Array{Float64}(undef, 3, npts)

    vx = Vector{Float64}(undef, npts)
    vy = Vector{Float64}(undef, npts)
    vmag = Vector{Float64}(undef, npts)
    ptype = Vector{Int}(undef, npts)   # 0=center, 1=c1, 2=c2

    for (i, p) in pairs(particles)
        # indices for this particle
        i_center = i
        i_c1 = N + i
        i_c2 = 2N + i

        # ----- center -----
        points[1, i_center] = p.x[1]
        points[2, i_center] = p.x[2]
        points[3, i_center] = 0.0

        vx[i_center] = p.v[1]
        vy[i_center] = p.v[2]
        vmag[i_center] = hypot(p.v[1], p.v[2])
        ptype[i_center] = 0

        # ----- left endpoint (c1) -----
        points[1, i_c1] = p.c1_x[1]
        points[2, i_c1] = p.c1_x[2]
        points[3, i_c1] = 0.0

        vx[i_c1] = p.c1_v[1]
        vy[i_c1] = p.c1_v[2]
        vmag[i_c1] = hypot(p.c1_v[1], p.c1_v[2])
        ptype[i_c1] = 1

        # ----- right endpoint (c2) -----
        points[1, i_c2] = p.c2_x[1]
        points[2, i_c2] = p.c2_x[2]
        points[3, i_c2] = 0.0

        vx[i_c2] = p.c2_v[1]
        vy[i_c2] = p.c2_v[2]
        vmag[i_c2] = hypot(p.c2_v[1], p.c2_v[2])
        ptype[i_c2] = 2
    end

    # Each point is a VTK_VERTEX cell
    cells = [MeshCell(VTKCellTypes.VTK_VERTEX, [i]) for i in 1:npts]

    fname = joinpath(outdir, @sprintf("particles_%06d", step))
    vtk = vtk_grid(fname, points, cells)

    vtk["vx", VTKPointData()] = vx
    vtk["vy", VTKPointData()] = vy
    vtk["vmag", VTKPointData()] = vmag
    vtk["ptype", VTKPointData()] = ptype   # 0=center,1=c1,2=c2

    vtk_save(vtk)
end

function write_pvd(pvd_path::String, prefix::String, folder::String)
    files = sort(readdir(folder))

    open(pvd_path, "w") do io

        println(io, """<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">""")
        println(io, "  <Collection>")

        timestep = 0
        for f in files
            if startswith(f, prefix)
                println(io, @sprintf("""    <DataSet timestep="%d" file="%s"/>""", timestep, f))
                timestep += 1
            end
        end

        println(io, "  </Collection>")
        println(io, "</VTKFile>")
    end

    #println("Wrote PVD: $pvd_path")
end

# Simple progress bar
function progress_bar(step, total; width=40)
    frac = step / total
    filled = round(Int, frac * width)
    empty = width - filled
    bar = repeat("█", filled) * repeat("─", empty)
    @printf("\r[%s] %3d%%", bar, round(Int, frac * 100))
    flush(stdout)
end


# Plots taking too much time to compile. instead outputting to excel and plotting it there
#=
function plot_tip_deflection(t_hist, w_tip; ylabel="Tip deflection (m)", xlabel="Time (s)", title="Cantilever Tip Deflection")
    plt = plot(
        t_hist,
        w_tip,
        lw=2,
        label="u₂y(t)",
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        legend=:topright,
        grid=true,
    )
    display(plt)
    return plt
end
=#

