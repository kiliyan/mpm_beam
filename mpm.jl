# main.jl

import Pkg
Pkg.activate(@__DIR__)

using Dates
using CSV, DataFrames, XLSX, LinearAlgebra, Tables
using StaticArrays
using Printf
using WriteVTK, VTKDataTypes
using WriteVTK: MeshCell, VTKCellTypes

include("helper.jl")
include("beam_helper.jl")


function main()

    println("Run at ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))

    # ----------------- Parameters -----------------
    n_particles = 1 # total number of beam elements
    h = 10.0       # grid spacing
    L = 10.0 - h / 10     # beam length 
    α = 0.0       # FLIP fraction
    clearance = Vec2(h / 20, -h / 2)  # initial clearance from origin
    fixed_edge = 0     # fixed BC at 0 but beam starts after clearance i.e h/20

    E = 1e6       # Young's modulus
    ρ0 = E / 1e5  # density

    Δt = 0.0001
    T_final = 200.0
    T_ramp = T_final * 0.85
    use_smoothstep = true

    h_rect = 0.1
    t_rect = 1.2
    A = 0.12        # cross-sectional area
    Iz = 1e-4       # second moment

    load = -1e-3    # vertical tip load (downwards)

    # Tip particle index (cantilever, load at free end)
    tip_id = n_particles

    # ----------------- Grid Initialisation -----------------
    grid, x_coords, y_coords = make_grid(h, L) #starting and ending postions of grid set in the function
    nx, ny = size(grid)

    # ----------------- Particles Initialisation -----------------
    particles = make_beam_particles(L, A, h_rect, ρ0, n_particles, clearance)

    K_local = beam_local_stiffness_new(E, A, Iz, L / n_particles)  ### local MATRIX DEFINED ONCE
    for p in particles
        p.m = ρ0 * p.l_initial * A
        p.I = p.m * ((p.l_initial^2 + h_rect^2) / 12.0)
        p.T, p.l = transform_local_to_global_2d_new(p.c1_x, p.c2_x)
        p.K_g = frame_K_global_2d_new(p.T, E, A, Iz, p.l, K_local)
    end

    #************************ LOGGING. CAN SKIP TO MAIN LOOP (line 130)************************************************************************
    # ----------------- Output directory + logs -----------------
    dir = "Cantilever_point_load_$(Dates.format(now(), "yyyy-mm-dd_HHMMSS"))"
    outdir = mkpath(joinpath("MPM_results", dir))
    logdir = mkpath(joinpath(outdir, "logs"))
    logfile = joinpath(logdir, "particle_logs.csv")
    if isfile(logfile)
        rm(logfile; force=true)
    end

    vtk_grid_dir = mkpath(joinpath(outdir, "Paraview"))
    vtk_part_dir = mkpath(joinpath(outdir, "Paraview"))

    # Parameter log
    param_logfile = joinpath(logdir, "simulation_parameters.log")
    if isfile(param_logfile)
        rm(param_logfile; force=true)
    end

    function log_parameters(filename::AbstractString; kwargs...)
        open(filename, "a") do io
            println(io, "-"^60)
            println(io, "Run at ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
            for (name, value) in kwargs
                println(io, rpad(String(name), 20), " = ", value)
            end
            println(io)
        end
    end

    # use first particle for reporting I_alpha
    p_ref = particles[1]

    CFL_t = round(h / sqrt(E / ρ0); sigdigits=3)
    axial_t = round(2L * sqrt(ρ0 / E); sigdigits=3)
    bending_t = round(2L^2 * sqrt(ρ0 * A / (E * Iz)) / 10; sigdigits=3)
    w_tip = load * L^3 / (3 * E * Iz)

    log_parameters(param_logfile;
        h=h, T=T_final, Δt=Δt, E=E, ρ⁰=ρ0, L=L, A=A, Iz=Iz, I_alpha=p_ref.I,
        α=α, T_ramp=T_ramp, CFL_t=CFL_t, axial_t=axial_t,
        bending_t=bending_t, load=load, length_particles=n_particles, w_tip=w_tip
    )

    io = open(logfile, "a")
    println(io,
        "step,time,particle_id,force_applied," *
        "c1_fx,c1_fy,c1_Mz," *
        "c2_fx,c2_fy,c2_Mz," *
        "c1_vx,c1_vy," *
        "c2_vx,c2_vy," *
        "c1_angvz,c2_angvz," *
        "c1_x,c1_y," *
        "c2_x,c2_y," *
        "c1_thz,c2_thz"
    )

    # Logging cadence
    n_steps = Int(round(T_final / Δt))
    per_record = max(1, n_steps ÷ 2000) # to make sure we log 2000 steps at most

    # ----------------- Time loop -----------------
    t = 0.0
    step = 0

    #********************************** MAIN LOOP STARTS HERE **************************************************

    while t < T_final + 1e-12

        # **************** Particle to grid (P2G) ****************
        # (1) Zero grid masses & forces
        for i in 1:nx, j in 1:ny
            node = grid[i, j]

            node.m = 0.0
            node.I = 0.0
            node.f = Vec2(0.0, 0.0)
            node.M = 0.0
            node.v = Vec2(0.0, 0.0)
            node.ang_v = 0.0
        end


        r = ramp_factor(t, T_ramp; use_smoothstep=use_smoothstep)
        # (2) Interpolate particles -> grid (P2G)
        for (pid, p) in pairs(particles)
            # precompute positions for weighting
            x_center = p.x
            x_c1 = p.c1_x
            x_c2 = p.c2_x

            for i in 1:nx, j in 1:ny
                node = grid[i, j]
                x_node = node.x

                w_center = hat2d(x_node, x_center, h)
                wI = hat2d(x_node, x_c1, h)
                wJ = hat2d(x_node, x_c2, h)

                if w_center == 0.0 && wI == 0.0 && wJ == 0.0
                    continue
                end

                node.m += w_center * p.m
                node.I += w_center * p.I
                node.v += w_center * p.m * p.v
                node.ang_v += w_center * p.I * p.ang_v
                node.f += wI * p.c1_f + wJ * p.c2_f
                node.M += wI * p.c1_M + wJ * p.c2_M

                # external (ramped) vertical load applied at tip particle center
                if pid == tip_id
                    node.f += wJ * r * Vec2(0.0, load)
                end
            end
        end

        #**************** Grid update ****************
        for i in 1:nx, j in 1:ny
            node = grid[i, j]

            if node.m > 0
                node.v_old = node.v / node.m
            else
                node.v_old = Vec2(0.0, 0.0)
                node.v = Vec2(0.0, 0.0)
            end

            if node.I > 0
                node.ang_v_old = node.ang_v / node.I
            else
                node.ang_v_old = 0.0
                node.ang_v = 0.0
            end

            # explicit time integration
            if node.m > 0
                node.v = node.v_old + (node.f / node.m) * Δt
            end
            if node.I > 0
                node.ang_v = node.ang_v_old + (node.M / node.I) * Δt
            end

            # Fixed BC at x = fixed_edge
            if node.x[1] <= fixed_edge
                node.v = Vec2(0.0, 0.0)
                node.ang_v = 0.0
                node.v_old = Vec2(0.0, 0.0)
                node.ang_v_old = 0.0
            end
        end

        #******************* Grid-to-particle (G2P): PIC + FLIP ***************
        for p in particles
            p.v_pic = Vec2(0.0, 0.0)
            p.ang_v_pic = 0.0
            p.v_flip = Vec2(0.0, 0.0)
            p.ang_v_flip = 0.0
            p.c1_v = Vec2(0.0, 0.0)
            p.c2_v = Vec2(0.0, 0.0)
            p.c1_ang_v = 0.0
            p.c2_ang_v = 0.0

            # recompute weights with updated positions
            x_center = p.x
            x_c1 = p.c1_x
            x_c2 = p.c2_x

            for i in 1:nx, j in 1:ny
                node = grid[i, j]
                x_node = node.x

                w_center = hat2d(x_node, x_center, h)
                wI = hat2d(x_node, x_c1, h)
                wJ = hat2d(x_node, x_c2, h)

                if w_center == 0.0 && wI == 0.0 && wJ == 0.0
                    continue
                end

                # PIC
                p.v_pic += w_center * node.v
                p.ang_v_pic += w_center * node.ang_v

                # FLIP increments
                p.v_flip += w_center * (node.v - node.v_old)
                p.ang_v_flip += w_center * (node.ang_v - node.ang_v_old)

                # corner velocities are PIC
                p.c1_v += wI * node.v
                p.c2_v += wJ * node.v
                p.c1_ang_v += wI * node.ang_v
                p.c2_ang_v += wJ * node.ang_v
            end

            # Blend PIC / FLIP
            p.v = (1 - α) * p.v_pic + α * p.v_flip
            p.ang_v = (1 - α) * p.ang_v_pic + α * p.ang_v_flip
        end

        # (5) Update endpoints (PIC) and centers for each particle
        for p in particles
            p.c1_x += p.c1_v * Δt
            p.c2_x += p.c2_v * Δt
            p.c1_theta += p.c1_ang_v * Δt
            p.c2_theta += p.c2_ang_v * Δt
            p.x = (p.c1_x + p.c2_x) / 2
        end

        # (6) Beam internal forces from K_g, using incremental displacements
        for p in particles
            dx1 = p.c1_v * Δt
            dx2 = p.c2_v * Δt
            dθ1 = p.c1_ang_v * Δt
            dθ2 = p.c2_ang_v * Δt

            p.T, p.l = transform_local_to_global_2d_new(p.c1_x, p.c2_x)
            p.K_g = frame_K_global_2d_new(p.T, E, A, Iz, p.l, K_local)  # This function reuses the local stiffness defined initially along with T

            p.Δq_g = pack_q_new(dx1, dθ1, dx2, dθ2)
            p.Δf_g = p.K_g * p.Δq_g

            p.Δc1_f, p.Δc1_M, p.Δc2_f, p.Δc2_M = unpack_end_forces_new(p.Δf_g)

            p.c1_f += p.Δc1_f
            p.c1_M += p.Δc1_M
            p.c2_f += p.Δc2_f
            p.c2_M += p.Δc2_M
        end

        #***************************** LOGGING RELATED AGAIN*****************************************************************************
        # (7) Logging – one row per particle
        if (step % (per_record) == 0) || step in (0:500) # per record (2000) + 500 steps recorded
            for (pid, p) in pairs(particles)
                force_applied = (pid == tip_id) ? (r * load) : 0.0

                f1 = p.c1_f
                f2 = p.c2_f
                M1 = p.c1_M
                M2 = p.c2_M
                v1 = p.c1_v
                v2 = p.c2_v
                x1 = p.c1_x
                x2 = p.c2_x
                angth1 = p.c1_ang_v
                angth2 = p.c2_ang_v
                th1 = p.c1_theta
                th2 = p.c2_theta

                println(io,
                    "$(step),$(t),$(pid),$(force_applied)," *
                    "$(f1[1]),$(f1[2])," *
                    "$(M1)," *
                    "$(f2[1]),$(f2[2])," *
                    "$(M2)," *
                    "$(v1[1]),$(v1[2])," *
                    "$(v2[1]),$(v2[2])," *
                    "$(angth1),$(angth2)," *
                    "$(x1[1]),$(x1[2])," *
                    "$(x2[1]),$(x2[2])," *
                    "$(th1)," *
                    "$(th2)"
                )
            end

            if (step % (per_record * 10) == 0 || step in (0:500)) # per record/10 (i.e 2000/10) + 500 steps recorded
                write_grid_vtk(step, grid, x_coords, y_coords, vtk_grid_dir)
                write_particles_vtk(step, particles, vtk_part_dir)
            end
        end

        progress_bar(step, n_steps)
        step += 1
        t += Δt
    end

    flush(io)
    close(io)

    # ----------------- CSV -> XLSX -----------------
    df = CSV.read(logfile, DataFrame)
    groups = groupby(df, :particle_id)
    triples = [
        begin
            names_g = names(g)
            cols = [Vector(g[!, nm]) for nm in names_g]
            hdrs = String.(names_g)
            ("p" * string(g.particle_id[1]), cols, hdrs)
        end
        for g in groups
    ]
    xlsx_out = replace(logfile, ".csv" => ".xlsx")
    XLSX.writetable(xlsx_out, triples; overwrite=true)


    # ----- Generate PVD index files -----
    write_pvd(joinpath(vtk_grid_dir, "grid.pvd"), "grid_", vtk_grid_dir)
    write_pvd(joinpath(vtk_part_dir, "particles.pvd"), "particles_", vtk_part_dir)

    println("\n finished logs and paraview files. Note : refer xlsx file with each sheet per particle (instead of csv) ")

    println("Expected tip displacement = ", (load * L^3) / (3 * E * Iz))


    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end


