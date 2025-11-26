
import Pkg
Pkg.activate(@__DIR__)

using Dates
using CSV, DataFrames, XLSX, LinearAlgebra, Tables
using StaticArrays
using Printf
#using Plots

include("helper.jl") # require plot_tip_deflection
include("beam_helper.jl") # require transform_local_to_global_2d_new,beam_local_stiffness_new,frame_K_global_2d_new,ramp_factor


function main(;
    E::Float64=1e6,
    ν::Float64=0.3,
    ρ⁰::Float64=10.0,
    L::Float64=9.0,  # 10 - 1/10 is what we have used in MPM (L - h/10)
    A::Float64=0.12,
    h_rect::Float64=0.1, # height of the rectangular cross section of beam
    Iz::Float64=1e-4,  # bending moment inertia . mass moment calculated later
    load::Float64=-1e-3,   # tip vertical load
    use_smoothstep::Bool=true,
    Δt::Float64=0.0001,
    T::Float64=200.0,
    T_ramp::Float64=200 *0.85
)

    c1_x0 = Vec(0.0, 0.0)
    c2_x0 = Vec(L, 0.0)

    # Reuse your transform_local_to_global_2d
    T0, L0 = transform_local_to_global_2d_new(c1_x0, c2_x0)

    # --------------- Stiffness: reuse your frame_K_global_2d ---------------
    K_local = beam_local_stiffness_new(E, A, Iz, L)
    K_g = frame_K_global_2d_new(T0, E, A, Iz, L0, K_local)
    #K_g = Matrix(K_raw)    # ensure it's a standard Matrix{Float64}

    # --------------- Lumped mass / inertia per DOF ---------------
    m_tot = ρ⁰ * A * L
    I_tot = m_tot * ((L^2 + h_rect^2) / 12.0)   #mass moment of inertia
    m_node = m_tot / 2
    I_node = I_tot / 2

    # DOF ordering:
    # q = [u1x, u1y, θ1z, u2x, u2y, θ2z]
    M_lumped = Diagonal([m_node, m_node, I_node, m_node, m_node, I_node])

    ndof = 6
    fixed_dofs = [1, 2, 3] # clamp node 1
    free_dofs = setdiff(1:ndof, fixed_dofs)

    # --------------- State vectors (mutable) ---------------
    u = zeros(ndof)   # displacements
    v = zeros(ndof)   # velocities
    a = zeros(ndof)   # accelerations

    nsteps = Int(floor(T / Δt)) + 1
    t_hist = [(i - 1) * Δt for i in 1:nsteps]
    tip_deflection = zeros(Float64, nsteps)  # u2y history (DOF 5)

    # scratch arrays
    f_ext = zeros(ndof)
    f_int = zeros(ndof)
    R = zeros(ndof)
    f_ext_hist = zeros(Float64, nsteps)
    stored_rows = NamedTuple[]

    # --------------- Time integration loop ---------------
    for (n, t) in enumerate(t_hist)

        # --- update geometry from current displacement ---
        c1_x = Vec(0.0, 0.0)
        c2_x = Vec(L + u[4], u[5])        # node 2 current position
        T, L_curr = transform_local_to_global_2d_new(c1_x, c2_x)

        # --- recompute stiffness with updated orientation ---
        K_g = frame_K_global_2d_new(T, E, A, Iz, L0, K_local)

        # external force (only tip vertical DOF u2y)
        fill!(f_ext, 0.0)
        r = ramp_factor(t, T_ramp; use_smoothstep=use_smoothstep)
        f_ext[5] = r * load   # vertical DOF at node 2
        f_ext_hist[n] = f_ext[5] # for logging

        # internal force
        mul!(f_int, K_g, u)   # f_int = K_g * u

        @. R = f_ext - f_int

        # BC: zero residual at fixed DOFs
        R[fixed_dofs] .= 0.0

        # acceleration: a = M_lumped \ R
        a .= M_lumped \ R

        # explicit update on free DOFs
        v[free_dofs] .+= Δt .* a[free_dofs]
        u[free_dofs] .+= Δt .* v[free_dofs]

        # store tip deflection (node 2, y DOF)
        tip_deflection[n] = u[5]

    end

    println("Expected tip displacement = ", (load * L^3) / (3 * E * Iz))
    println("Observed tip deflection at the final timestep = ", tip_deflection[end])

    # pick 2000 evenly-spaced indices between 1 and N_tot
    Nkeep = 2000
    N_tot = length(tip_deflection)
    idx = round.(Int, range(1, N_tot; length=min(Nkeep, N_tot)))

    df = DataFrame(time=t_hist[idx], tip_deflection=tip_deflection[idx],f_ext= f_ext_hist[idx])

    cols = [df.time, df.tip_deflection,df.f_ext]
    hdrs = ["time", "tip_deflection","f_ext"]
    sheet = ("tip_history", cols, hdrs)
    triples = [sheet]


    outdir = mkpath("fem_results")
    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    xlsx_file = joinpath(outdir, "fem_tip_deflection_" * timestamp * ".xlsx")

    XLSX.writetable(xlsx_file, triples; overwrite=true)
    println(" Tip deflection saved to XLSX: ", xlsx_file)

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

