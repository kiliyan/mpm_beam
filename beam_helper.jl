

# Rotation matrix 2x2
function transform_local_to_global_2d_new(c1::Vec2, c2::Vec2)
    t = c2 - c1
    L = norm(t)
    x̂ = t / L
    ŷ = Vec2(-x̂[2], x̂[1])   # perpendicular in plane
    # Rows are local basis in global coords: [x̂ᵀ; ŷᵀ]
    T = transpose(Mat22(hcat(x̂, ŷ)))
    return T, L
end

# Local Euler–Bernoulli beam stiffness (2D) in local coords
function beam_local_stiffness_new(E::Float64, A::Float64, Iz::Float64, L::Float64)
    EA_L = E * A / L
    EI = E * Iz
    L2 = L^2
    L3 = L^3

    K = zeros(6, 6)

    # Axial (u)
    K[1, 1] = EA_L
    K[1, 4] = -EA_L
    K[4, 1] = -EA_L
    K[4, 4] = EA_L

    # Bending (v, θ)
    K[2, 2] = 12EI / L3
    K[2, 3] = 6EI / L2
    K[2, 5] = -12EI / L3
    K[2, 6] = 6EI / L2

    K[3, 2] = 6EI / L2
    K[3, 3] = 4EI / L
    K[3, 5] = -6EI / L2
    K[3, 6] = 2EI / L

    K[5, 2] = -12EI / L3
    K[5, 3] = -6EI / L2
    K[5, 5] = 12EI / L3
    K[5, 6] = -6EI / L2

    K[6, 2] = 6EI / L2
    K[6, 3] = 2EI / L
    K[6, 5] = -6EI / L2
    K[6, 6] = 4EI / L

    return Mat66(K)
end

# Rotate local K to global using 2D rotation T
function frame_K_global_2d_new(T::Mat22, E::Float64, A::Float64, Iz::Float64, L::Float64, K_local::Mat66)
    # ex = first ROW of T (since we stored basis in rows)
    ex = Vec2(T[1, 1], T[1, 2])
    c = ex[1]
    s = ex[2]

    # Rotation matrix for [u1, v1, θ1, u2, v2, θ2]
    R = @SMatrix [c  s  0  0  0  0;
                  -s c  0  0  0  0;
                  0  0  1  0  0  0;
                  0  0  0  c  s  0;
                  0  0  0 -s  c  0;
                  0  0  0  0  0  1]

    return transpose(R) * K_local * R
end

# Pack/Unpack 6-DOF beam end resultants
@inline function pack_q_new(dx1::Vec2, dθ1::Float64, dx2::Vec2, dθ2::Float64)::Vec6
    return @SVector [dx1[1], dx1[2], dθ1, dx2[1], dx2[2], dθ2]
end

@inline function unpack_end_forces_new(fg::Vec6)
    c1_f = @SVector [fg[1], fg[2]]
    c1_M = fg[3]
    c2_f = @SVector [fg[4], fg[5]]
    c2_M = fg[6]
    return c1_f, c1_M, c2_f, c2_M
end
