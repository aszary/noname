module SolidBody

using LinearAlgebra
using Statistics

export EllipseFit, fit_ellipse, rotate_sparks!, rotate_sparks

"""
    EllipseFit

Pre-computed ellipse fit for a polar cap on a sphere.
"""
struct EllipseFit
    center_local::Vector{Float64}  # ellipse center in tangent plane
    a::Float64                      # semi-major axis
    b::Float64                      # semi-minor axis
    θ::Float64                      # orientation angle of major axis
    x_hat::Vector{Float64}          # local x-axis (tangent plane)
    y_hat::Vector{Float64}          # local y-axis (tangent plane)
    z_hat::Vector{Float64}          # local z-axis (normal to tangent plane)
    centroid::Vector{Float64}       # cap centroid on sphere
    R::Float64                      # sphere radius
end

"""
    fit_ellipse(boundary_points::Matrix{Float64}, R::Float64) -> EllipseFit

Fit an ellipse to polar cap boundary points lying on a sphere of radius R.
boundary_points: 3×N matrix of (x, y, z) on the sphere surface.

Method:
  1. Compute centroid and build local tangent plane frame
  2. Project boundary points onto tangent plane
  3. Algebraic conic fit: Au² + Buv + Cv² + Du + Ev = 1
  4. Extract semi-axes, orientation, center from conic matrix
"""
function fit_ellipse(boundary_points::Matrix{Float64}, R::Real)
    N = size(boundary_points, 2)
    @assert N >= 5 "Need at least 5 boundary points to fit an ellipse"

    # Centroid on sphere
    centroid = vec(mean(boundary_points, dims=2))
    centroid = centroid / norm(centroid) * R

    # Local tangent plane basis
    z_hat = centroid / R
    ref = abs(z_hat[1]) < 0.9 ? [1.0, 0.0, 0.0] : [0.0, 1.0, 0.0]
    x_hat = normalize(cross(ref, z_hat))
    y_hat = cross(z_hat, x_hat)

    # Project onto tangent plane
    u = zeros(N)
    v = zeros(N)
    for i in 1:N
        d = boundary_points[:, i] - centroid
        u[i] = dot(d, x_hat)
        v[i] = dot(d, y_hat)
    end

    # Algebraic ellipse fit: Au² + Buv + Cv² + Du + Ev = 1
    D_mat = hcat(u .^ 2, u .* v, v .^ 2, u, v)
    coeffs = D_mat \ ones(N)
    A, B, C, D_c, E_c = coeffs

    # Conic matrix and ellipse center
    M = [A B/2; B/2 C]
    center_local = -M \ [D_c / 2, E_c / 2]

    # Semi-axes and orientation
    evals, evecs = eigen(M)
    F = 1.0 + dot([D_c, E_c], center_local) / 2

    a = sqrt(abs(F / evals[1]))
    b = sqrt(abs(F / evals[2]))

    if a < b
        a, b = b, a
        evecs = evecs[:, [2, 1]]
    end

    θ = atan(evecs[2, 1], evecs[1, 1])

    return EllipseFit(center_local, a, b, θ, x_hat, y_hat, z_hat, centroid, R)
end

# ── Internal helpers ──

function _to_local_2d(point_3d, ef::EllipseFit)
    d = point_3d - ef.centroid
    return dot(d, ef.x_hat), dot(d, ef.y_hat)
end

function _to_3d_on_sphere(u, v, ef::EllipseFit)
    p = ef.centroid + u * ef.x_hat + v * ef.y_hat
    return p / norm(p) * ef.R
end

function _to_ellipse_frame(u, v, ef::EllipseFit)
    du = u - ef.center_local[1]
    dv = v - ef.center_local[2]
    cosθ = cos(ef.θ)
    sinθ = sin(ef.θ)
    return du * cosθ + dv * sinθ, -du * sinθ + dv * cosθ
end

function _from_ellipse_frame(ue, ve, ef::EllipseFit)
    cosθ = cos(ef.θ)
    sinθ = sin(ef.θ)
    du = ue * cosθ - ve * sinθ
    dv = ue * sinθ + ve * cosθ
    return du + ef.center_local[1], dv + ef.center_local[2]
end

"""
    _rotate_single_spark(point_3d, ef::EllipseFit, Δφ, cosφ, sinφ) -> Vector{Float64}

Full pipeline for a single spark:
  3D → tangent 2D → ellipse frame → normalize → rotate → unscale → tangent 2D → 3D on sphere
"""
function _rotate_single_spark(point_3d, ef::EllipseFit, cosφ, sinφ)
    # Project to tangent plane
    u, v = _to_local_2d(point_3d, ef)

    # To ellipse-aligned frame
    ue, ve = _to_ellipse_frame(u, v, ef)

    # Normalize to unit circle
    un = ue / ef.a
    vn = ve / ef.b

    # Rotate on circle
    ur = un * cosφ - vn * sinφ
    vr = un * sinφ + vn * cosφ

    # Unscale back to ellipse
    ue2 = ur * ef.a
    ve2 = vr * ef.b

    # Back to tangent plane frame
    u2, v2 = _from_ellipse_frame(ue2, ve2, ef)

    # Reproject onto sphere
    return _to_3d_on_sphere(u2, v2, ef)
end

"""
    rotate_sparks!(positions::Matrix{Float64}, ef::EllipseFit, Δφ::Float64)

In-place solid-body rotation of sparks on an elliptical polar cap.
positions: 3×N matrix, modified in place.
ef: pre-computed EllipseFit from fit_ellipse().
Δφ: rotation angle in radians.
"""
function rotate_sparks!(positions::Matrix{Float64}, ef::EllipseFit, Δφ::Float64)
    cosφ = cos(Δφ)
    sinφ = sin(Δφ)
    for i in axes(positions, 2)
        positions[:, i] = _rotate_single_spark(view(positions, :, i), ef, cosφ, sinφ)
    end
    return positions
end

"""
    rotate_sparks(positions::Matrix{Float64}, ef::EllipseFit, Δφ::Float64) -> Matrix{Float64}

Out-of-place version. Returns new matrix, original unchanged.
"""
function rotate_sparks(positions::Matrix{Float64}, ef::EllipseFit, Δφ::Float64)
    return rotate_sparks!(copy(positions), ef, Δφ)
end

"""
    rotate_sparks(positions::Matrix{Float64}, boundary::Matrix{Float64},
                  R::Float64, Δφ::Float64) -> Matrix{Float64}

Convenience: fit ellipse + rotate in one call.
Use the EllipseFit version for animation loops.
"""
function rotate_sparks(positions::Matrix{Float64}, boundary::Matrix{Float64},
                       R::Float64, Δφ::Float64)
    ef = fit_ellipse(boundary, R)
    return rotate_sparks(positions, ef, Δφ)
end

### added for variable conversion
function fit_ellipse(bpts::Vector{Vector{Float64}}, R::Real)
    mat = permutedims(hcat(bpts...))
    return fit_ellipse(mat, R)
end

function rotate_sparks!(positions::Vector{<:AbstractVector{<:Real}}, ef::EllipseFit, Δφ::Real)
    cosφ = cos(Δφ)
    sinφ = sin(Δφ)
    for i in eachindex(positions)
        positions[i] = _rotate_single_spark(positions[i], ef, cosφ, sinφ)
    end
    return positions
end

end # module