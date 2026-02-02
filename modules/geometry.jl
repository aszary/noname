"""
    Geometry

Module for converting rotational phase φ to emission point coordinates (θ, ψ)
in the magnetic axis reference frame.

Based on equations from Lorimer & Kramer (Handbook of Pulsar Astronomy, 2004).

# Coordinate System
- φ: rotational phase (longitude in the observer's frame)
- ρ: cone opening angle (angular distance from magnetic axis at the stellar surface)
- θ: colatitude from magnetic axis at emission height
- ψ: azimuthal angle around magnetic axis

# Key Equations Referenced
- Eq. (3.23): Relates θ to emission height via sin²θ = r/r_LC
- Eq. (3.27): Spherical triangle relation for beam geometry
- Eq. (3.28): Dipole field line geometry relating ρ and θ
"""
module Geometry

    using PhysicalConstants.CODATA2018 # in SI units


    # Speed of light [m/s]
    #const c = 299792458.0
    const c = SpeedOfLightInVacuum.val # no units hereafter


    """
        θ_from_ρ(ρ)

    Convert cone opening angle ρ to colatitude θ at emission height.

    For a dipolar magnetic field, the relationship between the footpoint angle ρ
    (at the stellar surface) and the colatitude θ (at emission height r_em)
    follows from the dipole field line equation r = r_0 sin²θ.

    From Lorimer & Kramer Eq. (3.28):
        tan(θ) = -3/(2tan(ρ)) + √[2 + (3/(2tan(ρ)))²]

    # Arguments
    - `ρ::Real`: cone opening angle at stellar surface [radians]

    # Returns
    - `θ::Float64`: colatitude from magnetic axis at emission height [radians]

    # Notes
    The solution corresponds to the smaller root, appropriate for emission
    within the open field line region.
    """
    function θ_from_ρ(ρ)
        if ρ ≈ 0
            return 0.0
        end
        tan_ρ = tan(ρ)
        x = 3 / (2 * tan_ρ)
        tan_θ = -x + sqrt(2 + x^2)
        return atan(tan_θ)
    end

    """
        ρ_from_θ(θ)

    Convert colatitude θ at emission height to cone opening angle ρ at stellar surface.

    Inverse of Eq. (3.28) from Lorimer & Kramer. Derived by algebraic inversion:
        tan(ρ) = 3·tan(θ) / (2 - tan²(θ))

    # Arguments
    - `θ::Real`: colatitude from magnetic axis at emission height [radians]

    # Returns
    - `ρ::Float64`: cone opening angle at stellar surface [radians]

    # Notes
    Valid for θ < arctan(√2) ≈ 54.7° where the denominator remains positive.
    """
    function ρ_from_θ(θ)
        if θ ≈ 0
            return 0.0
        end
        # Derivation from Eq. (3.28):
        # Let y = tan(θ), x = 3/(2·tan(ρ))
        # From: y = -x + √(2 + x²)
        # Rearranging: (y + x)² = 2 + x²
        # Expanding: y² + 2xy = 2
        # Solving for x: x = (2 - y²)/(2y)
        # Therefore: tan(ρ) = 3/(2x) = 3y/(2 - y²)
        
        y = tan(θ)
        tan_ρ = 3 * y / (2 - y^2)
        return atan(tan_ρ)
    end

    """
        θ_em_max(r_em, P)

    Calculate the maximum colatitude θ for the last open field line at emission height r_em.

    From Lorimer & Kramer Eq. (3.23), the last open field line satisfies:
        sin²θ = r_em / r_LC

    where r_LC = cP/(2π) is the light cylinder radius.

    # Arguments
    - `r_em::Real`: emission height above neutron star center [m]
    - `P::Real`: pulsar rotation period [s]

    # Returns
    - `θ_max::Float64`: maximum colatitude for open field lines [radians]

    # Throws
    - `ErrorException`: if r_em ≥ r_LC (emission height beyond light cylinder)

    # Example
    ```julia
    # For typical emission height of 500 km and 1 second period:
    θ_max = θ_em_max(500e3, 1.0)  # ≈ 0.032 rad ≈ 1.8°
    ```
    """
    function θ_em_max(r_em, P)
        r_LC = c * P / (2π)
        sin²θ = r_em / r_LC
        if sin²θ >= 1
            error("r_em ≥ r_LC: emission height beyond light cylinder")
        end
        return asin(sqrt(sin²θ))
    end

    """
        ρ_max(r_em, P)

    Calculate the maximum cone opening angle ρ for the last open field lines.

    Combines `θ_em_max` and `ρ_from_θ` to give the angular radius of the
    polar cap as seen from the stellar surface.

    # Arguments
    - `r_em::Real`: emission height [m]
    - `P::Real`: pulsar rotation period [s]

    # Returns
    - `ρ_max::Float64`: maximum cone opening angle [radians]

    # Notes
    For a typical pulsar with P = 1 s and r_em = 500 km, ρ_max ≈ 2.7°.
    The polar cap radius at the surface is approximately R_NS · sin(ρ_max).
    """
    function ρ_max(r_em, P)
        θ_max = θ_em_max(r_em, P)
        return ρ_from_θ(θ_max)
    end

    """
        φ_max(α, β, r_em, P)

    Calculate the maximum rotational phase φ where the line of sight
    intersects the open field line region (half-width of the pulse profile).

    From Lorimer & Kramer Eq. (3.27), using spherical trigonometry:
        cos(ρ) = cos(α)·cos(ζ) + sin(α)·sin(ζ)·cos(φ)

    where ζ = α + β is the angle between rotation axis and line of sight.

    # Arguments
    - `α::Real`: magnetic inclination angle (rotation axis to magnetic axis) [radians]
    - `β::Real`: impact parameter (closest approach of line of sight to magnetic axis) [radians]
    - `r_em::Real`: emission height [m]
    - `P::Real`: pulsar rotation period [s]

    # Returns
    - `φ_max::Union{Float64, Nothing}`: maximum phase [radians], or `nothing` if beam is invisible

    # Notes
    - Returns `nothing` if cos(φ_max) > 1 (line of sight never enters beam)
    - Returns `π` if cos(φ_max) < -1 (line of sight always inside beam, full rotation visible)
    - The pulse width W = 2·φ_max

    # Example
    ```julia
    α = deg2rad(45.0)  # 45° inclination
    β = deg2rad(5.0)   # 5° impact parameter
    φ_m = φ_max(α, β, 500e3, 1.0)
    pulse_width_deg = 2 * rad2deg(φ_m)
    ```
    """
    function φ_max(α, β, r_em, P)
        ρ_m = ρ_max(r_em, P)
        ζ = α + β  # angle between rotation axis and line of sight
        
        cos_φ_max = (cos(ρ_m) - cos(α) * cos(ζ)) / (sin(α) * sin(ζ))
        
        if cos_φ_max > 1
            # Line of sight never intersects the beam
            return nothing
        elseif cos_φ_max < -1
            # Line of sight always inside beam (interpulsar or aligned rotator)
            return π
        end
        
        return acos(cos_φ_max)
    end

    """
        ρ_from_φ(φ, α, β)

    Calculate cone opening angle ρ for a given rotational phase φ.

    Inverse application of Lorimer & Kramer Eq. (3.27). Given the observer's
    rotational phase, this computes which cone (defined by ρ) the line of sight
    is crossing.

    # Arguments
    - `φ::Real`: rotational phase [radians], where φ=0 is closest approach
    - `α::Real`: magnetic inclination angle [radians]
    - `β::Real`: impact parameter [radians]

    # Returns
    - `ρ::Float64`: cone opening angle [radians]

    # Notes
    The result is clamped to [-1, 1] before taking acos to handle
    numerical edge cases near the beam boundaries.
    """
    function ρ_from_φ(φ, α, β)
        ζ = α + β
        cos_ρ = cos(α) * cos(ζ) + sin(α) * sin(ζ) * cos(φ)
        return acos(clamp(cos_ρ, -1, 1))
    end

    """
        emission_points_from_phase(φ_array, α, β, r_em, P)

    Main function: Convert array of rotational phases φ to emission point
    coordinates (θ, ψ) in the magnetic axis reference frame.

    This is the primary interface for mapping observed pulse phases to
    physical locations in the magnetosphere.

    # Arguments
    - `φ_array::AbstractVector`: array of rotational phases [radians]
    - `α::Real`: magnetic inclination angle [radians]
    - `β::Real`: impact parameter [radians]
    - `r_em::Real`: emission height [m]
    - `P::Real`: pulsar rotation period [s]

    # Returns
    - `θ_array::Vector{Float64}`: colatitudes from magnetic axis [radians]
    - `ψ_array::Vector{Float64}`: azimuthal angles around magnetic axis [radians]

    # Notes
    The azimuthal angle ψ equals the rotational phase φ due to the geometry
    of the spherical triangle formed by the rotation axis, magnetic axis,
    and emission point. This assumes the conventional definition where ψ=0
    corresponds to the fiducial plane containing both axes.

    # Example
    ```julia
    α = deg2rad(30.0)
    β = deg2rad(3.0)
    r_em = 300e3  # 300 km
    P = 0.5       # 500 ms

    φ_array = generate_uniform_phase_array(100, α, β, r_em, P)
    θ_array, ψ_array = emission_points_from_phase(φ_array, α, β, r_em, P)
    ```
    """
    function emission_points_from_phase(φ_array, α, β, r_em, P)
        n = length(φ_array)
        θ_array = zeros(n)
        ψ_array = zeros(n)
        
        for i in 1:n
            φ = φ_array[i]
            ρ = ρ_from_φ(φ, α, β)
            θ_array[i] = θ_from_ρ(ρ)
            # Azimuthal angle equals rotational phase from spherical triangle geometry
            ψ_array[i] = φ
        end
        
        return θ_array, ψ_array
    end

    """
        generate_uniform_phase_array(n_points, α, β, r_em, P)

    Generate uniformly spaced array of phases within the open field line region.

    Creates a symmetric array of phases from -φ_max to +φ_max, suitable for
    sampling the entire visible portion of the emission beam.

    # Arguments
    - `n_points::Integer`: number of phase points to generate
    - `α::Real`: magnetic inclination angle [radians]
    - `β::Real`: impact parameter [radians]
    - `r_em::Real`: emission height [m]
    - `P::Real`: pulsar rotation period [s]

    # Returns
    - `φ_array::StepRangeLen`: array of phases [radians]

    # Throws
    - `ErrorException`: if line of sight doesn't intersect open field line region

    # Example
    ```julia
    # Generate 256 phase bins across the pulse window
    φ = generate_uniform_phase_array(256, deg2rad(45), deg2rad(5), 500e3, 1.0)
    ```
    """
    function generate_uniform_phase_array(n_points, α, β, r_em, P)
        φ_m = φ_max(α, β, r_em, P)
        
        if φ_m === nothing
            error("Line of sight does not intersect open field line region")
        end
        
        return range(-φ_m, φ_m, length=n_points)
    end

end # module Geometry
