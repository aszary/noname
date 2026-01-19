 module Signal
    using LinearAlgebra

    include("functions.jl")

    function generate_signal(psr; noise_level=0.1)
        # line of sight points at the polar cap
        los_points = []

        for line in psr.los_lines
            push!(los_points, [line[1][end], line[2][end], line[3][end]])
        end
        # pulsar radio signal
        signal_number = size(psr.sparks_locations)[1]
        bin_number = size(los_points)[1]
        psr.signal = zeros(signal_number, bin_number)
        sigma = psr.spark_radius / 3.72 # 2.355->FWHM, 3.03->1%, 3.72->0.1%   
        for (i, p) in enumerate(los_points)
            for (j, sparks) in enumerate(psr.sparks_locations)
                for s in sparks
                    dist = norm(p - s)
                    psr.signal[j, i] += exp(-dist^2 / (2 * sigma^2)) 
                end
            end
        end

        # adding gaussian noise
        noise = noise_level * randn(size(psr.signal))
        psr.signal .+= noise

    end

    function generate_pulses(psr; pulse_max=100)
        # use skip_steps in simulate_sparks to have single pulses

        # TODO generate full phase, change to longitude (calculate rho!)

        signal_number, bin_number = size(psr.signal)
        psr.pulses = zeros(pulse_max, bin_number)

        for i in 1:signal_number
            psr.pulses[i, :] = psr.signal[i, :]
            if i == pulse_max
                break
            end
        end

    end


    """
    Calculates pulse longitude based on coordinates of open field lines at the emission height 
    """
    function calculate_longitude(psr)
        alpha = deg2rad(psr.alpha)
        beta = deg2rad(psr.beta)

        # get emission points
        rs = []
        thetas = []
        phis = []
        npoints = size(psr.los_lines)[1]
        for p in 1:npoints
            x = psr.los_lines[p][1][1]
            y = psr.los_lines[p][2][1]
            z = psr.los_lines[p][3][1]
            r, theta, phi = Functions.cartesian2spherical([x, y, z])
            #println("p x=$x y=$y z=$z")
            #println("p r=$r theta=$theta phi=$phi")
            push!(rs, r)
            push!(thetas, theta)
            push!(phis, phi)
        end

        # calculate longitude
        for (i, r) in enumerate(rs)
            theta = thetas[i]
            phi = phis[i]
            #println("$i $r $theta")
            
            lm, lp = coords_to_longitude(r, theta, alpha, beta; exact=true)
            lmdeg = rad2deg(lm)
            lpdeg = rad2deg(lp)
            println("$i $r theta=$theta phi=$phi lm = $lmdeg [deg.] lp = $lpdeg [deg.]")

        end

        #println(size(psr.los_lines[1][1]))
        #for line in psr.los_lines


    end
    
"""
    theta_to_rho(θ; exact=true)

Compute the emission cone opening angle ρ from the polar coordinate θ
of a dipolar field line.

# Arguments
- `θ`: polar coordinate in radians
- `exact`: if true, use exact formula (3.28); if false, use approximation θ ≈ 2ρ/3

# Returns
- `ρ`: cone opening angle in radians
"""
function theta_to_rho(θ; exact=true)
    if !exact || θ < deg2rad(20)
        # Small angle approximation
        return 1.5 * θ
    else
        # Exact solution of equation (3.28)
        # tan(θ) = -3/(2tan(ρ)) + √(2 + (3/(2tan(ρ)))²)
        # Solve numerically using Newton-Raphson
        return _solve_theta_to_rho(θ)
    end
end

function _solve_theta_to_rho(θ; tol=1e-10, maxiter=50)
    # Newton-Raphson solver for the exact θ-ρ relation
    ρ = 1.5 * θ  # Initial guess from approximation
    
    for _ in 1:maxiter
        tanρ = tan(ρ)
        x = 3 / (2 * tanρ)
        f = tan(θ) - (-x + sqrt(2 + x^2))
        
        # Derivative df/dρ
        sec2ρ = 1 / cos(ρ)^2
        dx_dρ = -3 * sec2ρ / (2 * tanρ^2)
        df_dρ = -dx_dρ * (1 - x / sqrt(2 + x^2))
        
        δ = f / df_dρ
        ρ -= δ
        
        abs(δ) < tol && break
    end
    
    return ρ
end

"""
    rho_to_theta(ρ; exact=true)

Compute the polar coordinate θ from the cone opening angle ρ.
Inverse of theta_to_rho.

# Arguments
- `ρ`: cone opening angle in radians
- `exact`: if true, use exact formula; if false, use approximation

# Returns
- `θ`: polar coordinate in radians
"""
function rho_to_theta(ρ; exact=true)
    if !exact || ρ < deg2rad(30)
        # Small angle approximation
        return 2ρ / 3
    else
        # Exact formula from equation (3.28)
        tanρ = tan(ρ)
        return atan(-3/(2*tanρ) + sqrt(2 + (3/(2*tanρ))^2))
    end
end

"""
    rho_to_longitude(ρ, α, β)

Compute the rotational longitude (phase) φ for a given cone opening angle ρ.

Uses equation (3.27):
    cos(ρ) = cos(α)cos(α+β) + sin(α)sin(α+β)cos(W/2)

# Arguments
- `ρ`: cone opening angle in radians
- `α`: inclination angle (between rotation and magnetic axes) in radians
- `β`: impact parameter in radians

# Returns
- `φ`: longitude in radians (half pulse width W/2)
"""
function rho_to_longitude(ρ, α, β)
    numerator = cos(ρ) - cos(α) * cos(α + β)
    denominator = sin(α) * sin(α + β)
    
    # Check if line of sight intersects the cone
    cos_half_W = numerator / denominator
    
    if abs(cos_half_W) > 1
        return NaN  # Cone not intersected by line of sight
    end
    
    return acos(cos_half_W)
end


"""
    coords_to_longitude(r, θ, α, β; exact=true)

Main function: compute signal longitude from emission point coordinates (r, θ)
on a dipolar field line.

# Arguments
- `r`: distance from stellar center (arbitrary units, not used in calculation)
- `θ`: polar coordinate in radians
- `α`: inclination angle in radians
- `β`: impact parameter in radians
- `exact`: whether to use exact formulas

# Returns
- `(φ_leading, φ_trailing)`: tuple of longitudes for both sides of the profile
"""
function coords_to_longitude(r, θ, α, β; exact=true)
    ρ = theta_to_rho(θ; exact=exact)
    φ = rho_to_longitude(ρ, α, β)
    
    if isnan(φ)
        return (NaN, NaN)
    end
    
    return (-φ, φ)  # Leading edge (φ < 0) and trailing edge (φ > 0)
end

"""
    longitude_to_coords(φ, α, β, r_em; exact=true)

Inverse function: compute coordinates (r, θ) for a given longitude.

# Arguments
- `φ`: longitude in radians
- `α`: inclination angle in radians
- `β`: impact parameter in radians
- `r_em`: emission height (radius)
- `exact`: whether to use exact formulas

# Returns
- `(r, θ)`: emission point coordinates
"""
function longitude_to_coords(φ, α, β, r_em; exact=true)
    # From equation (3.27), solve for ρ
    cos_rho = cos(α) * cos(α + β) + sin(α) * sin(α + β) * cos(φ)
    ρ = acos(clamp(cos_rho, -1, 1))
    
    θ = rho_to_theta(ρ; exact=exact)
    
    return (r_em, θ)
end

"""
    pulse_width(ρ, α, β)

Compute the full pulse width W for given geometry.

# Arguments
- `ρ`: cone opening angle in radians
- `α`: inclination angle in radians
- `β`: impact parameter in radians

# Returns
- `W`: pulse width in radians
"""
function pulse_width(ρ, α, β)
    φ = rho_to_longitude(ρ, α, β)
    return isnan(φ) ? NaN : 2φ
end

"""
    emission_height_from_rho(ρ, P)

Estimate emission height from cone opening angle using equation (3.29).

# Arguments
- `ρ`: cone opening angle in radians
- `P`: pulsar period in seconds

# Returns
- `r_em`: emission height in km
"""
function emission_height_from_rho(ρ, P)
    # From equation (3.29): ρ = 1.24° × (r_em/10km)^(1/2) × (P/s)^(-1/2)
    # Solving for r_em:
    ρ_deg = rad2deg(ρ)
    return 10.0 * (ρ_deg / 1.24)^2 * P  # in km
end



end # module end