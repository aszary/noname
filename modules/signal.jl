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
        #psr.signal .+= noise

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
    Calculates pulse longitude based on coordinates of boundry open field lines at the emission height
    HACK!
    """
    function calculate_longitudes(psr)
        bins = size(psr.pulses)[2]
        #psr.longitudes = Vector{Float64}(undef, bins)
        alpha = deg2rad(psr.alpha)
        beta = deg2rad(psr.beta)

        # get boundry emission points
        x1 = psr.los_lines[1][1][1]
        y1 = psr.los_lines[1][2][1]
        z1 = psr.los_lines[1][3][1]

        x2 = psr.los_lines[end][1][1]
        y2 = psr.los_lines[end][2][1]
        z2 = psr.los_lines[end][3][1]

        l1 = point_to_longitude([x1, y1, z1], alpha, beta; exact=true)
        l2 = point_to_longitude([x2, y2, z2], alpha, beta; exact=true)

        l1deg = rad2deg(l1)
        l2deg = rad2deg(l2)

        psr.longitudes = range(l1deg, l2deg, length=bins)

    end


    """
    Calculates pulse longitude based on coordinates of open field lines at the emission height 
    GAP IN THE MIDDLE?! but this seems to be real
    """
    function calculate_longitudes_real(psr)
        psr.longitudes = Vector{Float64}(undef, size(psr.pulses)[2])
        alpha = deg2rad(psr.alpha)
        beta = deg2rad(psr.beta)

        # get emission points
        xs = []
        ys = []
        zs = []
        npoints = size(psr.los_lines)[1]
        for p in 1:npoints
            x = psr.los_lines[p][1][1]
            y = psr.los_lines[p][2][1]
            z = psr.los_lines[p][3][1]
            push!(xs, x)
            push!(ys, y)
            push!(zs, z)
        end

        # calculate longitude
        for (i, x) in enumerate(xs)
            ll = point_to_longitude([xs[i], ys[i], zs[i]], alpha, beta; exact=true)
            lldeg = rad2deg(ll)
            psr.longitudes[i] = lldeg
            #println("$i long = $lldeg [deg.] $(xs[i]) $(ys[i]) $(zs[i]) ")
        end
    end


    """
        point_to_longitude(point, α, β; exact=true)

    Compute longitude for a 3D emission point given in magnetic frame coordinates.

    # Arguments
    - `point`: (x, y, z) in magnetic frame where z is along magnetic axis
    - `α`: inclination angle in radians
    - `β`: impact parameter in radians

    # Returns
    - `φ`: longitude in radians
    """
    function point_to_longitude(point, α, β; exact=true)
        x, y, z = point
        r = sqrt(x^2 + y^2 + z^2)
        θ = acos(z / r)  # Colatitude from magnetic axis
        ϕ_mag = atan(y, x)  # Azimuth in magnetic frame
        
        ρ = theta_to_rho(θ; exact=exact)
        φ_abs = rho_to_longitude(ρ, α, β)
        
        isnan(φ_abs) && return NaN
        
        return sign(ϕ_mag) * φ_abs
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

end # module end