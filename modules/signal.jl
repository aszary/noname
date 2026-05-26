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

    function generate_signal_radii(psr; noise_level=0.1, v_scale=0.3)
        los_points = []
        for line in psr.los_lines
            push!(los_points, [line[1][end], line[2][end], line[3][end]])
        end
        signal_number = size(psr.sparks_locations)[1]
        bin_number = size(los_points)[1]
        psr.signal = zeros(signal_number, bin_number)
        for (i, p) in enumerate(los_points)
            for (j, sparks) in enumerate(psr.sparks_locations)
                for (k, s) in enumerate(sparks)
                    if isnothing(psr.spark_radii)
                        sigma = psr.spark_radius / 3.72 # 2.355->FWHM, 3.03->1%, 3.72->0.1%
                    else
                        sigma = psr.spark_radii[j][k] / 3.72
                    end
                    dist = norm(p - s)
                    psr.signal[j, i] += exp(-dist^2 / (2 * sigma^2))
                end
            end
        end
        signal_clean = copy(psr.signal)  # save noise-free for polarization
        noise = noise_level * randn(size(psr.signal))
        psr.signal .+= noise

        # position angle calculation
        psr.pa = zeros(bin_number)
        # psr.line_of_sight stores the A/R-corrected emission point (ψ shifted by -Δψ).
        # For the LOS direction used in PA, we need the observed direction (before A/R shift).
        Δψ_ar = 2 * psr.r_em / Functions.rlc(psr.p)
        for i in 1:length(psr.longitudes)
            # Get the field line assigned to this phase bin
            line = psr.los_lines[i]

            # Point 1 is the emission point (topmost), Point 2 is slightly lower on the same line
            p1 = [line[1][1], line[2][1], line[3][1]]
            p2 = [line[1][2], line[2][2], line[3][2]]

            # Local magnetic field vector is the direction between these two consecutive line points
            B_local = p1 .- p2

            # Undo the A/R azimuthal shift to recover the true observed LOS direction
            los_ar = psr.line_of_sight[i]
            los_current = [los_ar[1]*cos(Δψ_ar) - los_ar[2]*sin(Δψ_ar),
                           los_ar[1]*sin(Δψ_ar) + los_ar[2]*cos(Δψ_ar),
                           los_ar[3]]
            rot_vec = Functions.spherical2cartesian(psr.rotation_axis)

            # Calculate numerical PA in radians, then convert to degrees
            psr.pa[i] = rad2deg(calculate_numerical_pa(B_local, los_current, rot_vec))
        end

        # Full Stokes: Q, U, V per pulse
        # V(φ) ∝ dI/dφ per pulse; L² + V² = I² (100% total polarization)
        psr.stokes_q = zeros(signal_number, bin_number)
        psr.stokes_u = zeros(signal_number, bin_number)
        psr.stokes_v = zeros(signal_number, bin_number)
        for j in 1:signal_number
            pulse = @view signal_clean[j, :]
            # central differences for dI/dφ
            dI = zeros(bin_number)
            for k in 2:bin_number-1
                dI[k] = (pulse[k+1] - pulse[k-1]) / 2
            end
            dI[1]   = pulse[2] - pulse[1]
            dI[end] = pulse[end] - pulse[end-1]
            # scale V so max|V| = v_scale * max|I|
            max_dI = maximum(abs.(dI))
            max_I  = maximum(abs.(pulse))
            V = (max_dI > 0 && max_I > 0) ? v_scale * max_I * dI / max_dI : zeros(bin_number)
            # clamp so |V| ≤ |I| everywhere, then L from remainder
            V = clamp.(V, -abs.(pulse), abs.(pulse))
            L = sqrt.(max.(0.0, pulse .^ 2 .- V .^ 2))
            psr.stokes_v[j, :] = V
            psr.stokes_q[j, :] = L .* cos.(2 .* deg2rad.(psr.pa))
            psr.stokes_u[j, :] = L .* sin.(2 .* deg2rad.(psr.pa))
        end

    end

    function generate_pulses(psr)

        signal_number, bin_number = size(psr.signal)
        psr.pulses = zeros(psr.npulse, bin_number)

        for i in 1:min(signal_number, psr.npulse)
            psr.pulses[i, :] = psr.signal[i, :]
        end

    end


   """
    calculate_numerical_pa(B_vec, los_vec, rot_vec)

    Calculates the Polarization Angle (PA) directly from the local magnetic field.
    Projects the B vector from the emission point onto the observer's Plane of the Sky.
    """
    function calculate_numerical_pa(B_vec, los_vec, rot_vec)
        # 1. Ensure all vectors have a length of 1 (normalize)
        n = normalize(los_vec)       # Line of Sight vector (towards the telescope)
        omega = normalize(rot_vec)   # Rotation Axis vector (Red axis)
        B = normalize(B_vec)         # Local magnetic field vector at the emission point
        
        # 2. Build a 2D coordinate system on the Plane of the Sky (as seen by the telescope)
        # Y-axis on the sky (North) is the projection of the rotation axis onto the viewing plane
        y_sky = omega .- dot(omega, n) .* n
        y_sky = normalize(y_sky)
        
        # X-axis on the sky (East) is perpendicular to the Line of Sight and North
        x_sky = cross(n, y_sky)
        
        # 3. Project our local magnetic field vector (B) onto this sky plane
        B_x = dot(B, x_sky)
        B_y = dot(B, y_sky)
        
        # 4. Calculate the angle using atan2 (yields a result from -pi to pi)
        pa = atan(B_x, B_y)
        
        # Normalize the angle to the [-pi/2, pi/2] range, standard for RVM
        return atan(tan(pa))
    end



"""
    pulse_width(α, β, ρ)

Oblicza szerokość profilu pulsara W na podstawie geometrii wiązki.

# Argumenty
- `α`: kąt inklinacji między osią rotacji a osią magnetyczną [radiany]
- `β`: parametr uderzenia - najbliższe przejście linii widzenia od osi magnetycznej [radiany]
- `ρ`: kąt otwarcia wiązki emisyjnej [radiany]

# Zwraca
- `W`: szerokość profilu w radianach (lub `NaN` jeśli wiązka nie jest przecięta)

# Uwagi
Używa równania: sin²(W/4) = [sin²(ρ/2) - sin²(β/2)] / [sin(α)·sin(α+β)]
Warunek |β| ≤ ρ musi być spełniony, aby wiązka była widoczna.
"""
function pulse_width(α, β, ρ)
    # Sprawdź czy wiązka jest przecięta
    if abs(β) > ρ
        return NaN
    end

    # align rotator
    if abs(α) < 1e-10
        return 2π  # 360°
    end
    
    numerator = sin(ρ/2)^2 - sin(β/2)^2
    denominator = sin(α) * sin(α + β)
    
    # Sprawdź poprawność wartości
    if denominator ≤ 0 || numerator < 0
        return NaN
    end
    
    sin2_W4 = numerator / denominator
    
    # sin² nie może przekroczyć 1
    if sin2_W4 > 1
        return 2π
    end
    
    W = 4 * asin(sqrt(sin2_W4))
    return W
end

# Wersja z kątami w stopniach
function pulse_width_deg(α_deg, β_deg, ρ_deg)
    W_rad = pulse_width(deg2rad(α_deg), deg2rad(β_deg), deg2rad(ρ_deg))
    return rad2deg(W_rad)
end


"""
    rho_from_theta(θ)

Numeric solution for ρ from θ 
"""
function rho_from_theta(θ; tol=1e-12, max_iter=50)

    ρ = 1.5 * θ
    for _ in 1:max_iter
        x = 3 / (2 * tan(ρ))
        θ_calc = atan(-x + sqrt(2 + x^2))
        Δ = θ - θ_calc
        
        if abs(Δ) < tol
            return ρ
        end
        
        ρ += 0.5 * Δ * 1.5
    end
    
    return ρ

end


end # module end