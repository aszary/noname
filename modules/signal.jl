 module Signal
    using LinearAlgebra
    using Statistics

    include("functions.jl")
    import ..Geometry
    

    function generate_signal(psr; noise_level=0.0)
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

    function generate_signal_radii(psr; noise_level=0.1)
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
        noise = noise_level * randn(size(psr.signal))
        psr.signal .+= noise

        # position angle calculation
        psr.pa = zeros(bin_number)
        for i in 1:length(psr.longitudes)
            # Get the field line assigned to this phase bin
            line = psr.los_lines[i]

            # Point 1 is the emission point (topmost), Point 2 is slightly lower on the same line
            p1 = [line[1][1], line[2][1], line[3][1]]
            p2 = [line[1][2], line[2][2], line[3][2]]

            # Local magnetic field vector is the direction between these two consecutive line points
            B_local = p1 .- p2

            # Line of sight and rotation axis vectors for this bin
            los_current = psr.line_of_sight[i]
            rot_vec = Functions.spherical2cartesian(psr.rotation_axis)

            # Calculate numerical PA in radians, then convert to degrees
            pa_rad = calculate_numerical_pa(B_local, los_current, rot_vec)
            psr.pa[i] = rad2deg(pa_rad)
        end

    end

   """
        generate_signal_solid_body(psr; noise_level=0.1)

    Generates a signal using discrete frozen sparks. It preserves the classic 
    Gaussian peaks but deforms each spark into an ellipse proportional to the 
    actual simulated polar cap boundary (psr.pc).
    """
    function generate_signal_solid_body(psr; noise_level=0.0)
        los_points = []
        for line in psr.los_lines
            push!(los_points, [line[1][end], line[2][end], line[3][end]])
        end
        signal_number = size(psr.sparks_locations)[1]
        bin_number = length(los_points)
        psr.signal = zeros(signal_number, bin_number)

        #Calculate the geometric center of the polar cap based on its true boundaries.
        pc_x, pc_y, pc_z = psr.pc[1], psr.pc[2], psr.pc[3]
        cx, cy, cz = sum(pc_x)/length(pc_x), sum(pc_y)/length(pc_y), sum(pc_z)/length(pc_z)
        center = [cx, cy, cz]

        # Push the center exactly onto the spherical surface of the star
        center = (center ./ norm(center)) .* psr.r 
        
        # Group boundary points into 3D vectors for easier math
        pc_boundary = [[pc_x[i], pc_y[i], pc_z[i]] for i in 1:length(pc_x)]

        # Calculate the average internal radius of this specific polar cap.
        # This makes the function completely immune to scale mismatches between theoretical dipole radii and the actual simulated NSField dimensions.
        radii_list = [norm(b .- center) for b in pc_boundary]
        avg_pc_radius = sum(radii_list) / length(radii_list)

        for (i, p) in enumerate(los_points)
            for (j, sparks) in enumerate(psr.sparks_locations)
                for (k, s) in enumerate(sparks)
                    # Base standard deviation (sigma) for a perfectly circular spark# Baza - idealne koło
                    if isnothing(psr.spark_radii)
                        base_sigma = psr.spark_radius / 3.72 
                    else
                        base_sigma = psr.spark_radii[j][k] / 3.72
                    end
                    
                    # Vector pointing from the spark center directly to the telescope's current position
                    vec_ps = p .- s
                    dist_ps = norm(vec_ps)
                    
                    # Default stretch factor (assumes a perfect circle)
                    ellipse_stretch = 1.0
                    
                    # If the telescope isn't exactly dead-center on the spark, calculate the shape distortion
                    if dist_ps > 1e-6
                        dir_ps = vec_ps ./ dist_ps
                        best_cos = -1.0
                        dist_edge = avg_pc_radius # Fallback to average radius
                        
                        # Find the exact distance to the polar cap boundary in the specific direction the telescope is currently slicing through.
                        for b in pc_boundary
                            vec_cb = b .- center
                            norm_cb = norm(vec_cb)
                            if norm_cb > 1e-6
                                dir_cb = vec_cb ./ norm_cb
                                cos_angle = dot(dir_ps, dir_cb)
                                # Maximize the dot product to find the boundary point lying perfectly on this ray
                                if cos_angle > best_cos
                                    best_cos = cos_angle
                                    dist_edge = norm_cb # This is the length of the ellipse's "arm" in this direction
                                end
                            end
                        end
                        
                        # Calculate stretch: Length of the ellipse arm in this direction / Average polar cap radius
                        # e.g., if the ellipse is stretched by 20% in this direction, the stretch factor = 1.2
                        ellipse_stretch = dist_edge / avg_pc_radius
                    end
                    
                    # Apply the polar cap's elliptical stretch directly to the individual spark's width
                    local_sigma = base_sigma * ellipse_stretch
                    
                    # Calculate the classic Gaussian intensity using the newly deformed local_sigma
                    psr.signal[j, i] += exp(-dist_ps^2 / (2 * local_sigma^2))
                end
            end
        end

        # Add noise
        noise = noise_level * randn(size(psr.signal))
        psr.signal .+= noise

        #  Position angle calculation
        psr.pa = zeros(bin_number)
        for i in 1:length(psr.longitudes)
            line = psr.los_lines[i]
            p1 = [line[1][1], line[2][1], line[3][1]]
            p2 = [line[1][2], line[2][2], line[3][2]]
            B_local = p1 .- p2
            
            los_current = psr.line_of_sight[i]
            rot_vec = Functions.spherical2cartesian(psr.rotation_axis)
            
            pa_rad = calculate_numerical_pa(B_local, los_current, rot_vec)
            psr.pa[i] = rad2deg(pa_rad)
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

    
    #stare PA
    function generate_signal_with_polarization(psr; noise_level=0.1)
        los_points = []
        for line in psr.los_lines
            push!(los_points, [line[1][end], line[2][end], line[3][end]])
        end

        signal_number = size(psr.sparks_locations)[1]
        bin_number = size(los_points)[1]
        
        # --- NEW: Calculate Polarization Angle (PA) ---
        # We calculate one PA value for each bin in the profile
        psr.pa = zeros(bin_number)
        alpha_rad = deg2rad(psr.alpha)
        beta_rad = deg2rad(psr.beta)

        for i in 1:length(psr.longitudes)
            # Pobieramy linię pola przypisaną do tej fazy
            line = psr.los_lines[i] 
            
            # Punkt 1 to punkt emisji (najwyższy), Punkt 2 to punkt minimalnie niżej na tej samej linii
            p1 = [line[1][1], line[2][1], line[3][1]]
            p2 = [line[1][2], line[2][2], line[3][2]]
           
            
            # Wektor lokalnego pola magnetycznego to po prostu kierunek linii między tymi punktami!
            B_local = p1 .- p2
            
            # Wektor widzenia i rotacji
            los_current = psr.line_of_sight[i]
            rot_vec = Functions.spherical2cartesian(psr.rotation_axis)
            
            # Obliczamy numeryczne PA w radianach, a potem na stopnie
            pa_rad = calculate_numerical_pa(B_local, los_current, rot_vec)
            psr.pa[i] = rad2deg(pa_rad)
        end
        # -----------------------------------------------

        psr.signal = zeros(signal_number, bin_number)
        
        for (i, p) in enumerate(los_points)
            for (j, sparks) in enumerate(psr.sparks_locations)
                for (k, s) in enumerate(sparks)
                    if isnothing(psr.spark_radii)
                        sigma = psr.spark_radius / 3.72 
                    else
                        sigma = psr.spark_radii[j][k] / 3.72
                    end
                    dist = norm(p - s)
                    # Adding the spark's intensity to the specific bin[cite: 8]
                    psr.signal[j, i] += exp(-dist^2 / (2 * sigma^2))
                end
            end
        end
        
        noise = noise_level * randn(size(psr.signal))
        psr.signal .+= noise
    end

    function generate_pulses(psr)
        # use skip_steps in simulate_sparks to have single pulses

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