module Lines
    using LinearAlgebra
    import ..NSField

    include("field.jl")
    include("functions.jl")
    include("transformations.jl")
    include("signal.jl")
    include("geometry.jl")
    include("solid_body.jl")


    """
    calculate_polarcaps!(psr; phi_num=100)

    Calculates the polar cap boundaries by tracing the last open magnetic field lines 
    from height rmax down to the stellar surface (r = psr.r), taking into account 
    the influence of magnetic anomalies (NSField).
    """
    function calculate_polarcaps!(psr; phi_num=100)
        nf = psr.nsfield
        step = nf.rmax / nf.size

        # Opening angle at maximum height (where the field is dipole-dominated)
        z_rmax = nf.rmax / psr.r  
        theta_rmax = Functions.theta_max(z_rmax, psr)
        
        # Generate azimuthal points
        phis = range(0, 2*pi, length=phi_num+1)[1:phi_num]

        # Arrays for coordinates of the northern and southern hemispheres
        x_n, y_n, z_n = Float64[], Float64[], Float64[]
        x_s, y_s, z_s = Float64[], Float64[], Float64[]

        for phi in phis
            # ==========================================
            # NORTHERN HEMISPHERE
            # ==========================================
            pos_n = Functions.spherical2cartesian([nf.rmax, theta_rmax, phi])
            pos_sph_n = Functions.cartesian2spherical(pos_n)

            # Integrate the field line downwards to the surface
            st_n = [0.0, 0.0, 0.0]
            while pos_sph_n[1] > psr.r
                b_sph_n = NSField.BSph(nf, pos_sph_n[1]/psr.r, pos_sph_n[2], pos_sph_n[3])
                b_n = Functions.vec_spherical2cartesian(pos_sph_n, [b_sph_n[1], b_sph_n[2], b_sph_n[3]])
                st_n = -b_n / norm(b_n) * step
                pos_n += st_n
                pos_sph_n = Functions.cartesian2spherical(pos_n)
            end
            
            # Interpolate exactly to the surface r = psr.r (avoiding overshoot)
            r_last_n = norm(pos_n)
            pos_prev_n = pos_n - st_n
            r_prev_n = norm(pos_prev_n)
            t_n = (r_prev_n - psr.r) / (r_prev_n - r_last_n)
            surf_n = pos_prev_n + t_n * (pos_n - pos_prev_n)

            push!(x_n, surf_n[1])
            push!(y_n, surf_n[2])
            push!(z_n, surf_n[3])

            # ==========================================
            # SOUTHERN HEMISPHERE
            # ==========================================
            # Start from angle pi - theta_rmax
            pos_s = Functions.spherical2cartesian([nf.rmax, pi - theta_rmax, phi])
            pos_sph_s = Functions.cartesian2spherical(pos_s)

            st_s = [0.0, 0.0, 0.0]
            while pos_sph_s[1] > psr.r
                b_sph_s = NSField.BSph(nf, pos_sph_s[1]/psr.r, pos_sph_s[2], pos_sph_s[3])
                b_s = Functions.vec_spherical2cartesian(pos_sph_s, [b_sph_s[1], b_sph_s[2], b_sph_s[3]])
                st_s = -b_s / norm(b_s) * step
                pos_s += st_s
                pos_sph_s = Functions.cartesian2spherical(pos_s)
            end
            
            # Interpolate exactly to the surface r = psr.r
            r_last_s = norm(pos_s)
            pos_prev_s = pos_s - st_s
            r_prev_s = norm(pos_prev_s)
            t_s = (r_prev_s - psr.r) / (r_prev_s - r_last_s)
            surf_s = pos_prev_s + t_s * (pos_s - pos_prev_s)

            push!(x_s, surf_s[1])
            push!(y_s, surf_s[2])
            push!(z_s, surf_s[3])
        end

        # Assign the calculated anomalous shapes to the pulsar object
        psr.polar_caps = [[x_n, y_n, z_n], [x_s, y_s, z_s]]
        psr.pc = [x_n, y_n, z_n]  # Use the northern polar cap by default
    end

    """
    

    """
    function generate_open_obsolete!(psr, step=10)
        fv = psr.fields
        stepsnum = div(fv.rmax, step)

        # first polar cap only
        pc = psr.polar_caps[1]
        xs = pc[1]
        ys = pc[2]
        zs = pc[3]
        for (j,x) in enumerate(xs)
            pos = [xs[j], ys[j], zs[j]]
            pos_sph = Functions.cartesian2spherical(pos)
            b_sph = Field.bvac(pos_sph, psr.r, fv.beq)
            b = Functions.vec_spherical2cartesian(pos_sph, b_sph)
            push!(psr.open_lines, [[pos[1]], [pos[2]], [pos[3]]]) # adding initial position
            ml = psr.open_lines[j] # magnetic line add values to ml[1] is x ml[2] is y ml[3] is z
            for k in 1:stepsnum
                pos_sph = Functions.cartesian2spherical(pos)
                b_sph = Field.bvac(pos_sph, psr.r, fv.beq)
                b = Functions.vec_spherical2cartesian(pos_sph, b_sph)
                st = b / norm(b) * step
                pos += st # new position for magnetic line
                push!(ml[1], pos[1])
                push!(ml[2], pos[2])
                push!(ml[3], pos[3])
            end
        end
    end


    """
    generate_open!(psr; num=100)

    Traces the last open magnetic field lines from height rmax (defined in psr.nsfield)
    down to the stellar surface using the non-dipolar field (NSField).

    At rmax the last-open-line polar angle is estimated from the dipole formula.
    num uniformly spaced azimuthal points are placed on that circle and each line
    is integrated inward. Results are appended to psr.open_lines.

    # Arguments
    - num: number of points (field lines) distributed around the circle at rmax
    """
    function generate_open!(psr; num=100)
        nf = psr.nsfield
        step = nf.rmax / nf.size

        z_rmax = nf.rmax / psr.r  # rmax in stellar radii
        theta_rmax = Functions.theta_max(z_rmax, psr)
        phis = range(0, 2*pi, length=num+1)[1:num]

        for phi in phis
            pos = Functions.spherical2cartesian([nf.rmax, theta_rmax, phi])
            pos_sph = Functions.cartesian2spherical(pos)

            push!(psr.open_lines, [[pos[1]], [pos[2]], [pos[3]]])
            ml = psr.open_lines[end]

            while pos_sph[1] > psr.r
                b_sph = NSField.BSph(nf, pos_sph[1]/psr.r, pos_sph[2], pos_sph[3])
                b = Functions.vec_spherical2cartesian(pos_sph, [b_sph[1], b_sph[2], b_sph[3]])
                st = -b / norm(b) * step
                pos += st
                pos_sph = Functions.cartesian2spherical(pos)
                push!(ml[1], pos[1])
                push!(ml[2], pos[2])
                push!(ml[3], pos[3])
            end

            # Interpolate last point to the stellar surface
            n = length(ml[1])
            if n >= 2
                pos_prev = [ml[1][n-1], ml[2][n-1], ml[3][n-1]]
                pos_last = [ml[1][n],   ml[2][n],   ml[3][n]]
                r_prev = norm(pos_prev)
                r_last = norm(pos_last)
                t = (r_prev - psr.r) / (r_prev - r_last)
                pos_surface = pos_prev + t * (pos_last - pos_prev)
                ml[1][n] = pos_surface[1]
                ml[2][n] = pos_surface[2]
                ml[3][n] = pos_surface[3]
            end
        end

        # Fit ellipse to the last points of open lines (lying on the stellar surface)
        xs = [ml[1][end] for ml in psr.open_lines]
        ys = [ml[2][end] for ml in psr.open_lines]
        zs = [ml[3][end] for ml in psr.open_lines]
        psr.pc = [xs, ys, zs]
        psr.ellipse_fit = SolidBody.fit_ellipse([xs, ys, zs], psr.r)
    end


    function init_line_of_sight(psr; num=100)

         # calculates line of sight
        psr.line_of_sight = []
        psr.longitudes = zeros(num)

        α = deg2rad(psr.alpha)
        β = deg2rad(psr.beta)

        φ_s = Geometry.generate_uniform_phase_array(num, α ,β , psr.r_em, psr.p)
        #θ_array, ψ_array = Geometry.emission_points_from_phase(φ_s, α, β, psr.r_em, psr.p)
        θ_array, ψ_array = Geometry.emission_points_with_ar(φ_s, α, β, psr.r_em, psr.p)

        for i in eachindex(θ_array)
            ψ = ψ_array[i]
            θ = θ_array[i]
            push!(psr.line_of_sight, Functions.spherical2cartesian([psr.r_em, θ, ψ]))
            psr.longitudes[i] = rad2deg(φ_s[i])
        end

    end



    function calculate_line_of_sight_dipole(psr, step=10)

        if isnothing(psr.line_of_sight)
            println("Init line of sight first!")
            return
        end

        fv = psr.fields

        for point in psr.line_of_sight
            pos = copy(point)
            pos_sph = Functions.cartesian2spherical(pos)
            # new line with 
            push!(psr.los_lines, [Float64[], Float64[], Float64[]]) # push!(los[end][1], x) etc.
            push!(psr.los_lines[end][1], pos[1]) # x coordinate
            push!(psr.los_lines[end][2], pos[2]) # y coordinate
            push!(psr.los_lines[end][3], pos[3]) # z coordinate
            while (pos_sph[1] >= psr.r) 
                b_sph = Field.bvac(pos_sph, psr.r, fv.beq)
                b = Functions.vec_spherical2cartesian(pos_sph, b_sph)
                st = - b / norm(b) * step # negative step towards the surface
                pos += st # new position for magnetic line
                pos_sph = Functions.cartesian2spherical(pos)
                push!(psr.los_lines[end][1], pos[1])
                push!(psr.los_lines[end][2], pos[2])
                push!(psr.los_lines[end][3], pos[3])
            end
            
            # Correct last point - interpolate to surface
            line = psr.los_lines[end]
            n = length(line[1])
            
            # Second to last point (above surface)
            pos_prev = [line[1][n-1], line[2][n-1], line[3][n-1]]
            # Last point (below surface)
            pos_last = [line[1][n], line[2][n], line[3][n]]
            
            r_prev = norm(pos_prev)
            r_last = norm(pos_last)
            
            # Interpolation factor
            t = (r_prev - psr.r) / (r_prev - r_last)
            pos_surface = pos_prev + t * (pos_last - pos_prev)
            
            # Overwrite last point
            line[1][n] = pos_surface[1]
            line[2][n] = pos_surface[2]
            line[3][n] = pos_surface[3]

        end
        #println(size(psr.los_lines))
    end



    function calculate_line_of_sight(psr)

        if isnothing(psr.line_of_sight)
            println("Init line of sight first!")
            return
        end

        nf = psr.nsfield
        step = nf.rmax / nf.size

        for point in psr.line_of_sight
            pos = copy(point)
            pos_sph = Functions.cartesian2spherical(pos)

            # new line 
            push!(psr.los_lines, [Float64[], Float64[], Float64[]])
            push!(psr.los_lines[end][1], pos[1]) # x coordinate
            push!(psr.los_lines[end][2], pos[2]) # y coordinate
            push!(psr.los_lines[end][3], pos[3]) # z coordinate
            while (pos_sph[1] > psr.r) 
                #b_sph = Field.bvac(pos_sph, psr.r, fv.beq)
                b_sph = NSField.BSph(nf, pos_sph[1]/psr.r, pos_sph[2], pos_sph[3])
                b = Functions.vec_spherical2cartesian(pos_sph, b_sph)
                st = - b / norm(b) * step # negative step towards the surface
                pos += st # new position for magnetic line
                pos_sph = Functions.cartesian2spherical(pos)
                push!(psr.los_lines[end][1], pos[1])
                push!(psr.los_lines[end][2], pos[2])
                push!(psr.los_lines[end][3], pos[3])
            end
           
             # Correct last point - interpolate to surface
            line = psr.los_lines[end]
            n = length(line[1])
            
            # Second to last point (above surface)
            pos_prev = [line[1][n-1], line[2][n-1], line[3][n-1]]
            # Last point (below surface)
            pos_last = [line[1][n], line[2][n], line[3][n]]
            
            r_prev = norm(pos_prev)
            r_last = norm(pos_last)
            
            # Interpolation factor
            t = (r_prev - psr.r) / (r_prev - r_last)
            pos_surface = pos_prev + t * (pos_last - pos_prev)
            
            # Overwrite last point
            line[1][n] = pos_surface[1]
            line[2][n] = pos_surface[2]
            line[3][n] = pos_surface[3]


        end
    end




"""
    Traces the last open magnetic field lines using the full anomalous magnetic field 
    (with proper scaling from Field.b_anomalous).
    """
    function generate_open_anomalous!(psr; num=100)
        nf = psr.nsfield
        step = nf.rmax / nf.size

        z_rmax = nf.rmax / psr.r  # rmax in stellar radii
        theta_rmax = Functions.theta_max(z_rmax, psr)
        phis = range(0, 2*pi, length=num+1)[1:num]

        for phi in phis
            pos = Functions.spherical2cartesian([nf.rmax, theta_rmax, phi])
            pos_sph = Functions.cartesian2spherical(pos)

            push!(psr.open_lines, [[pos[1]], [pos[2]], [pos[3]]])
            ml = psr.open_lines[end]

            while pos_sph[1] > psr.r
                # Using the dedicated anomaly function with real units
                b_sph = NSField.b_anomalous(pos_sph, psr)
                b = Functions.vec_spherical2cartesian(pos_sph, [b_sph[1], b_sph[2], b_sph[3]])
                st = -b / norm(b) * step
                pos += st
                pos_sph = Functions.cartesian2spherical(pos)
                
                push!(ml[1], pos[1])
                push!(ml[2], pos[2])
                push!(ml[3], pos[3])
            end

            # Interpolate last point to the stellar surface
            n = length(ml[1])
            if n >= 2
                pos_prev = [ml[1][n-1], ml[2][n-1], ml[3][n-1]]
                pos_last = [ml[1][n],   ml[2][n],   ml[3][n]]
                r_prev = norm(pos_prev)
                r_last = norm(pos_last)
                t = (r_prev - psr.r) / (r_prev - r_last)
                pos_surface = pos_prev + t * (pos_last - pos_prev)
                
                ml[1][n] = pos_surface[1]
                ml[2][n] = pos_surface[2]
                ml[3][n] = pos_surface[3]
            end
        end
    end

    """
    Calculates line of sight using the full anomalous magnetic field.
    """
    function calculate_line_of_sight_anomalous!(psr)
        if isnothing(psr.line_of_sight)
            println("Init line of sight first!")
            return
        end

        nf = psr.nsfield
        step = nf.rmax / nf.size

        for point in psr.line_of_sight
            pos = copy(point)
            pos_sph = Functions.cartesian2spherical(pos)

            push!(psr.los_lines, [Float64[], Float64[], Float64[]])
            push!(psr.los_lines[end][1], pos[1]) 
            push!(psr.los_lines[end][2], pos[2]) 
            push!(psr.los_lines[end][3], pos[3]) 
            
            while (pos_sph[1] > psr.r) 
                # Using the dedicated anomaly function with real units
                b_sph = NSField.b_anomalous(pos_sph, psr)
                b = Functions.vec_spherical2cartesian(pos_sph, [b_sph[1], b_sph[2], b_sph[3]])
                st = - b / norm(b) * step 
                pos += st 
                pos_sph = Functions.cartesian2spherical(pos)
                
                push!(psr.los_lines[end][1], pos[1])
                push!(psr.los_lines[end][2], pos[2])
                push!(psr.los_lines[end][3], pos[3])
            end
           
            # Correct last point - interpolate to surface
            line = psr.los_lines[end]
            n = length(line[1])
            
            pos_prev = [line[1][n-1], line[2][n-1], line[3][n-1]]
            pos_last = [line[1][n], line[2][n], line[3][n]]
            
            r_prev = norm(pos_prev)
            r_last = norm(pos_last)
            
            t = (r_prev - psr.r) / (r_prev - r_last)
            pos_surface = pos_prev + t * (pos_last - pos_prev)
            
            line[1][n] = pos_surface[1]
            line[2][n] = pos_surface[2]
            line[3][n] = pos_surface[3]
        end
    end
"""
    generate_closed_anomalous!(psr; shells=3, num_lines=10, r_max_factor=3.0)

    Generates closed magnetic field lines considering magnetic anomalies.
    Saves the results directly into psr.closed_lines.
    """
    function generate_closed_anomalous!(psr; shells=3, num_lines=10, r_max_factor=3.0)
        nf = psr.nsfield
        psr.closed_lines = [] # Clear previous lines if any
        step_val = nf.rmax / nf.size

        # Starting points at the magnetic equator (theta = pi/2)
        reqs = range(psr.r * 1.2, psr.r * r_max_factor, length=shells)
        phis = range(0, 2pi, length=num_lines+1)[1:end-1]

        for req in reqs
            for phi in phis
                pos_sph = [req, pi/2, phi]
                pos = Functions.spherical2cartesian(pos_sph)

                # --- 1. Integrate "along" the B vector (towards one pole) ---
                pos1 = copy(pos)
                pos_sph1 = copy(pos_sph)
                line1 = [[pos1[1]], [pos1[2]], [pos1[3]]]

                # Limit to 5000 points to prevent infinite loops
                while pos_sph1[1] > psr.r && length(line1[1]) < 5000
                    b_sph = NSField.BSph(nf, pos_sph1[1]/psr.r, pos_sph1[2], pos_sph1[3])
                    b = Functions.vec_spherical2cartesian(pos_sph1, [b_sph[1], b_sph[2], b_sph[3]])
                    st = b / norm(b) * step_val
                    pos1 += st
                    pos_sph1 = Functions.cartesian2spherical(pos1)
                    
                    push!(line1[1], pos1[1])
                    push!(line1[2], pos1[2])
                    push!(line1[3], pos1[3])
                end

                # --- 2. Integrate "against" the B vector (towards the other pole) ---
                pos2 = copy(pos)
                pos_sph2 = copy(pos_sph)
                line2 = [Float64[], Float64[], Float64[]]

                while pos_sph2[1] > psr.r && length(line2[1]) < 5000
                    b_sph = NSField.BSph(nf, pos_sph2[1]/psr.r, pos_sph2[2], pos_sph2[3])
                    b = Functions.vec_spherical2cartesian(pos_sph2, [b_sph[1], b_sph[2], b_sph[3]])
                    st = -b / norm(b) * step_val # NEGATIVE step to go the other way
                    pos2 += st
                    pos_sph2 = Functions.cartesian2spherical(pos2)
                    
                    push!(line2[1], pos2[1])
                    push!(line2[2], pos2[2])
                    push!(line2[3], pos2[3])
                end

                # --- 3. Merge both halves into a single continuous line ---
                xs = vcat(reverse(line2[1]), line1[1])
                ys = vcat(reverse(line2[2]), line1[2])
                zs = vcat(reverse(line2[3]), line1[3])
                
                push!(psr.closed_lines, [xs, ys, zs])
            end
        end
    end
end # module end