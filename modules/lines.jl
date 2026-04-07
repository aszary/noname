module Lines
    using LinearAlgebra
    import Main.NoName.NSField
    import Main.NoName.Field
    import Main.NoName.Functions
    import Main.NoName.Transformations
    import Main.NoName.Geometry
    import Main.NoName.SolidBody


    """
    Simple function assuming dipolar configuration of magnetic field at the surface. To be OBSOLETE.  

    """
    function calculate_polarcaps!(psr; phi_num=100)
        theta = Functions.theta_max(1, psr)
        phis = range(0, 2*pi, length=phi_num)
        x = Array{Float64}(undef, phi_num)
        y = Array{Float64}(undef, phi_num)
        z = Array{Float64}(undef, phi_num)
        x2 = Array{Float64}(undef, phi_num)
        y2 = Array{Float64}(undef, phi_num)
        z2 = Array{Float64}(undef, phi_num)
    
        for (i,ph) in enumerate(phis)
            ca = Functions.spherical2cartesian([psr.r, theta, ph])
            x[i] = ca[1]
            y[i] = ca[2]
            z[i] = ca[3]
            ca2 = Functions.spherical2cartesian([psr.r, pi - theta, ph]) # south pole
            x2[i] = ca2[1]
            y2[i] = ca2[2]
            z2[i] = ca2[3]
        end
        psr.polar_caps = [[x, y, z], [x2, y2, z2]]
        psr.pc = [x, y, z]
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






end # module end