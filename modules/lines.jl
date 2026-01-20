module Lines
    using LinearAlgebra

    include("field.jl")
    include("functions.jl")
    include("transformations.jl")
    include("signal.jl")

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

    function generate_open!(psr, step=10)
        fv = psr.fields
        stepsnum = div(fv.rmax, step)

        # two polar caps
        for (i,pc) in enumerate(psr.polar_caps)
            # points at the polar cap
            xs = pc[1]
            ys = pc[2]
            zs = pc[3]
            # goint other direction at south pole
            if i == 2
                step = -step
            end
            # all points
            for (j,x) in enumerate(xs) 
                pos = [xs[j], ys[j], zs[j]]
                pos_sph = Functions.cartesian2spherical(pos)
                b_sph = Field.bvac(pos_sph, psr.r, fv.beq)
                b = Functions.vec_spherical2cartesian(pos_sph, b_sph)
                push!(psr.open_lines[i], [[pos[1]], [pos[2]], [pos[3]]]) # adding initial position
                ml = psr.open_lines[i][j] # magnetic line add values to ml[1] is x ml[2] is y ml[3] is z 
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
    end

    function init_line_of_sight(psr; num=10)

        theta_max = Functions.theta_max(psr.r_em/psr.r, psr) # at the emission height

        # determine length to have ~num points
        init_length = 1000
        phis = range(0, 2π, length=init_length)
        points_num = 0
        for phi in phis
            vec = Transformations.beaming(Functions.spherical2cartesian(psr.rotation_axis), deg2rad(psr.alpha+psr.beta), phi)
            if vec[2] <= theta_max
                points_num += 1
            end
        end
        if points_num == 0
            println("No points in open field line region! Change beta?")
            return
        end

        len = floor(Int, num / points_num * init_length)
        # calculates line of sight
        psr.line_of_sight = []
        longs = []
        phis = range(0, 2π, length=len)

        for phi in phis
            vec = Transformations.beaming(Functions.spherical2cartesian(psr.rotation_axis), deg2rad(psr.alpha+psr.beta), phi)
            if vec[2] <= theta_max
                push!(psr.line_of_sight, Functions.spherical2cartesian(vec)/psr.r* psr.r_em)
                push!(longs, rad2deg(vec[3])) # TODO I am not sure!
            end
        end

        # Odwiń żeby było ciągłe
        for i in 2:length(longs)
            while longs[i] - longs[i-1] > π
                longs[i] -= 2π
            end
            while longs[i] - longs[i-1] < -π
                longs[i] += 2π
            end
        end
    
        # center at zero
        lon_center = (minimum(longs) + maximum(longs)) / 2
        longs .-= lon_center
        psr.longitudes = longs

    end

    function calculate_line_of_sight(psr, step=10)

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


end # module end