module LinesNew
    include("functions.jl")
    include("geometry.jl")
    import ..NSField


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
            push!(psr.line_of_sight, [psr.r_em, θ, ψ]) # spherical coordinates
            psr.longitudes[i] = rad2deg(φ_s[i])
        end

    end


    function calculate_line_of_sight(psr)

        if isnothing(psr.line_of_sight)
            println("Init line of sight first!")
            return
        end

        nf = psr.nsfield
        step = nf.rmax / nf.size

        for point in psr.line_of_sight
            pos_sph = point
            pos = Functions.spherical2cartesian(pos_sph)

            #println(NSField.BSph(nf, pos_sph[1], pos_sph[2], pos_sph[3]))
            println(pos_sph[1]," ", pos_sph[2]," ", pos_sph[3])
            # TODO work here...
            return

            # new line 
            push!(psr.los_lines, [Float64[], Float64[], Float64[]]) # push!(los[end][1], r) etc.
            push!(psr.los_lines[end][1], pos_sph[1]) # r coordinate
            push!(psr.los_lines[end][2], pos_sph[2]) # theta coordinate
            push!(psr.los_lines[end][3], pos_sph[3]) # phi coordinate
            while (pos_sph[1] > psr.r) 
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