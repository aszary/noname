module Lines
    using LinearAlgebra # Niezbędne dla: norm, dot, cross
    include("functions.jl")
    include("field.jl")
    include("transformations.jl")
    include("geometry.jl")

    """
    Core geometry calculator.
    Calculates the path of the Line of Sight (LOS).
    Returns a list of tuples: (phase, p_surface, p_height)
    only for the moments when the pulsar is visible (active).
    """
    function calculate_los(psr; height=500000.0, points=1000)
        # 1. Wyczyść tablicę w obiekcie pulsara
        psr.line_of_sight = [] 

        # 2. Geometria
        zeta = deg2rad(psr.alpha) + deg2rad(psr.beta)
        
        k = Functions.spherical2cartesian(psr.rotation_axis)
        k = k / norm(k)
        
        aux = [0.0, 0.0, 1.0]
        if abs(dot(k, aux)) > 0.99; aux = [1.0, 0.0, 0.0]; end
        u = cross(k, aux); u = u/norm(u)
        v = cross(k, u)

        m_vec = Functions.spherical2cartesian(psr.magnetic_axis)
        m_vec = m_vec / norm(m_vec)

        # 3. Limity (szerokość wiązki na wysokości)
        theta_pc = Functions.theta_max(1, psr)
        scaling_factor = sqrt(psr.r / (psr.r + height))
        beam_width_at_height = asin(sin(theta_pc) / scaling_factor)

        # 4. Pętla obliczeniowa
        phis = range(0, 2pi, length=points)
        cz, sz = cos(zeta), sin(zeta)
        
        for phi in phis
            # Wektor LOS na wysokości
            dir = k * cz + (u * cos(phi) + v * sin(phi)) * sz
            
            # Sprawdzenie widoczności (czy jesteśmy w wiązce)
            cos_theta = clamp(dot(dir, m_vec), -1.0, 1.0)
            theta_los = acos(cos_theta)
            
            in_north = theta_los < beam_width_at_height
            in_south = (pi - theta_los) < beam_width_at_height

            if in_north || in_south
                # Rzutowanie z powrotem na powierzchnię (dipol)
                sin_theta_surf = sin(theta_los) * scaling_factor
                
                perp = dir - m_vec * cos_theta
                if norm(perp) > 1e-10; perp = perp / norm(perp); else; perp = [0.0, 0.0, 0.0]; end

                cos_theta_surf = sqrt(1 - sin_theta_surf^2)
                if in_south; cos_theta_surf = -cos_theta_surf; end

                dir_surf = m_vec * cos_theta_surf + perp * sin_theta_surf
                
                # Współrzędne
                p_surf = dir_surf * psr.r        # Punkt na powierzchni (Iskry)
                p_end = dir * (psr.r + height)   # Punkt na wysokości (Obserwator)
                
                # --- KLUCZOWA ZMIANA ---
                # Kolejność: [Punkt_Wysokości, Punkt_Powierzchni]
                # Dzięki temu line[1][end] to p_surf[1]
                xs = [p_end[1], p_surf[1]]
                ys = [p_end[2], p_surf[2]]
                zs = [p_end[3], p_surf[3]]
                
                push!(psr.line_of_sight, [xs, ys, zs])
            end
        end
    end
    function init_line_of_sight(psr; num=10)

         # calculates line of sight
        psr.line_of_sight = []
        psr.longitudes = zeros(num)

        α = deg2rad(psr.alpha)
        β = deg2rad(psr.beta)

        φ_s = Geometry.generate_uniform_phase_array(num, α ,β , psr.r_em, psr.p)
        θ_array, ψ_array = Geometry.emission_points_from_phase(φ_s, α, β, psr.r_em, psr.p)

        for i in eachindex(θ_array)
            ψ = ψ_array[i]
            θ = θ_array[i]
            push!(psr.line_of_sight, Functions.spherical2cartesian([psr.r_em, θ, ψ]))
            psr.longitudes[i] = rad2deg(φ_s[i])
        end

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

end
#przerzucić wszystko co związane z liniami pola magnetycznego i los tutaj