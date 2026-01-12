module Lines
    using LinearAlgebra # Niezbędne dla: norm, dot, cross
    include("functions.jl")

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
end
#przerzucić wszystko co związane z liniami pola magnetycznego i los tutaj