module Field
    
    include("functions.jl")
    using LinearAlgebra


    mutable struct Test
        rmax 
        size
        beq # Magnetic field strength at the stellar equator
        locations
        magnetic_lines
        magnetic_fields
        
        function Test()
            rmax = 50000 # in meters
            size = 15
            new(rmax, size, nothing, [], [], [])
        end 
    end
    
    



    """
    dipole(r, theta)

    Dipole magnetic field

    # Arguments

    - r: in stellar radius units
    - theta: in radians

    """
    function dipole(r, theta)
        b_r = 2. * cos(theta) / (r ^ 3)
        b_theta = sin(theta) / (r ^ 3)
        return (b_r, b_theta, 0.)
    end
    function beq(p, pdot)
        return 3.2e19*sqrt(p * pdot) # 3.2e19 usunąłem mnożenie przez to
    end


    """
    Spherical components of magnetic field for an aligned rotator (Cerutti, 2016 p.3) for Q=0
    """
    function bvac(pos_sph, rstar, beq)
        #println(beq)
        r = pos_sph[1]
        theta = pos_sph[2]
        phi = pos_sph[3]
        #println(rstar, r, beq)
        br = beq * (rstar / r) ^ 3 * 2 * cos(theta)
        btheta = beq * (rstar / r) ^ 3 * sin(theta)
        bphi = 0
        return(br, btheta, bphi)
    end

   
    """
    Calculate the dipole magnetic field on a grid of spherical coordinates
    and store the field vectors (converted to Cartesian) in psr.fields.
    """
function calculate_dipole!(psr)
    fv = psr.fields
    fv.beq = beq(psr.p, psr.pdot)   # magnetic field at equator from P and Pdot

    # define a 3D grid in spherical coordinates
    rs = LinRange(psr.r, fv.rmax, fv.size)
    thetas = LinRange(0, pi, fv.size)
    phis = LinRange(0, 2pi, fv.size)

    # loop over all grid points
    for i in 1:fv.size
        for j in 1:fv.size
            for k in 1:fv.size
                pos_sph = [rs[i], thetas[j], phis[k]]

                # magnetic field in spherical coordinates
                b_sph = bvac(pos_sph, psr.r, fv.beq)

                # save position and field vector in Cartesian coordinates
                push!(fv.locations, Functions.spherical2cartesian(pos_sph))
                push!(fv.magnetic_fields, Functions.vec_spherical2cartesian(pos_sph, b_sph))
            end
        end
    end
end


"""
Generate magnetic field lines starting from the stellar surface
for many (theta, phi) positions. Each line is followed step by step
by integrating along the field direction.
"""
function generate_lines!(psr; step=10, stepsnum=20000, phi=nothing)
    fv = psr.fields
    r = psr.r  # start radius at stellar surface

    # sampling in theta and phi
    thetas = LinRange(0, pi, fv.size)
    if phi === nothing
        phis = LinRange(0, 2pi, fv.size+1)[1:end-1]  # avoid duplicate point at 2pi
    else
        phis = [phi, phi+pi]   # symmetric pair if single phi requested
    end

    # loop over starting positions
    for i in 1:length(thetas)
        for j in 1:length(phis)
            pos_sph = [r, thetas[i], phis[j]]
            b_sph = Field.bvac(pos_sph, psr.r, fv.beq)

            pos = Functions.spherical2cartesian(pos_sph)
            b = Functions.vec_spherical2cartesian(pos_sph, b_sph)

            # initialize a new magnetic line with starting point
            push!(fv.magnetic_lines, [[pos[1]], [pos[2]], [pos[3]]])
            ml = fv.magnetic_lines[end]

            posb = copy(pos)
            step = abs(step)  # step size along the field line

            # follow the line for a fixed number of steps
            for k in 1:stepsnum
                posb_sph = Functions.cartesian2spherical(posb)

                # if we cross inside the star
                if posb_sph[1] < psr.r
                    if size(ml[1], 1) > 2
                        break   # stop if we already have a valid line
                    else
                        step = -step  # otherwise reverse direction
                    end
                end

                # recalculate field at new position
                b_sph = Field.bvac(posb_sph, psr.r, fv.beq)
                b = Functions.vec_spherical2cartesian(posb_sph, b_sph)

                # move along normalized field vector
                st = b / norm(b) * step
                posb += st

                # append new point to the magnetic line
                push!(ml[1], posb[1])
                push!(ml[2], posb[2])
                push!(ml[3], posb[3])
            end
        end
    end
end
function pc(psr; phi_num=100)
        theta =Functions.theta_max(1, psr)
        phis = range(0, 2*pi, length=phi_num)
        x=Array{Float64}(undef, phi_num)
        y=Array{Float64}(undef, phi_num)
        z=Array{Float64}(undef, phi_num)

        for(i, phi) in enumerate(phis)
            car = Functions.spherical2cartesian([psr.r, theta, phi])
            x[i] = car[1]
            y[i] = car[2]
            z[i] = car[3]
        end
        psr.pc = [x, y, z]
        
end

"""
Calculate the polar cap contour on the stellar surface.
Only one (northern) polar cap is computed here.
"""
function calculate_polarcap!(psr; phi_num=100)
    theta = Functions.theta_max(1, psr)   # maximum opening angle for open field lines
    phis = range(0, 2*pi, length=phi_num)

    # store spherical coordinates of the polar cap rim
    psr.pc = Functions.spherical2cartesian.([[psr.r, theta, phi] for phi in phis])
end



"""
Generate magnetic field lines starting only from the edge of the polar caps
(north and south). Each starting point lies on the polar cap rim, and
field lines are traced step by step along the field direction.
"""
function generate_polarcap_lines!(psr; phi_num=10, step=10, stepsnum=200000)
    fv = psr.fields
    r = psr.r
    
    # polar cap angle
    theta = Functions.theta_max(1, psr)
    phis = range(0, 2pi, length=phi_num)

    for ph in phis
        # start lines from both north and south polar caps
        for theta0 in (theta, pi - theta)
            pos_sph = [r, theta0, ph]
            b_sph = Field.bvac(pos_sph, psr.r, fv.beq)
            pos = Functions.spherical2cartesian(pos_sph)

            # initialize new line
            push!(fv.magnetic_lines, [[pos[1]], [pos[2]], [pos[3]]])
            ml = fv.magnetic_lines[end]

            posb = copy(pos)
            step_dir = abs(step)

            # follow the field line
            for k in 1:stepsnum
                posb_sph = Functions.cartesian2spherical(posb)

                # if we enter inside the star
                if posb_sph[1] < psr.r
                    if size(ml[1], 1) > 2
                        break   # stop if the line is already long enough
                    else
                        step_dir = - step_dir  # reverse integration direction
                    end
                end

                # recalculate field
                b_sph = Field.bvac(posb_sph, psr.r, fv.beq)
                b = Functions.vec_spherical2cartesian(posb_sph, b_sph)

                # step along field
                st = b / norm(b) * step_dir
                posb += st

                # save point to the line
                push!(ml[1], posb[1])
                push!(ml[2], posb[2])
                push!(ml[3], posb[3])
            end
        end
    end
end

    """
    dipolar component of magnetic field at the surface based on x, y components
    |B_d| = 2 at the pole
    """
    function bd(x, y, psr)
        d = sqrt(x ^ 2 + y ^ 2)
        theta = asin(d / psr.r)
        bd_sph = dipole(1, theta)
        #println(dipole(1, 0))
        bd = Functions.spherical2cartesian(bd_sph)
        return bd
    end

function bc(psr, height_m; phi_num=100)
        # 1. Calculate the opening angle (theta) for this height
        # Function expects ratio of radius to stellar radius
        theta = Functions.theta_max(height_m / psr.r, psr)
        
        # 2. Prepare the Rotation
        # We need to rotate from the Standard Z-axis to the Pulsar's Magnetic Axis
        z_axis = [0.0, 0.0, 1.0]
        
        # Get magnetic axis direction (normalized)
        mag_axis = Functions.spherical2cartesian(psr.magnetic_axis)
        mag_axis = mag_axis / norm(mag_axis)
        
        # Calculate rotation axis (k) and angle (alpha) to get from Z to Magnetic Axis
        # Axis k is perpendicular to both Z and Magnetic Axis
        k = cross(z_axis, mag_axis)
        sin_alpha = norm(k)
        cos_alpha = dot(z_axis, mag_axis)
        
        # Safety check: if magnetic axis is already Z, k is zero length
        rotation_needed = sin_alpha > 1e-6
        if rotation_needed
            k = k / sin_alpha # Normalize the rotation axis
        end

        # 3. Generate and Rotate Points
        phis = range(0, 2*pi, length=phi_num)
        xs = Float64[]
        ys = Float64[]
        zs = Float64[]
        
        for phi in phis
            # A. Generate point in standard Z-aligned spherical coordinates
            # This is the "fashion" you asked for: using standard spherical definitions
            p_standard = Functions.spherical2cartesian([height_m, theta, phi])
            
            # B. Rotate the point if necessary
            if rotation_needed
                # Rodrigues' Rotation Formula
                # v_rot = v*cos(a) + (k x v)*sin(a) + k*(k.v)*(1-cos(a))
                
                term1 = p_standard * cos_alpha
                term2 = cross(k, p_standard) * sin_alpha
                term3 = k * dot(k, p_standard) * (1 - cos_alpha)
                
                p_final = term1 + term2 + term3
            else
                p_final = p_standard
            end
            
            push!(xs, p_final[1])
            push!(ys, p_final[2])
            push!(zs, p_final[3])
        end
        
        return [xs, ys, zs]
    end
function generate_los!(ax, psr; height=500000.0, points=500)
    # 1. Przygotowanie geometrii rotacji (tak jak w plot_los!)
    # Upewniamy się, że kąty są w radianach
    zeta = deg2rad(psr.alpha) + deg2rad(psr.beta)
    
    k = Functions.spherical2cartesian(psr.rotation_axis)
    k = k / norm(k)

    # Baza prostopadła do osi rotacji
    aux = [0.0, 0.0, 1.0]
    if abs(dot(k, aux)) > 0.99; aux = [1.0, 0.0, 0.0]; end
    u = cross(k, aux); u = u/norm(u)
    v = cross(k, u)

    # 2. Definicja Osi Magnetycznej
    # W Twoim kodzie (generate_polarcap_lines) linie pola są generowane sferycznie od theta,
    # co implikuje, że oś magnetyczna to oś Z [0, 0, 1].
    m_axis = [0.0, 0.0, 1.0] 

    # 3. Limit Pola Magnetycznego (Szerokość wiązki)
    # Kąt czapy polarnej na powierzchni
    theta_pc = Functions.theta_max(1, psr)
    
    # Fizyka: Linie pola otwierają się proporcjonalnie do pierwiastka z promienia.
    # Na wysokości 'height' stożek pola jest szerszy.
    # width(r) ~ theta_pc * sqrt(r / R_star)
    beam_width_at_height = theta_pc * sqrt((psr.r + height) / psr.r)

    # 4. Generowanie Linii
    # Używamy NaN, aby oddzielić od siebie poszczególne linie w jednej tablicy (szybkie rysowanie)
    xs, ys, zs = Float64[], Float64[], Float64[]

    cz, sz = cos(zeta), sin(zeta)
    phis = range(0, 2pi, length=points)

    for phi in phis
        # Wektor jednostkowy Linii Wzroku (LOS) w danej fazie
        dir = k * cz + (u * cos(phi) + v * sin(phi)) * sz
        
        # Oblicz kąt między LOS a Osią Magnetyczną
        # clamp służy do uniknięcia błędów numerycznych (gdyby dot wyszedł 1.0000000002)
        angle_to_mag = acos(clamp(dot(dir, m_axis), -1, 1))
        
        # Sprawdzamy czy LOS wpada w "lejek" otwartych linii pola (Północ lub Południe)
        # Używamy szerszego kąta (na wysokości), aby pokazać całą aktywną objętość
        is_active_north = angle_to_mag < beam_width_at_height
        is_active_south = (pi - angle_to_mag) < beam_width_at_height

        if is_active_north || is_active_south
            # Punkt startowy (na powierzchni gwiazdy)
            p_start = dir * psr.r
            # Punkt końcowy (na zadanej wysokości)
            p_end = dir * (psr.r + height)
            
            # Dodajemy punkty do tablicy, oddzielając je NaN (tworzy przerwę)
            push!(xs, p_start[1], p_end[1], NaN)
            push!(ys, p_start[2], p_end[2], NaN)
            push!(zs, p_start[3], p_end[3], NaN)
        end
    end
    
    # 5. Rysowanie
    # color=:red - aby odróżnić od zielonego okręgu LOS i niebieskich linii pola
    if !isempty(xs)
        lines!(ax, xs, ys, zs, color=:red, linewidth=2, label="Active Beam")
    else
        println("Warning: LOS does not intersect the open field lines (Pulse not visible).")
    end
end
end # module end