module Field
    using LinearAlgebra
    include("functions.jl")


    mutable struct Test
        rmax # Maximum radius in stellar radius units
        size 
        beq # Magnetic field strength at the stellar equator
        locations
        magnetic_lines 
        magnetic_fields
        magnetic_lines_from_cap   
        

         
        function Test()
           
            rmax = 50e3
            size= 12

            new( rmax, size, nothing, [], [], [], [])
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



    """
    Magnetic field strength at the stellar equator (Handbook p. 267)
    """
    function beq(p, pdot)
        return 3.2e19 * sqrt(p * pdot)
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
    Calculates magnetic fields using Test.
    """
    function calculate_dipole!(psr)
 
        fv = psr.fields
        fv.beq = beq(psr.p, psr.pdot)
        #println(fv)

        rs = LinRange(psr.r, fv.rmax, fv.size)
        thetas = LinRange(0, pi, fv.size)
        phis = LinRange(0, 2pi, fv.size)

        for i in 1:fv.size
            for j in 1:fv.size
                for k in 1:fv.size
                    #println(rs[i], " ", thetas[j], " ", phis[k])
                    pos_sph = [rs[i], thetas[j], phis[k]]
                    b_sph = bvac(pos_sph, psr.r, fv.beq)
                    #e_sph = evac(pos_sph, psr.r, fv.beq, psr.omega)
                    push!(fv.locations, Functions.spherical2cartesian(pos_sph))
                    push!(fv.magnetic_fields, Functions.vec_spherical2cartesian(pos_sph, b_sph))
                    #push!(fv.electric, Functions.vec_spherical2cartesian(pos_sph, e_sph))
                end
            end
        end
    end


    """
    Generates magnetic and electric field lines for vacuum around neutron star
    step in meters
    """
    function generate_vacuum!(psr; step=10, stepsnum=20000, phi=nothing)
        fv = psr.fields

        # starting points
        r = psr.r
        #thetas = LinRange(0, pi/2, fv.size)
        #thetas = LinRange(0, pi, fv.size+1)[1:end-1]
        thetas = LinRange(0, pi, fv.size)
        if phi === nothing
            phis = LinRange(0, 2pi, fv.size+1)[1:end-1] # get rid of last point
        else
            phis = [phi, phi+pi]
        end

        # TODO add the second half! done? but too many lines?
        
        for i in 1:size(thetas)[1]
            for j in 1:size(phis)[1]
                pos_sph = [r, thetas[i], phis[j]]
                b_sph = Field.bvac(pos_sph, psr.r, fv.beq)
                #e_sph = Field.evac(pos_sph, psr.r, fv.beq, psr.omega)
                pos = Functions.spherical2cartesian(pos_sph)                
                b = Functions.vec_spherical2cartesian(pos_sph, b_sph)
                #e = Functions.vec_spherical2cartesian(pos_sph, e_sph)
                push!(fv.magnetic_lines, [[pos[1]], [pos[2]], [pos[3]]]) # adding initial position
                ml = fv.magnetic_lines[end]
                posb = copy(pos)
                step = abs(step) # start with positive step
                for k in 1:stepsnum
                    # new fields
                    posb_sph = Functions.cartesian2spherical(posb)
                    if posb_sph[1] < psr.r
                        # going the other direction if needed (e.g. southern hemisphere) or break
                        if size(ml[1], 1) > 2
                            #println("$i $j $k - break")
                            break
                        else
                            step = - step
                            #println("$i $j $k - minus")
                        end
                    end
                    #println(k)
                    b_sph = Field.bvac(posb_sph, psr.r, fv.beq)
                    b = Functions.vec_spherical2cartesian(posb_sph, b_sph)
                    st = b / norm(b) * step
                    posb += st # new position for magnetic
                    push!(ml[1], posb[1])
                    push!(ml[2], posb[2])
                    push!(ml[3], posb[3])
                end
                #push!(fv.electric_lines, [[pos[1]], [pos[2]], [pos[3]]]) # adding initial position
                #el = fv.electric_lines[end]
                #pose = copy(pos)
                #step = abs(step) # start with positive step
                #for k in 1:stepsnum
                #    pose_sph = Functions.cartesian2spherical(pose)
                #    if pose_sph[1] < psr.r
                #        # going the other direction if needed (e.g. southern hemisphere) or break
                #        if size(el[1], 1) > 2
                #            break
                #        else
                #            step = - step
                #        end
                #    end
                #    e_sph = Field.evac(pose_sph, psr.r, fv.beq, psr.omega)
                #    e = Functions.vec_spherical2cartesian(pose_sph, e_sph)
                #    st = e / norm(e) * step
                #    pose += st # new position for magnetic
                #    push!(el[1], pose[1])
                #    push!(el[2], pose[2])
                #    push!(el[3], pose[3])
                #end
            end
        end
    end

    function pc(psr; phi_num=10)
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

    function generate_polar_cap!(psr)
        θ_pc = asin(sqrt(psr.r / psr.r_lc))   # polar cap opening angle
        φs = range(0, 2π; length=200)         # evenly spaced azimuths

            # North cap (around +z)
        north = [(psr.r, θ_pc, φ) for φ in φs]

            # South cap (around -z, just mirror across equator)
        south = [(psr.r, π - θ_pc, φ) for φ in φs]

            # Store both
        psr.polar_cap = vcat(north, south)

    end


    function generate_magnetic_lines_from_polar_cap!(psr; num_lines_from_cap=50, step=100, stepsnum=20000)

        f = psr.fields
        empty!(f.magnetic_lines_from_cap)    # clear previous cap lines

        for cap_sign in (-1, 1)  # -1 = south, 1 = north
            for i in 1:num_lines_from_cap
                # --- starting point on the star surface ---
                θ_pc = asin(sqrt(psr.r / psr.r_lc))
                φ = 2π * rand()                      # random azimuth
                r0 = psr.r
                θ0 = cap_sign == 1 ? θ_pc : π - θ_pc
                φ0 = φ

                pos = Functions.spherical2cartesian([r0, θ0, φ0])
                line_x, line_y, line_z = [pos[1]], [pos[2]], [pos[3]]

                # integrate along the field
                step_sign = cap_sign                  # north +1, south -1
                for k in 1:stepsnum
                    pos_sph = Functions.cartesian2spherical(pos)
                    if pos_sph[1] > f.rmax; break; end

                    b_sph = Field.bvac(pos_sph, psr.r, f.beq)
                    b = Functions.vec_spherical2cartesian(pos_sph, b_sph)
                    pos += step_sign * step * b / norm(b)

                    push!(line_x, pos[1])
                    push!(line_y, pos[2])
                    push!(line_z, pos[3])
                end

                push!(f.magnetic_lines_from_cap, (line_x, line_y, line_z))
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
    


    


end # module end