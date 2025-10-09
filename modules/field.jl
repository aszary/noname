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

"""
Calculate the polar cap contour on the stellar surface.
Only one (northern) polar cap is computed here.
"""
function calculate_polarcap!(psr; phi_num=100)
    theta = Functions.theta_max(1, psr)   # maximum opening angle for open field lines
    phis = range(0, 2*pi, length=phi_num)

    # store spherical coordinates of the polar cap rim
    psr.pc = [[psr.r, theta, ph] for ph in phis]
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


end # module end