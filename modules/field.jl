module Field
    using LinearAlgebra

    include("functions.jl")


    mutable struct Test
        rmax # maximum radius for which we will calculate magnetic field lines
        size # number of points to calculate
        beq # Magnetic field strength at the stellar equator
        locations # locations of points to calculate magnetic field
        magnetic_lines # magnetic field lines
        magnetic_fields # magnetic field vectors 
        function Test()
            rmax = 50000 # 50 km 
            size = 10
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
                    push!(fv.locations, Functions.spherical2cartesian(pos_sph))
                    push!(fv.magnetic_fields, Functions.vec_spherical2cartesian(pos_sph, b_sph))
                end
            end
        end
    end


    function generate_lines!(psr; step=10, stepsnum=20000, phi=nothing)
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
                pos = Functions.spherical2cartesian(pos_sph)                
                b = Functions.vec_spherical2cartesian(pos_sph, b_sph)
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