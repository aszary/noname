module Field



    mutable struct TestClass
        locations
        magnetic_lines 
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
    Calculates magnetic fields using TestClass.
    """
    function calculate_dipole!(psr)

        return 
        fv = psr.field_vacuum
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
                    e_sph = evac(pos_sph, psr.r, fv.beq, psr.omega)
                    push!(fv.locations, Functions.spherical2cartesian(pos_sph))
                    push!(fv.magnetic, Functions.vec_spherical2cartesian(pos_sph, b_sph))
                    push!(fv.electric, Functions.vec_spherical2cartesian(pos_sph, e_sph))
                end
            end
        end
    end



end # module end