module Functions 
    using PhysicalConstants.CODATA2018 # in SI units


    """
    Calculates rlc in meters
    """
    function rlc(p, pdot)
        c = SpeedOfLightInVacuum.val # no units hereafter
        return c * p / (2 * pi)
    end

    """
    spherical2cartesian(spherical)

    Converts spherical cordinates to cartesian ones...
    """
    function spherical2cartesian(spherical)
        x = spherical[1] * sin(spherical[2]) * cos(spherical[3])
        y = spherical[1] * sin(spherical[2]) * sin(spherical[3])
        z = spherical[1] * cos(spherical[2])
        return [x, y, z]
    end


    """
    Converts vector in spehrical coordinates (vec_sph) at position in pos_sph to cartesian coordiantes
    """
    function vec_spherical2cartesian(pos_sph, vec_sph)
        r = pos_sph[1]
        theta = pos_sph[2]
        phi = pos_sph[3]

        v_r = vec_sph[1]
        v_theta = vec_sph[2]
        v_phi = vec_sph[3]

        v_x = v_r * sin(theta) * cos(phi) + v_theta * cos(theta) * cos(phi) - v_phi * sin(phi)
        v_y = v_r * sin(theta) * sin(phi) + v_theta * cos(theta) * sin(phi) + v_phi * cos(phi)
        v_z = v_r * cos(theta) - v_theta * sin(theta)
        return [v_x, v_y, v_z]
    end    

    """
    Converts cartesian cordinates to spherical ones...
    """
    function cartesian2spherical(cartesian)
        x = cartesian[1]
        y = cartesian[2]
        z = cartesian[3]
        r = sqrt(x^2 + y^2 + z^2)
        if z > 0
            theta = atan(sqrt(x^2 + y^2) / z)
        elseif z < 0
            theta = pi + atan(sqrt(x^2 + y^2) / z)
        elseif (z == 0) && (x*y != 0)
            theta = pi
        else
            theta = 0 # undefined changed to zero
        end
        if x > 0
            phi = atan(y / x)
        elseif (x < 0) && (y >=0)
            phi = atan(y / x) + pi
        elseif (x < 0) && (y < 0)
            phi = atan(y / x) - pi
        elseif (x == 0) && (y > 0)
            phi = pi
        elseif (x == 0) && (y < 0)
            phi = -pi
        elseif (x == 0) && (y ==0) # undefined changed to zero
            phi = 0
        end
        return [r, theta, phi]
    end


end # module end

