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

end # module end

