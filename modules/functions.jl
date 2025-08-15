module Functions 
    using PhysicalConstants.CODATA2018 # in SI units


    """
    Calculates rlc in meters
    """
    function rlc(p, pdot)
        c = SpeedOfLightInVacuum.val # no units hereafter
        return c * p / (2 * pi)
    end


end # module end

