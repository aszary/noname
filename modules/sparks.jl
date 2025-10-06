module Sparks
    using LinearAlgebra
    include("functions.jl")


    """

    Random sparks at the polar cap (grids will be generated later!)

    # Arguments

    - min_dist: minimum distant in meters

    """
    function random_sparks!(psr; min_dist=40, trials=1000)
        sp = []

        # maximum theta
        thm = asin((psr.r_pc - min_dist) / psr.r) # plus distance from pc boundry condition
        for i in 1:trials
            phi = rand() * 2pi
            theta = rand() * thm
            car = Functions.spherical2cartesian([psr.r, theta, phi])
            md = 2 * min_dist
            for s in sp
                dist = norm([s[1], s[2], s[3]] - [car[1], car[2], car[3]])
                if dist < md
                    md = dist
                end
            end
            if md > min_dist
                push!(sp, [car[1], car[2], car[3]])
            end
        end
        psr.sparks = sp
        println("Number of sparks added: ", size(sp)[1])
    end







end # module end