module NoName

    include("modules/functions.jl")
    include("modules/plot.jl") 

    mutable struct Pulsar
        r # pulsar radiuis in [m]
        p # pulsar period in [s]
        pdot # pulsar period derivative in [s/s]
        r_lc # light cylinder radius in [m]
        alpha # inclination angle in [deg.]
        magnetic_axis # in spherical coordinates
        rotation_axis # in spherical coordinates

        function Pulsar()
            r = 10000 # 10 km in merters
            p = 1 # period in seconds
            pdot = 1e-15 # period derivative in s/s
            r_lc = Functions.rlc(p, pdot)
            alpha = 30 # 30 deg by default
            return new(r, p, pdot, r_lc, alpha, (r, 0, 0), (r, deg2rad(alpha), 0))
        end
    end


    function main()

        psr = Pulsar()
        println(fieldnames(Pulsar))

        Plot.pulsar(psr)
        println("Bye")
    end


end # module end

NoName.main()