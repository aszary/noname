module NoName

    include("modules/functions.jl")
    include("modules/plot.jl") 
    include("modules/field.jl")
    include("modules/lines.jl")

    mutable struct Pulsar
        r # pulsar radiuis in [m]
        p # pulsar period in [s]
        pdot # pulsar period derivative in [s/s]
        r_lc # light cylinder radius in [m]
        alpha # inclination angle in [deg.]
        magnetic_axis # in spherical coordinates
        rotation_axis # in spherical coordinates
        fields # magnetic and electric fields
        polar_caps # polar caps boundries (xs, ys, zs)
        open_lines # magnetic lines at polar cap boundries
        function Pulsar()
            r = 10000 # 10 km in merters
            p = 1 # period in seconds
            pdot = 1e-15 # period derivative in s/s
            r_lc = Functions.rlc(p)
            alpha = 30 # 30 deg by default
            fields = Field.Test() # using test class for now
            return new(r, p, pdot, r_lc, alpha, (r, 0, 0), (r, deg2rad(alpha), 0), fields, nothing, [[], []])
        end
    end


    function main()

        psr = Pulsar()
        Field.calculate_dipole!(psr)
        Field.generate_lines!(psr)
        Lines.calculate_polarcaps!(psr)
        Lines.generate_open!(psr)
        #println(fieldnames(Pulsar))
        #println(psr.r_lc / 1e3, " km")

        Plot.pulsar(psr)
        println("Bye")
    end


end # module end

NoName.main()