module NoName

    include("modules/functions.jl")
    include("modules/plot.jl") 
    include("modules/field.jl")
    include("modules/lines.jl")
    include("modules/sparks.jl")


    mutable struct Pulsar
        r # pulsar radiuis in [m]
        p # pulsar period in [s]
        pdot # pulsar period derivative in [s/s]
        r_pc # polar cap radius [in m]       
        r_lc # light cylinder radius in [m]
        alpha # inclination angle in [deg.]
        magnetic_axis # in spherical coordinates
        rotation_axis # in spherical coordinates
        fields # magnetic and electric fields
        polar_caps # two polar caps boundries (xs, ys, zs)
        pc # single polar cap
        open_lines # magnetic lines at polar cap boundries
        sparks # sparks locations
        grid # grid at the polar cap to calculate potential
        function Pulsar()
            r = 10000 # 10 km in merters
            p = 1 # period in seconds
            pdot = 1e-15 # period derivative in s/s
            r_pc = Functions.rdp(p, r)            
            r_lc = Functions.rlc(p)
            alpha = 30 # 30 deg by default
            magnetic_axis = (r, 0, 0)
            rotation_axis = (r, deg2rad(alpha), 0)
            fields = Field.Test() # using test class for now
            polar_caps = nothing
            pc = nothing
            open_lines = [[], []]
            sparks = nothing
            grid = nothing
            return new(r, p, pdot, r_pc, r_lc, alpha, magnetic_axis, rotation_axis, fields, polar_caps, pc, open_lines, sparks, grid)
        end
    end


    function main()

        psr = Pulsar()
        Field.calculate_dipole!(psr)
        Field.generate_lines!(psr)
        Lines.calculate_polarcaps!(psr)
        Lines.generate_open!(psr)
        #Sparks.random_sparks!(psr) # cannot calculate potential (points beyond grid. do not use it, just for show) 
        Sparks.create_grid!(psr)
        Sparks.random_sparks_grid!(psr)
        
        #Sparks.calculate_potential!(psr)
        #println(fieldnames(Pulsar))
        #println(psr.r_lc / 1e3, " km")

        Plot.pulsar(psr)
        println("Bye")
    end


end # module end

NoName.main()