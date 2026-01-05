module NoName

    include("modules/functions.jl")
    include("modules/plot.jl") 
    include("modules/field.jl")
    include("modules/lines.jl")
    include("modules/sparks.jl")
    include("modules/signal.jl")


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
        grid # grid at the polar cap to calculate potential # multiple grids for simulation
        potential
        electric_field
        drift_velocity
        pot_minmax # what is this? do we need it here?
        sparks_locations # locations in simulation # locations in drift2
        sparks_velocity # step in simulation)
        potential_simulation # potential for simulation step
        spark_radius # spark radius in meters
        line_of_sight # line of sight points
        r_em # emission height
        beta # impact parameter
        los_lines # magnetic lines defined by the line of sight points
        signal # radio intensity
        function Pulsar()
            r = 10_000 # 10 km in merters
            p = 1 # period in seconds
            pdot = 1e-15 # period derivative in s/s
            r_pc = Functions.rdp(p, r)            
            r_lc = Functions.rlc(p)
            alpha = 30 # 30 deg by default
            magnetic_axis = (r, 0, 0)
            rotation_axis = (r, deg2rad(alpha), 0)
            fields = Field.Test() # using test class for now
            fields.beq = Field.beq(p, pdot)
            polar_caps = nothing
            pc = nothing
            open_lines = [[], []]
            sparks = nothing
            grid = nothing
            potential = nothing
            electric_field = nothing
            drift_velocity = nothing
            pot_minmax = nothing
            sparks_locations = []
            sparks_velocity = nothing
            potential_simulation = []
            spark_radius = 15
            line_of_sight = nothing
            r_em = 500_000 # TODO 20 km for tests change it to 500km
            beta = 0 # -0.3 # deg by default
            los_lines = Vector{Vector{Vector{Float64}}}() # zamiast [], szybsze
            signal = nothing
            return new(r, p, pdot, r_pc, r_lc, alpha, magnetic_axis, rotation_axis, fields, polar_caps, pc, open_lines, sparks, grid, potential, electric_field, drift_velocity, pot_minmax, sparks_locations, sparks_velocity, potential_simulation, spark_radius, line_of_sight, r_em, beta, los_lines, signal)
        end
    end


    function full_grid()
        psr = Pulsar()
        Field.calculate_dipole!(psr)
        Field.generate_lines!(psr)
        Lines.calculate_polarcaps!(psr)
        Lines.generate_open!(psr)

        #Sparks.random_sparks!(psr) # cannot calculate potential (points beyond grid. do not use it, just for show) 
        Sparks.create_grid!(psr)
        Sparks.random_sparks_grid!(psr)
        
        Sparks.calculate_potential!(psr)

        Lines.init_line_of_sight(psr)
        Lines.calculate_line_of_sight(psr)
        Plot.pulsar(psr)
        #Plot.potential2D(psr)
        #Plot.potential2Dv2(psr)
    end

    function small_grids()
        psr = Pulsar()
        Field.calculate_dipole!(psr)
        Field.generate_lines!(psr)
        Lines.calculate_polarcaps!(psr)
        Lines.generate_open!(psr)

        #Sparks.random_sparks!(psr) 
        Sparks.init_sparks1!(psr ;num=5)
        #Sparks.init_sparks2!(psr ;num=5)
        #Sparks.init_sparks3!(psr ;num=10, rfmax=0.7)
        
        # 3D simulation starts here 
        Sparks.create_grids!(psr)
        Sparks.calculate_potentials!(psr) # calculates step in sparks_velocity
        
        #Plot.small_grids(psr) # plots small grids
        Plot.steps(psr) # moving sparks

       
    end

    function full_plus_smallgrids()
        psr = Pulsar()
        #Field.calculate_dipole!(psr)
        #Field.generate_lines!(psr)
        Lines.calculate_polarcaps!(psr)
        #Lines.generate_open!(psr)

        #Sparks.random_sparks!(psr) 
        Sparks.init_sparks1!(psr ;num=5)
        #Sparks.init_sparks2!(psr ;num=5)
        #Sparks.init_sparks3!(psr ;num=30, rfmax=0.7)

        #Sparks.generate_potentials # TODO
        Sparks.simulate_sparks(psr)
        Plot.steps2D(psr)
    end

    function generate_signal()
        psr = Pulsar()

        Lines.calculate_polarcaps!(psr)

        #Field.calculate_dipole!(psr)

        Lines.init_line_of_sight(psr, num=100)
        Lines.calculate_line_of_sight(psr)

        Sparks.init_sparks1!(psr ;num=5)

        Sparks.simulate_sparks(psr; n_steps=1000)
        Signal.generate(psr)
        Plot.signal(psr)
        
    end

    function main()

        #full_grid()
        #small_grids()
        #full_plus_smallgrids()

        generate_signal()

        println("Bye")
    end


end # module end

NoName.main()