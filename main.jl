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
        beta # impact parameter in [deg.]
        los_lines # magnetic lines defined by the line of sight points
        signal # radio intensity for continous signal
        pulses # single pulses generated from signal
        longitudes # single pulse longitudes
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
            spark_radius = 20
            line_of_sight = nothing
            r_em = 500_000 # TODO 20 km for tests change it to 500km
            beta = 3.0 # deg by default
            los_lines = Vector{Vector{Vector{Float64}}}() # instead [], faster
            signal = nothing
            pulses = nothing
            longitudes = nothing
            return new(r, p, pdot, r_pc, r_lc, alpha, magnetic_axis, rotation_axis, fields, polar_caps, pc, open_lines, sparks, grid, potential, electric_field, drift_velocity, pot_minmax, sparks_locations, sparks_velocity, potential_simulation, spark_radius, line_of_sight, r_em, beta, los_lines, signal, pulses, longitudes)
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
        Plot.steps(psr) # moving sparks not moving? repair..

       
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

        Lines.init_line_of_sight(psr, num=200)
        Lines.calculate_line_of_sight(psr)

        # TODO work on skip_steps to ensure only single pulses
        #Sparks.init_sparks1!(psr ;num=5)
        #Sparks.simulate_sparks(psr; n_steps=2000, skip_steps=20, speedup=10)
        #Sparks.save_sparks(psr; num=1)

        Sparks.load_sparks(psr; num=1)
        Signal.generate_signal(psr; noise_level=0.05)
        Signal.generate_pulses(psr)

        
        # CHECKING
        x,y,z = psr.los_lines[end][1][1], psr.los_lines[end][2][1], psr.los_lines[end][3][1]
        sph = Functions.cartesian2spherical([x,y,z])
        rho = rad2deg(Signal.rho_from_theta(sph[2]))
        println("rho: $rho")
        println("theta=", rad2deg(sph[2]))
        W = Signal.pulse_width_deg(psr.alpha, psr.beta, rho)
        println("alpha = $(psr.alpha) beta = $(psr.beta)")
        println("Szerokość profilu teoretyczna: $(round(W, digits=2))°")
        println("Szerokość profilu z symukacji ", psr.longitudes[end]-psr.longitudes[1], " deg.")

        #Plot.signal(psr)
        Plot.pulses(psr)
        Plot.pulses0(psr)
        Plot.pulses1(psr)
        
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