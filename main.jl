module NoName

    using JSON3
    include("modules/functions.jl")
    include("modules/plot.jl")
    include("modules/field.jl")
    include("modules/nsfield.jl")
    include("modules/lines.jl")
    include("modules/sparks.jl")
    include("modules/signal.jl")
    include("modules/lbc.jl")


    mutable struct Pulsar
        r # pulsar radius in [m]
        p # pulsar period in [s]
        pdot # pulsar period derivative in [s/s]
        r_pc # polar cap radius [in m]       
        r_lc # light cylinder radius in [m]
        alpha # inclination angle in [deg.]
        magnetic_axis # in spherical coordinates
        rotation_axis # in spherical coordinates
        nsfield # non-dipolar magneitc field structure 
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
        spark_radii # spark radii in meters for e.g. LBC model [n_step][spark_num]?
        line_of_sight # line of sight points
        r_em # emission height
        beta # impact parameter in [deg.]
        los_lines # magnetic lines defined by the line of sight points
        signal # radio intensity for continous signal
        pulses # single pulses generated from signal
        longitudes # single pulse longitudes
        ellipse_fit # ellipse fit to the polar cap points
        function Pulsar()
            r = 10_000 # 10 km in merters
            p = 1 # period in seconds
            pdot = 1e-15 # period derivative in s/s
            r_pc = Functions.rdp(p, r)            
            r_lc = Functions.rlc(p)
            alpha = 30 # 30 deg by default
            magnetic_axis = (r, 0, 0)
            rotation_axis = (r, deg2rad(alpha), 0)
            nsfield = NSField.Field()
            fields = Field.Test() # using test class for now
            fields.beq = Field.beq(p, pdot)
            polar_caps = nothing
            pc = nothing
            open_lines = []
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
            spark_radii = nothing
            line_of_sight = nothing
            r_em = 500_000  # 500 km
            beta = 4.0 # deg by default
            los_lines = Vector{Vector{Vector{Float64}}}() # instead [], faster
            signal = nothing
            pulses = nothing
            longitudes = nothing
            ellipse_fit = nothing
            return new(r, p, pdot, r_pc, r_lc, alpha, magnetic_axis, rotation_axis, nsfield, fields, polar_caps, pc, open_lines, sparks, grid, potential, electric_field, drift_velocity, pot_minmax, sparks_locations, sparks_velocity, potential_simulation, spark_radius, spark_radii, line_of_sight, r_em, beta, los_lines, signal, pulses, longitudes, ellipse_fit)
        end
        function Pulsar(json_file)
            d = JSON3.read(json_file)
            #open("input/test.json", "w") do io
            #    JSON3.pretty(io, JSON3.write(d))
            #end
            r = d.psr.R
            p = d.psr.P0
            pdot = d.psr.PDOT
            alpha = d.psr.alpha
            beta = d.psr.beta
            r_em = d.psr.R_em
            spark_radius = d.psr.spark_radius
            spark_radii = nothing
            r_pc = Functions.rdp(p, r)            
            r_lc = Functions.rlc(p)
            magnetic_axis = (r, 0, 0)
            rotation_axis = (r, deg2rad(alpha), 0)
            nsfield = NSField.Field(d)
            fields = Field.Test() # using test class for now
            fields.beq = Field.beq(p, pdot)
            polar_caps = nothing
            pc = nothing
            open_lines = []
            sparks = nothing
            grid = nothing
            potential = nothing
            electric_field = nothing
            drift_velocity = nothing
            pot_minmax = nothing
            sparks_locations = []
            sparks_velocity = nothing
            potential_simulation = []
            line_of_sight = nothing
            los_lines = Vector{Vector{Vector{Float64}}}() # instead [], faster
            signal = nothing
            pulses = nothing
            longitudes = nothing
            ellipse_fit = nothing
            return new(r, p, pdot, r_pc, r_lc, alpha, magnetic_axis, rotation_axis, nsfield, fields, polar_caps, pc, open_lines, sparks, grid, potential, electric_field, drift_velocity, pot_minmax, sparks_locations, sparks_velocity, potential_simulation, spark_radius, spark_radii, line_of_sight, r_em, beta, los_lines, signal, pulses, longitudes, ellipse_fit)
        end
    end


    function full_grid()
        psr = Pulsar()
        Field.calculate_dipole!(psr)
        Field.generate_lines!(psr)
        Lines.calculate_polarcaps!(psr)
        Lines.generate_open_obsolete!(psr)

        #Sparks.random_sparks!(psr) # cannot calculate potential (points beyond grid. do not use it, just for show) 
        Sparks.create_grid!(psr)
        Sparks.random_sparks_grid!(psr)
        
        Sparks.calculate_potential!(psr)

        Lines.init_line_of_sight(psr)
        Lines.calculate_line_of_sight_dipole(psr)
        Plot.pulsar(psr)
        #Plot.potential2D(psr)
        #Plot.potential2Dv2(psr)
    end

    function small_grids()
        psr = Pulsar()
        Field.calculate_dipole!(psr)
        Field.generate_lines!(psr)
        Lines.calculate_polarcaps!(psr)
        Lines.generate_open_obsolete!(psr)

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
        #Lines.generate_open_obsolete!(psr)

        #Sparks.random_sparks!(psr) 
        Sparks.init_sparks1!(psr ;num=5)
        #Sparks.init_sparks2!(psr ;num=5)
        #Sparks.init_sparks3!(psr ;num=30, rfmax=0.7)

        #Sparks.generate_potentials # TODO
        Sparks.simulate_sparks_mc(psr)
        Plot.steps2D(psr)
    end


    function generate_signal_dipole()
        psr = Pulsar("input/1.json") # works only with no anomalies!!!

        Lines.calculate_polarcaps!(psr)

        #Field.calculate_dipole!(psr)

        Lines.init_line_of_sight(psr, num=100)
        Lines.calculate_line_of_sight_dipole(psr)

        Lines.generate_open!(psr, num=10)

        # TODO work on n_steps + save_every for single pulses
        #Sparks.init_sparks1!(psr ;num=5)
        #Sparks.simulate_sparks_mc(psr; n_steps=2000, save_every=20, speedup=10)
        #Sparks.simulate_sparks_solidbody(psr; n_steps=100)
        Sparks.simulate_sparks_lbc(psr; n_steps=500, co_angl=-90.0, h_drft=0.5, save_every=1)
        Sparks.save_sparks(psr; num=2)

        #Plot.sparks(psr)
        Sparks.load_sparks(psr; num=2)

        #Signal.generate_signal(psr; noise_level=0.05) # old same sizes!
        Signal.generate_signal_radii(psr; noise_level=0.05) # new
        Signal.generate_pulses(psr, pulse_max=500)
        
        Plot.signal(psr)
        #Plot.pulses(psr, number=500)
        #Plot.pulses0(psr)
        #Plot.pulses1(psr)
        
    end


    function generate_signal()
        #psr = Pulsar("input/1.json")
        psr = Pulsar("input/3.json")

        Lines.init_line_of_sight(psr, num=100)
        Lines.calculate_line_of_sight(psr)

        Lines.generate_open!(psr, num=20)

        # TODO work on init sparks for non-dipolar polar cap!
        #Sparks.init_sparks1!(psr ;num=5) # dipolar 
        Sparks.init_sparks1_ellipse!(psr ;rfs=[0.2, 0.5, 0.79], num=3) # non-dipolar 
        # TODO work on n_steps + save_every for single pulses
        #Sparks.simulate_sparks_mc(psr; n_steps=2000, save_every=40, speedup=0.1)
        Sparks.simulate_sparks_solidbody(psr; n_steps=500)
        #Sparks.simulate_sparks_lbc(psr; n_steps=5000, co_angl=90.0, h_drft=0.1, save_every=10)
        Sparks.save_sparks(psr; num=2)

        #Plot.sparks(psr)
        Sparks.load_sparks(psr; num=2)

        #Signal.generate_signal(psr; noise_level=0.05) # old same sizes!
        Signal.generate_signal_radii(psr; noise_level=0.05) # new
        Signal.generate_pulses(psr, pulse_max=500)
        
        Plot.signal(psr)
        #Plot.pulses(psr, number=500)
        #Plot.pulses0(psr)
        #Plot.pulses1(psr)
        
    end


    function model_field()
        #psr = Pulsar("input/1.json")
        psr = Pulsar("input/2.json")
        #psr = Pulsar("input/3.json")

        Lines.init_line_of_sight(psr, num=5)
        Lines.calculate_line_of_sight(psr)

        Lines.generate_open!(psr, num=10)

        #println(psr.nsfield)

        #Plot.anomalies(psr)
        #Plot.anomalies2D(psr)
        Plot.polar_cap2D(psr)

       
    end



    function main()

        #full_grid()
        #small_grids()
        #full_plus_smallgrids()

        #generate_signal_dipole()
        #generate_signal()

        model_field()

        println("Bye")
    end


end # module end

NoName.main()