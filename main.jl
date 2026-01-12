module NoName

    include("modules/functions.jl")
    include("modules/plot.jl")
    include("modules/field.jl")
    include("modules/sparks.jl")
    include("modules/transformations.jl")
    include("modules/lines.jl")
    include("modules/signal.jl")
    

    mutable struct Pulsar
        r # pulsar radiuis in [m]
        p # pulsar period in [s]
        pdot # pulsar period derivative in [s/s]
        r_lc # light cylinder radius in [m]
        alpha # inclination angle in [deg.]
        magnetic_axis # in spherical coordinates
        rotation_axis # in spherical coordinates
        pc #polar cap in spherical coordinates
        r_pc # polar cap radius [in m]
        fields  # magnetic and electric fields
        grid #idk od oli
        sparks # sparks positions in cartesian coordinates
        locations # for the simulation
        sparks_velocities # for the simulation
        potential
        pot_minmax
        electric_field
        drift_velocity 
        sparks_locations # locations in simulation # locations in drift2
        sparks_velocity # step in simulation)
        potential_simulation # potential for simulation step
        line_of_sight #line of sight points at the polar cap
        r_em #emission height
        beta #impact parameter
        signal # radio intensity
        spark_radius # spark radius in meters
        pulses
        function Pulsar()
            r = 10000 # 10 km in merters
            p = 1 # period in seconds
            pdot = 1e-15 # period derivative in s/s
            r_lc = Functions.rlc(p)
            alpha = 30 # 30 deg by default
        
            magnetic_axis = (r, 0, 0) # in spherical coordinates
            rotation_axis= (r, deg2rad(alpha), 0) # in spherical coordinates
            pc= nothing # 
        
            r_pc = Functions.rdp(p, r)
            fields = Field.Test() # usingt testt class for now
            grid = nothing
            sparks = nothing
            locations = []
            sparks_velocities = [] 
            potential = nothing
            pot_minmax = nothing
            electric_field = nothing
            drift_velocity = nothing
            sparks_locations = []
            sparks_velocity = nothing
            potential_simulation = []
            line_of_sight = Vector{Vector{Vector{Float64}}}()
            r_em = 500000
            beta = 1
            signal = nothing
            spark_radius = 20
            pulses = nothing
            return new(r, p, pdot, r_lc, alpha, magnetic_axis, rotation_axis, pc, r_pc, fields, grid, sparks, locations, sparks_velocities, potential, pot_minmax, electric_field, drift_velocity, sparks_locations, sparks_velocities, potential_simulation, line_of_sight, r_em, beta, signal, spark_radius, pulses)
        end
    end
    function full_grid()
        psr = Pulsar()
        Field.calculate_dipole!(psr)
        #Field.generate_lines!(psr)
        #Field.calculate_polarcap!(psr)
        Field.pc(psr; phi_num=100)
        Field.generate_polarcap_lines!(psr)
        Sparks.create_grid!(psr; size=500)
        Sparks.random_sparks_grid!(psr; min_dist=20, trials=10)
        Sparks.calculate_potential!(psr)
        #println(fieldnames(Pulsar))
        #println(psr.r_lc / 1e3, " km")
        #Plot.potential2D(psr)
        #Plot.potential2Dv2(psr)
        Plot.pulsar(psr)
        println("Bye")
    end
    function small_grid()
        psr = Pulsar()
        #Field.calculate_dipole!(psr)
        #Field.generate_lines!(psr)
        #Field.calculate_polarcap!(psr)
        Field.pc(psr; phi_num=100)
        #Field.generate_polarcap_lines!(psr)
        #Sparks.create_grid!(psr; size=100)
        #Sparks.random_sparks_grid!(psr; min_dist=20, trials=20)
        #Sparks.calculate_potential!(psr)
        #println(fieldnames(Pulsar))
        #println(psr.r_lc / 1e3, " km")
        Sparks.init_sparks1!(psr, num=5)
        #Sparks.init_sparks2!(psr)
        #Sparks.init_sparks3!(psr)
        Sparks.create_grids!(psr)
        Sparks.calculate_potentials!(psr)
        #Sparks.calculate_potential_sparks!(psr)
        #Plot.potential2D(psr)
        #Plot.potential2Dv2(psr)
        Plot.steps(psr)
        #Plot.steps2D(psr)
        #Plot.small_grid(psr)
        #Plot.pulsar(psr)
        println("Bye")
    end
    function full_plus_small_grid()
        psr = Pulsar()
        Field.calculate_dipole!(psr)
        Field.generate_lines!(psr)
        Field.pc(psr; phi_num=100)
        Sparks.init_sparks1!(psr, num=5)
        #Sparks.create_grids!(psr)
        #Sparks.calculate_potentials!(psr)
        #Field.generate_lines!(psr)
        #Field.calculate_polarcap!(psr)
        Sparks.simulate_sparks!(psr)
        #Field.generate_polarcap_lines!(psr)
        #Sparks.create_grid!(psr; size=100)
        #Sparks.random_sparks_grid!(psr; min_dist=20, trials=20)
        #Sparks.calculate_potential!(psr)
        #println(fieldnames(Pulsar))
        #println(psr.r_lc / 1e3, " km")
        
        #Sparks.init_sparks2!(psr)
        #Sparks.init_sparks3!(psr)
        #Sparks.calculate_potential_sparks!(psr)
        #Plot.potential2D(psr)
        #Plot.potential2Dv2(psr)
        #Plot.steps(psr)
        Plot.steps2D(psr)
        #Plot.small_grid(psr)
        #Plot.pulsar(psr)
        println("Bye")
    end
    function generate_signal()
        psr = Pulsar()
        Field.calculate_dipole!(psr)
        Field.pc(psr; phi_num=100)
        Field.generate_polarcap_lines!(psr)
        Sparks.create_grid!(psr; size=500) # Still needed for potential map calculation
        Sparks.init_sparks1!(psr; rfs=[0.3, 0.6], num=6, center=true)
        Sparks.simulate_sparks!(psr)
        Sparks.calculate_potential!(psr)
        Sparks.create_grids!(psr)
        Sparks.calculate_potentials!(psr; save=true) 
        Lines.calculate_los(psr)
        Signal.generate_signal(psr)
        Signal.generate_pulses(psr)
        #Plot.signal(psr)
        #Plot.pulses(psr)
        Plot.pulses0(psr)

        println("Bye")
    end
    function main()
        generate_signal()
        #full_grid()
        #small_grid()
        #full_plus_small_grid()
    end


end # module end

NoName.main()