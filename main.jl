module NoName

    include("modules/functions.jl")
    include("modules/plot.jl")
    include("modules/field.jl")
    include("modules/sparks.jl")

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
        line_of_sight #line of sight points
        r_em #emission height
        beta #impact parameter
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
            line_of_sight = nothing
            r_em = 500000
            beta = deg2rad(1)
            return new(r, p, pdot, r_lc, alpha, magnetic_axis, rotation_axis, pc, r_pc, fields, grid, sparks, locations, sparks_velocities, potential, pot_minmax, electric_field, drift_velocity, sparks_locations, sparks_velocities, potential_simulation, line_of_sight, r_em, beta)
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

    function main()
        full_grid()
        #small_grid()
        #full_plus_small_grid()
    end


end # module end

NoName.main()