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
        dalpha
        magnetic_axis # in spherical coordinates
        rotation_axis # in spherical coordinates
        fields # magnetic and electric fields
        r_pc # polar cap radius [in m]
        pc
        polar_cap # polar cap radius in spherical coordinates
        grid
        sparks # sparks positions in cartesian coordinates
        potential
        electric_field
        drift_velocity
        pot_minmax # what is this? do we need it here?
        sparks_locations # locations in simulation # locations in drift2
        sparks_velocity # step in simulation)
        potential_simulation
        line_of_sight 
        r_em # emission height in meters
        beta
        signal
        spark_radius
        pulses
       

        
        

        function Pulsar()
            r = 10000 # 10 km in merters
            p = 1 # period in seconds
            pdot = 1e-15 # period derivative in s/s
            r_lc = Functions.rlc(p)
            alpha = 30 # 30 deg by default
            dalpha = deg2rad(5)
            magnetic_axis = (r, 0, 0) # in spherical coordinates
            rotation_axis= (r, deg2rad(alpha), 0) # in spherical coordinates
            fields = Field.Test() # using Test struct from Field module
            r_pc = Functions.rdp(p, r)
            pc = nothing
            polar_cap = []
            grid = nothing
            sparks = nothing
            potential = nothing
            electric_field = nothing
            drift_velocity = nothing
            pot_minmax = nothing
            sparks_locations = []
            sparks_velocity = nothing
            potential_simulation = []
            line_of_sight = nothing
            r_em = 500_000 # emission height in meters
            beta = 1 # impact angle in deg.
            signal = nothing
            spark_radius = 20
            pulses = nothing
            return new(r, p, pdot, r_lc, alpha, dalpha, magnetic_axis, rotation_axis, fields, r_pc, pc, polar_cap, grid, sparks, potential, electric_field, drift_velocity, pot_minmax, sparks_locations, sparks_velocity, potential_simulation, line_of_sight, r_em, beta, signal, spark_radius, pulses)
        end
    end



    function full_grid()

        psr = Pulsar()
        Field.calculate_dipole!(psr)
        Field.generate_vacuum!(psr, step=100, stepsnum=20000)
        Field.pc(psr; phi_num=10)
        Field.generate_polar_cap!(psr)
        Field.generate_magnetic_lines_from_polar_cap!(psr; num_lines_from_cap=50, step=100, stepsnum=20000)
        Sparks.create_grid!(psr; size=50)
        Sparks.random_sparks_grid!(psr; min_dist=20, trials=10000)

        Sparks.calculate_potential!(psr)

        Plot.pulsar(psr)
        #Plot.potential2D(psr)
        #Plot.potential2Dv2(psr) #nie działa
    end

    function small_gridS()

        psr = Pulsar()
        #Field.calculate_dipole!(psr)
        #Field.generate_vacuum!(psr, step=100, stepsnum=20000)
        Field.pc(psr; phi_num=10)
        Field.generate_polar_cap!(psr)
        #Field.generate_magnetic_lines_from_polar_cap!(psr; num_lines_from_cap=50, step=100, stepsnum=20000)

        #Sparks.create_grid!(psr; size=50)

        #Sparks.random_sparks_grid!(psr; min_dist=20, trials=10000)
        #Sparks.random_sparks!(psr; num=10, min_dist=20, trials=1000)

        Sparks.init_sparks1!(psr; rfs=[0.295, 0.5], num=3)
        #Sparks.init_sparks2!(psr; rfs=[0.7, 0.9], num=3)
        #Sparks.init_sparks3!(psr; num=10, rfmax=0.7)


        #simulation
        Sparks.create_grids!(psr)
        Sparks.calculate_potentials!(psr)
        
        
        Plot.steps(psr)
        #Plot.potential2D(psr)
        #Plot.potential2Dv2(psr) 
    end

    function full_plus_smallgrids()
        psr = Pulsar()
        Field.calculate_dipole!(psr)
        Field.generate_vacuum!(psr, step=100)
        Field.pc(psr; phi_num=10)
        Field.generate_polar_cap!(psr)
        Field.generate_magnetic_lines_from_polar_cap!(psr; num_lines_from_cap=50, step=100, stepsnum=20000)

        #Sparks.random_sparks!(psr) 
        Sparks.init_sparks1!(psr ;num=5)
        #Sparks.init_sparks2!(psr ;num=5)
        #Sparks.init_sparks3!(psr ;num=10, rfmax=0.7)

        Sparks.simulate_sparks(psr)
        #Sparks.calculate_potential_custom!(psr)
        Plot.init_line_of_sight!(psr; thresh=Functions.theta_max(psr.r_pc, psr), w=deg2rad(10))
        Plot.steps2DD(psr)
        #Plot.steps(psr)
    end

    function test()

        psr = Pulsar()
        Field.calculate_dipole!(psr)
        Field.generate_vacuum!(psr, step=100)
        Field.pc(psr; phi_num=10)
        Field.generate_polar_cap!(psr)
        Field.generate_magnetic_lines_from_polar_cap!(psr; num_lines_from_cap=40, step=100, stepsnum=20000)
        Lines.init_line_of_sight!(psr)
        Field.generate_magnetic_lines_from_los!(psr)
        Sparks.create_grid!(psr; size=50)
        Sparks.random_sparks_grid!(psr; min_dist=20, trials=10000)
        Sparks.calculate_potential!(psr)
        
        Signal.signal_from_sparks(psr; sigma=50.0)
        
        Plot.steps_with_los_signal(psr)
        #Plot.pulsar2(psr)

        #Plot.potential2D(psr)
        #Plot.potential2Dv2(psr) #nie działa
    end

     function generate_signal()
        
        psr = Pulsar()
        Field.calculate_dipole!(psr)
        Field.generate_vacuum!(psr, step=100)
        Field.pc(psr; phi_num=10)
        Field.generate_polar_cap!(psr)
        Field.generate_magnetic_lines_from_polar_cap!(psr; num_lines_from_cap=40, step=100, stepsnum=20000)
        Sparks.create_grid!(psr; size=50)
        Sparks.init_sparks1!(psr; rfs=[0.3, 0.6], num=6, center=true)
        Sparks.simulate_sparks(psr, n_steps=500, skip_steps=10, speedup=10)
        #Sparks.random_sparks_grid!(psr; min_dist=20, trials=10000)
        #Sparks.calculate_potential!(psr)
        Lines.init_line_of_sight!(psr)
        Field.generate_magnetic_lines_from_los!(psr)
        
        Signal.signal_from_sparks(psr)
        Plot.signal(psr)
     end
        
        
        



    function main()
        generate_signal()
        #test()
        #small_gridS()
        #full_grid()
        #full_plus_smallgrids()
    end


end # module end


NoName.main()