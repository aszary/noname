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
        
        

        function Pulsar()
            r = 10000 # 10 km in merters
            p = 1 # period in seconds
            pdot = 1e-15 # period derivative in s/s
            r_lc = Functions.rlc(p)
            alpha = 30 # 30 deg by default
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
            return new(r, p, pdot, r_lc, alpha, magnetic_axis, rotation_axis, fields, r_pc, pc, polar_cap, grid, sparks, potential, electric_field, drift_velocity, pot_minmax, sparks_locations, sparks_velocity, potential_simulation)
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
        Field.generate_vacuum!(psr, step=100, stepsnum=20000)
        Field.pc(psr; phi_num=10)
        Field.generate_polar_cap!(psr)
        Field.generate_magnetic_lines_from_polar_cap!(psr; num_lines_from_cap=50, step=100, stepsnum=20000)

        #Sparks.random_sparks!(psr) 
        Sparks.init_sparks1!(psr ;num=5)
        #Sparks.init_sparks2!(psr ;num=5)
        #Sparks.init_sparks3!(psr ;num=10, rfmax=0.7)

        Sparks.simulate_sparks(psr)
        #Sparks.calculate_potential_custom!(psr)

        Plot.steps2DD(psr)
        #Plot.steps(psr)
    end

    function test()

        psr = Pulsar()
        Field.calculate_dipole!(psr)
        Field.generate_vacuum!(psr, step=100, stepsnum=20000)
        Field.pc(psr; phi_num=10)
        Field.generate_polar_cap!(psr)
        Field.generate_magnetic_lines_from_polar_cap!(psr; num_lines_from_cap=50, step=100, stepsnum=20000)
        Sparks.create_grid!(psr; size=50)
        Sparks.random_sparks_grid!(psr; min_dist=20, trials=10000)

        #Sparks.calculate_potential!(psr)

        Plot.pulsar(psr)
        #Plot.potential2D(psr)
        #Plot.potential2Dv2(psr) #nie działa
    end


    function main()
        #test()
        #small_gridS()
        #full_grid()
        full_plus_smallgrids()

    end


end # module end


NoName.main()