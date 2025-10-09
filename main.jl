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
        potential
        electric_field
        drift_velocity 

        function Pulsar()
            r = 10000 # 10 km in merters
            p = 1 # period in seconds
            pdot = 1e-15 # period derivative in s/s
            r_lc = Functions.rlc(p)
            alpha = 30 # 30 deg by default
        
            magnetic_axis = (r, 0, 0) # in spherical coordinates
            rotation_axis= (r, deg2rad(alpha), 0) # in spherical coordinates
            pc= nothing
            r_pc = Functions.rdp(p, r)
            fields = Field.Test() # usingt testt class for now
            grid = nothing
            sparks = nothing
            potential = nothing
            electric_field = nothing
            drift_velocity = nothing
            return new(r, p, pdot, r_lc, alpha, magnetic_axis, rotation_axis, pc, r_pc, fields, grid, sparks, potential, electric_field, drift_velocity)
        end
    end


    function main()

        psr = Pulsar()
        Field.calculate_dipole!(psr)
        Field.generate_lines!(psr)
        Field.pc(psr; phi_num=10)
        #Field.calculate_polarcap!(psr)
        Field.generate_polarcap_lines!(psr)
        Sparks.create_grid!(psr; size=50)
        Sparks.random_sparks_grid!(psr; min_dist=20, trials=10000)
        Sparks.calculate_potential!(psr)
        #println(fieldnames(Pulsar))
        #println(psr.r_lc / 1e3, " km")
        Plot.pulsar(psr)
        println("Bye")
    end


end # module end

NoName.main()