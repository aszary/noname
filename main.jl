module NoName

    include("modules/functions.jl")
    include("modules/plot.jl") 
    include("modules/field.jl") 

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
        polar_cap # polar cap radius in spherical coordinates

        

        function Pulsar()
            r = 10000 # 10 km in merters
            p = 1 # period in seconds
            pdot = 1e-15 # period derivative in s/s
            r_lc = Functions.rlc(p, pdot)
            alpha = 30 # 30 deg by default
            fields = Field.Test() # using Test struct from Field module
            r_pc = Functions.rdp(p, r)
            polar_cap = []
            return new(r, p, pdot, r_lc, alpha, (r, 0, 0), (r, deg2rad(alpha), 0), fields, r_pc, polar_cap)
        end
    end


    function main()

        psr = Pulsar()
        Field.calculate_dipole!(psr)
        Field.generate_vacuum!(psr, step=100, stepsnum=20000)
        Field.generate_polar_cap!(psr)
        Field.generate_magnetic_lines_from_polar_cap!(psr; num_lines_from_cap=50, step=100, stepsnum=20000)

        println(fieldnames(Pulsar))
        println(psr.r_lc / 1e3, " km")
        println("Polar cap radius: ", psr.r_pc / 1e3, " km")

        Plot.pulsar(psr)
        println("Bye")
    end


end # module end

NoName.main()