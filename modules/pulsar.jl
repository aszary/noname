    const DEFAULT_SPARKS_CONFIG = (
        model     = "solidbody",
        mc        = (n_steps = 2000, save_every = 40, speedup = 10.1),
        lbc       = (co_angl = 0.0,),
        init      = (method = "ellipse", rfs = [0.2, 0.5, 0.79], num = 3, spacing = "t"),
    )

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
        pa # position angle
        stokes_q # Stokes Q [npulse × nbins], same shape as signal
        stokes_u # Stokes U [npulse × nbins]
        stokes_v # Stokes V [npulse × nbins], model: V ∝ dI/dφ per pulse
        pulses # single pulses generated from signal
        longitudes # single pulse longitudes
        ellipse_fit # ellipse fit to the polar cap points
        p3 # drift repetation time, one value per pulse (Vector{Float64}, length npulse); constant P3 is stored as a filled vector
        npulse # number of single pulses
        noise_level # noise level in single pulses
        output_num # output directory number for save_sparks/load_sparks
        sparks_config # spark simulation model and its parameters
        amplitudes # spark amplitudes
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
            pa = nothing
            stokes_q = nothing
            stokes_u = nothing
            stokes_v = nothing
            pulses = nothing
            longitudes = nothing
            ellipse_fit = nothing
            npulse = 500
            p3 = fill(10.0, npulse)
            noise_level = 0.05
            output_num = 1
            sparks_config = DEFAULT_SPARKS_CONFIG
            amplitudes = nothing
            return new(r, p, pdot, r_pc, r_lc, alpha, magnetic_axis, rotation_axis, nsfield, fields, polar_caps, pc, open_lines, sparks, grid, potential, electric_field, drift_velocity, pot_minmax, sparks_locations, sparks_velocity, potential_simulation, spark_radius, spark_radii, line_of_sight, r_em, beta, los_lines, signal, pa, stokes_q, stokes_u, stokes_v, pulses, longitudes, ellipse_fit, p3, npulse, noise_level, output_num, sparks_config, amplitudes)
        end
        function Pulsar(json_file)
            d = JSON3.read(read(json_file, String))
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
            pa = nothing
            stokes_q = nothing
            stokes_u = nothing
            stokes_v = nothing
            pulses = nothing
            longitudes = nothing
            ellipse_fit = nothing
            npulse = d.psr.npulse

            # P3 can be given as a single constant value, or as an array of
            # values (one per pulse/period) to model a changing P3. Either way
            # it is stored as a Vector{Float64} of length npulse.
            raw_p3 = d.psr.P3
            p3 = isa(raw_p3, AbstractArray) ? Float64.(collect(raw_p3)) : fill(Float64(raw_p3), npulse)

            noise_level = d.psr.noise_level
            output_num = d.psr.output_num
            sparks_config = haskey(d, :sparks) ? d.sparks : DEFAULT_SPARKS_CONFIG
            amplitudes = haskey(sparks_config, :amplitudes) ? collect(sparks_config.amplitudes) : nothing
            return new(r, p, pdot, r_pc, r_lc, alpha, magnetic_axis, rotation_axis, nsfield, fields, polar_caps, pc, open_lines, sparks, grid, potential, electric_field, drift_velocity, pot_minmax, sparks_locations, sparks_velocity, potential_simulation, spark_radius, spark_radii, line_of_sight, r_em, beta, los_lines, signal, pa, stokes_q, stokes_u, stokes_v, pulses, longitudes, ellipse_fit, p3, npulse, noise_level, output_num, sparks_config, amplitudes)
        end
    end
