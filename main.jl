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


    # Pulsar type and its JSON constructor
    include("modules/pulsar.jl")


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
        Sparks.simulate_sparks_mc(psr;n_steps=5000)
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
        Sparks.simulate_sparks_lbc(psr; n_steps=500, co_angl=-90.0)
        Sparks.save_sparks(psr; num=psr.output_num)

        #Plot.sparks(psr)
        Sparks.load_sparks(psr; num=psr.output_num)

        #Signal.generate_signal(psr; noise_level=psr.noise_level) # old same sizes!
        Signal.generate_signal_radii(psr; noise_level=psr.noise_level) # new
        Signal.generate_pulses(psr)
        
        #Plot.signal(psr)
        Plot.pulses(psr)
        #Plot.pulses0(psr)
        #Plot.pulses1(psr)
        
    end


    function generate_signal()
        #psr = Pulsar("input/1.json")
        #psr = Pulsar("input/2.json")
        #psr = Pulsar("input/3.json")
        #psr = Pulsar("input/4.json")
        #psr = Pulsar("input/15.json")
        psr = Pulsar("input/B1839-04.json")
        #psr = Pulsar("input/B1839-04_fit.json")
        

        Lines.init_line_of_sight(psr, num=psr.nsfield.nlos)
        Lines.calculate_line_of_sight(psr)

        Lines.generate_open!(psr, num=psr.nsfield.nopen)

        sc = psr.sparks_config
        si = sc.init
        if si.method == "ellipse"
            Sparks.init_sparks1_ellipse!(psr; rfs=collect(si.rfs), num=si.num, spacing=get(si, :spacing, "t"), phase=get(si, :phase, 0.0), num_mode=get(si, :num_mode, "arc"))
        elseif si.method == "dipolar"
            Sparks.init_sparks1!(psr; rfs=collect(si.rfs), num=si.num)
        elseif si.method == "dipolar2"
            Sparks.init_sparks2!(psr; rfs=collect(si.rfs), num=si.num)
        elseif si.method == "none"
            # skip spark initialization
        else
            error("Unknown spark init method: $(si.method). Use \"ellipse\", \"dipolar\", \"dipolar2\", or \"none\".")
        end
        if sc.model == "mc"
            Sparks.simulate_sparks_mc(psr; n_steps=sc.mc.n_steps, save_every=sc.mc.save_every, speedup=sc.mc.speedup)
        elseif sc.model == "solidbody"
            Sparks.simulate_sparks_solidbody(psr)
        elseif sc.model == "lbc"
            Sparks.simulate_sparks_lbc(psr; n_steps=psr.npulse, co_angl=sc.lbc.co_angl)
        else
            error("Unknown spark model: $(sc.model). Use \"mc\", \"solidbody\", or \"lbc\".")
        end

        #Sparks.save_sparks(psr; num=psr.output_num)
        #Plot.sparks(psr)
        #Sparks.load_sparks(psr; num=psr.output_num)


        #Signal.generate_signal(psr; noise_level=psr.noise_level) # old  obsolete same sizes! NO PA
        #Signal.generate_signal_radii(psr; noise_level=psr.noise_level, v_scale=0.3) # new with full Stokes and radii
        #Signal.generate_signal_ola(psr; noise_level=psr.noise_level, v_scale=0.3) # new with full Stokes and elliptical sparks
        Signal.generate_signal_new(psr; noise_level=psr.noise_level, v_scale=0.3) # new with full Stokes and elliptical sparks based on ellipse fit
        Signal.generate_pulses(psr)


        Plot.signal(psr)
        Plot.pulses(psr, number=psr.npulse)
        #Plot.pulses0(psr)
        #Plot.pulses1(psr)
        #Plot.average_stokes(psr)
        #Plot.polarization_vector_study(psr)
        Plot.lrfs(psr, darkness=0.3)
        #Plot.two_dfs(psr, darkness=0.3)
    end


    function model_field()
        #psr = Pulsar("input/1.json")
        #psr = Pulsar("input/2.json")
        #psr = Pulsar("input/3.json")
        psr = Pulsar("input/B1839-04.json")

        Lines.init_line_of_sight(psr, num=psr.nsfield.nlos)
        Lines.calculate_line_of_sight(psr)

        Lines.generate_open!(psr, num=psr.nsfield.nopen)

        Lines.generate_closed!(psr)

        #println(psr.nsfield)

        #Plot.closed_lines(psr)
        #Plot.anomalies(psr)
        Plot.anomalies2D(psr)
        #Plot.polar_cap2D(psr)

       
    end



    function main()

        #full_grid()
        #small_grids()
        #full_plus_smallgrids()

        #generate_signal_dipole()
        generate_signal()

        #model_field()

        println("Bye")
    end


end # module end

NoName.main()