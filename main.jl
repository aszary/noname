module NoName

    include("modules/functions.jl")
    include("modules/plot.jl")
    

    mutable struct Pulsar
        r # pulsar radiuis [in meters]
        p # pulsar period [in s]
        pdot # pulsar period derivative in [s/s]
        r_lc # light cylinder radius [in meters]


        function Pulsar()
            r = 10000
            p = 1
            pdot = 1e-15
            r_lc = Functions.rlc(p, pdot) 
            return new(r, p, pdot, r_lc)
        end
    end


    function main()

        psr = Pulsar()
        println(fieldnames(Pulsar))

        Plot.pulsar(psr)
        println("Bye")
    end


end # module end

NoName.main()