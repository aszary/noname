module Signal
    include("functions.jl")
    include("transformations.jl")
    include("field.jl")
    include("lines.jl")
    include("sparks.jl")

    using LinearAlgebra

    function calculate_spark_distances(psr; height=500000.0, points=1000)
        # 1. Get Sparks
        spark_coords = Vector{Vector{Float64}}()
        if psr.sparks !== nothing
            gr = psr.grid
            for s in psr.sparks
                if length(s) == 2; push!(spark_coords, [gr[1][Int(s[1])], gr[2][Int(s[2])], gr[3][Int(s[1]), Int(s[2])]]);
                elseif length(s) == 3; push!(spark_coords, [s[1], s[2], s[3]]); end
            end
        end
        if isempty(spark_coords); return []; end

        # 2. Get Geometry from shared function
        # To wywołanie wykonuje całą ciężką matematykę z Functions.jl
        traj_data = Lines.calculate_los(psr; height=height, points=points)

        # 3. Calculate Distances
        spark_dist = [] 
        
        for (phi, p_surf, p_end) in traj_data
            dists = Float64[]
            for sp in spark_coords
                push!(dists, norm(p_surf - sp))
            end
            push!(spark_dist, (phi, dists))
        end
        println(spark_dist)
        return spark_dist
    end
end # module end