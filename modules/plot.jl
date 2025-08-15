module Plot
    using GLMakie
    using GeometryBasics

    function pulsar(psr)
        f = Figure()
        ax = Axis3(f[1, 1])
        # Draw a sphere centered at (0,0,0) with radius r
        mesh!(ax, Sphere(Point3f(0, 0, 0), psr.r), color = :teal, shading = true)
        display(f)

    end

end # module end