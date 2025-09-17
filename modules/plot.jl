module Plot
    using GLMakie
    using GeometryBasics
    using LinearAlgebra

    include("functions.jl")
    include("field.jl")


    function pulsar(psr)
        fi = psr.fields


        fig, ax, p = mesh(Sphere(Point3f(0, 0, 0), psr.r), color = (:teal, 0.7), transparency = true)
        
        #f = Figure()
        
        #ax = Axis3(f[1, 1], aspect = :equal, perspectiveness = 0)
        # Draw a sphere centered at (0,0,0) with radius r
        #mesh!(ax, Sphere(Point3f(0, 0, 0), psr.r), color = (:teal, 0.7), transparency = true)
        # rotation axis
        rot_vec = Functions.spherical2cartesian(psr.rotation_axis)
        # magnetic axis
        mag_vec = Functions.spherical2cartesian(psr.magnetic_axis)

        #println(psr.rotation_axis)
        #println(psr.magnetic_axis)

        #println(rot_vec)
        #println(mag_vec)

        arrows3d!(ax,[0,], [0,0], [0,0], [rot_vec[1], mag_vec[1]], [rot_vec[2], mag_vec[2]], [rot_vec[3], mag_vec[3]], color = [:red, :blue])#,  shaftradius = 0.01, tipradius = 0.01, tiplength=0.01)
        #arrows3d!(ax, [10000], [10], [10], [1, 0.5], color = [:red, :blue]) #,  shaftradius = 0.01, tipradius = 0.01, tiplength=0.01)

        r = 1.2e4
        theta = 30 # w stopniach
        phi = 0
        p_sp = [r, deg2rad(theta), phi]
        p_car = Functions.spherical2cartesian(p_sp)

        #draw points
        scatter!(ax, [0, p_car[1]], [0, p_car[2]], [1.2e4, p_car[3]], markersize=10, color=:red)

        mx = 2e4
        #limits!(ax, -mx, mx, -mx, mx, -mx, mx)

        
        plot_magnetic_field!(ax, psr.fields)
        #plot_magnetic_lines!(ax, psr.fields)
        plot_polar_cap!(ax, psr; color=:red, linewidth=2.0)
        plot_magnetic_lines_from_polar_cap!(ax, psr.fields)

        
        # Draw light cylinder
        #light_cylinder(psr, ax)
        display(fig)
    end

    function plot_magnetic_field!(ax, fi; marker_size=10, colormap=:plasma)
        # Convert positions to flat Float32 arrays
        xs = Float32[p[1] for p in fi.locations]
        ys = Float32[p[2] for p in fi.locations]
        zs = Float32[p[3] for p in fi.locations]

        # Magnitude of B for coloring
        Bmag = Float32[norm(v) for v in fi.magnetic_fields]

        # Plot as scatter (dots)
        scatter!(ax, xs, ys, zs;
                markersize=marker_size,
                color=Bmag,
                colormap=colormap,
                transparency=true,
                alpha=0.5)
    end



    
    function plot_magnetic_lines!(ax, fv)
        # fv = psr.fields
        for ml in fv.magnetic_lines
            lines!(ax, ml[1], ml[2], ml[3], color=:blue, linewidth=1.0, transparency=true)
        end
    end

    function plot_polar_cap!(ax, psr; color=:red, linewidth=2.0)
        # group points by cap (north/south)
        caps = Dict(:north => (Float64[], Float64[], Float64[]),
                    :south => (Float64[], Float64[], Float64[]))

        for (r_pc, θ, φ) in psr.polar_cap
            r = psr.r
            x = r * sin(θ) * cos(φ)
            y = r * sin(θ) * sin(φ)
            z = r * cos(θ)

            if z ≥ 0
                push!(caps[:north][1], x)
                push!(caps[:north][2], y)
                push!(caps[:north][3], z)
            else
                push!(caps[:south][1], x)
                push!(caps[:south][2], y)
                push!(caps[:south][3], z)
            end
        end

        # plot separately
        lines!(ax, caps[:north]..., color=color, linewidth=linewidth)
        lines!(ax, caps[:south]..., color=color, linewidth=linewidth)
    end


   function plot_magnetic_lines_from_polar_cap!(ax, f; color=:blue, linewidth=1.0)
        for ml in f.magnetic_lines_from_cap
            lines!(ax, ml[1], ml[2], ml[3], color=color, linewidth=linewidth)
        end
    end

     
    function light_cylinder(psr, ax)
        # Create wireframe cylinder with lines
        θ = range(0, 2π, length=200)
        height = 3 * psr.r
        radius = 2* psr.r
        #radius = 2* psr.r_lc
        
        # Vertical lines
        for t in θ[1:end-1]  # Skip last point to avoid duplication
            x_line = [radius * cos(t), radius * cos(t)]
            y_line = [radius * sin(t), radius * sin(t)]
            z_line = [-2 * psr.r, height]
            lines!(ax, Point3f.(x_line, y_line, z_line), color = :red, alpha = 0.5, linewidth = 2)
        end
        
        # Horizontal circles at different heights
        z_levels = range(0, height, length=5)
        for z in z_levels
            x_circle = [radius * cos(t) for t in θ]
            y_circle = [radius * sin(t) for t in θ]
            z_circle = fill(z, length(θ))
            lines!(ax, Point3f.(x_circle, y_circle, z_circle), color = :red, alpha = 0.5, linewidth = 2)
        end
        # Set equal limits for all axes
        max_extent = max(height, radius)  # Assuming cylinder height is 2*psr.r
        limits!(ax, -max_extent, max_extent, -max_extent, max_extent, -max_extent, max_extent)

    end

end # module end