module Plot
    using GLMakie
    using GeometryBasics
    include("functions.jl")



    function pulsar(psr)
        #fig = Figure()
        #ax = Axis3(f[1, 1], aspect = :equal)

        fig, ax, p = mesh(Sphere(Point3f(0, 0, 0), psr.r), color = (:teal, 0.7), transparency = true) # better camera control (Scene), but zlims does not work

        # Draw a sphere centered at (0,0,0) with radius r
        #mesh!(ax, Sphere(Point3f(0, 0, 0), psr.r), color = (:teal, 0.7), transparency = true)
        
        # rotation axis
        rot_vec = Functions.spherical2cartesian(psr.rotation_axis)
        # magnetic axis
        mag_vec = Functions.spherical2cartesian(psr.magnetic_axis)

        arrows3d!(ax,[0,], [0,0], [0,0], [rot_vec[1], mag_vec[1]], [rot_vec[2], mag_vec[2]], [rot_vec[3], mag_vec[3]], color = [:red, :blue])#,  shaftradius = 0.01, tipradius = 0.01, tiplength=0.01)

        #draw test points
        r = 1.2e4
        theta = 30 # w stopniach
        phi = 0
        p_sp = [r, deg2rad(theta), phi]
        p_car = Functions.spherical2cartesian(p_sp)
        scatter!(ax, [0, p_car[1]], [0, p_car[2]], [1.2e4, p_car[3]], markersize=10, color=:red)

        magnetic_field!(ax, psr)

        magnetic_lines!(ax, psr)

        # draw polar caps
        for pc in psr.polar_caps
            lines!(ax, pc[1], pc[2], pc[3], color=:red, linewidth=1)
        end

        # draw open lines
        for i in eachindex(psr.open_lines)
            for j in eachindex(psr.open_lines[i])
                lines!(ax, psr.open_lines[i][j][1], psr.open_lines[i][j][2], psr.open_lines[i][j][3], color=:green)         
            end
        end


        mx = 2e4
        #limits!(ax, -mx, mx, -mx, mx, -mx, mx)

        # Draw light cylinder
        #light_cylinder(psr, ax)

        
        cam3d!(ax.scene, eyeposition=[5000, 5000, 5000])

        #= 
       cam = camera(ax.scene) 
        # Try accessing the scene's camera directly
        println("Scene camera type: ", typeof(cam))
        println("Scene camera fields: ", fieldnames(typeof(cam)))
        println("Scene camera properties: ", propertynames(cam))

        # Add a button to print camera state
        #button = Button(f[7, 1], label = "Print Camera State")
        button = Button(fig[1, 1], label = "Print", 
               width = 80, height = 25,
               halign = :right, valign = :top,
               tellwidth = false, tellheight = false)


        on(button.clicks) do n
            println("\n--- Camera State (Click $n) ---")
            println("Camera type: ", typeof(cam))
            println("Eye position: ", cam.eyeposition[])
            println("View direction: ", cam.view_direction[])
            println("Up vector: ", cam.upvector[])
        end
        
        cam.eyeposition[]= [21818.455, 23425.26, 25757.416]
        notify(cam.eyeposition)


        # Update the camera
        update_cam!(ax.scene)
        =#
        display(fig)



    end


    """
        Plots light cylinder
    """
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


    function magnetic_field!(ax, psr)
        fv = psr.fields

        positions = fv.locations
        magnitudes = [Functions.norm(B) for B in fv.magnetic_fields]

        xs = [p[1] for p in positions]
        ys = [p[2] for p in positions]
        zs = [p[3] for p in positions]

        scatter!(ax, xs, ys, zs, color=magnitudes, colormap=:viridis, markersize=5)
    end 

    function magnetic_lines!(ax, psr)
        fv = psr.fields
        for line in fv.magnetic_lines
            xs, ys, zs = line[1], line[2], line[3]
            lines!(ax, xs, ys, zs, color=:blue, linewidth=1)
        end
    end


end # module end