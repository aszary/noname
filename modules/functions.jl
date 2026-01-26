module Functions 
    using PhysicalConstants.CODATA2018 # in SI units
    using LinearAlgebra


    mutable struct LineOfSight
        phase::Vector{Float64}
        phase_deg::Vector{Float64}
        sph::Vector{Vector{Float64}}
        cart::Vector{Vector{Float64}}
        open_flags::Vector{Bool}

        # empty constructor
        function LineOfSight()
            new([], [], [], [], [])
        end

        # full constructor
        function LineOfSight(phase, phase_deg, sph, cart, open_flags)
            new(phase, phase_deg, sph, cart, open_flags)
        end
    end

    """
    Calculates rlc in meters
    """
    function rlc(p)
        c = SpeedOfLightInVacuum.val # no units hereafter
        return c * p / (2 * pi)
    end

    """
    spherical2cartesian(spherical)

    Converts spherical cordinates to cartesian ones...
    """
    function spherical2cartesian(spherical)
        x = spherical[1] * sin(spherical[2]) * cos(spherical[3])
        y = spherical[1] * sin(spherical[2]) * sin(spherical[3])
        z = spherical[1] * cos(spherical[2])
        return [x, y, z]
    end

    
    """

    Converts vector in spehrical coordinates (vec_sph) at position in pos_sph to cartesian coordiantes
    """
    function vec_spherical2cartesian(pos_sph, vec_sph)
        r = pos_sph[1]
        theta = pos_sph[2]
        phi = pos_sph[3]

        v_r = vec_sph[1]
        v_theta = vec_sph[2]
        v_phi = vec_sph[3]

        v_x = v_r * sin(theta) * cos(phi) + v_theta * cos(theta) * cos(phi) - v_phi * sin(phi)
        v_y = v_r * sin(theta) * sin(phi) + v_theta * cos(theta) * sin(phi) + v_phi * cos(phi)
        v_z = v_r * cos(theta) - v_theta * sin(theta)
        return [v_x, v_y, v_z]
    end
    

    """
    Converts cartesian cordinates to spherical ones...
    """
    function cartesian2spherical(cartesian)
        x = cartesian[1]
        y = cartesian[2]
        z = cartesian[3]
        r = sqrt(x^2 + y^2 + z^2)
        if z > 0
            theta = atan(sqrt(x^2 + y^2) / z)
        elseif z < 0
            theta = pi + atan(sqrt(x^2 + y^2) / z)
        elseif (z == 0) && (x*y != 0)
            theta = pi
        else
            theta = 0 # undefined changed to zero
        end
        if x > 0
            phi = atan(y / x)
        elseif (x < 0) && (y >=0)
            phi = atan(y / x) + pi
        elseif (x < 0) && (y < 0)
            phi = atan(y / x) - pi
        elseif (x == 0) && (y > 0)
            phi = pi
        elseif (x == 0) && (y < 0)
            phi = -pi
        elseif (x == 0) && (y ==0) # undefined changed to zero
            phi = 0
        end
        return [r, theta, phi]
    end


    """
    returns radius of the polar in meters
    """
    function rdp(p, r)
        #https://juliaphysics.github.io/PhysicalConstants.jl/stable/reference/
        c = SpeedOfLightInVacuum.val # no units hereafter
        #println((2. * pi * r ^ 3. / (SpeedOfLightInVacuum * p)) ^ 0.5)
        return (2. * pi * r ^ 3. / (c * p)) ^ 0.5
    end

       """
    theta_max(z, psr)

    Eq. 3.23 in the handbook

    Note this is not an opening angle! (Note 1.5 differecence E.g 3.29 in the handbook)

    # Arguments

    - z: distance from the star's center [in stellar radius]
    - psr: pulsar class (struct) - to get period and star r

    returns the theta component of the last open magnetic field line for a given distance from the star's
    center and pulsar period [in radians]
    """
    function theta_max(z, psr)
        return asin(sqrt(z * psr.r / rlc(psr.p)))
    end

       """
    numpy_gradient_2d(A::AbstractMatrix) -> (grad_y, grad_x)

    Computes the gradient of a 2D array using finite differences, consistent with `numpy.gradient` implementation.

    # Arguments
    - `A::AbstractMatrix`: Input 2D array (matrix)

    # Returns
    - A tuple `(grad_y, grad_x)` containing:
    - `grad_y`: Gradient along the first dimension (rows, Y-axis)
    - `grad_x`: Gradient along the second dimension (columns, X-axis)

    # Algorithm Description
    The function computes gradients using:
    - **One-sided differences** at array boundaries:
    - First row/column: `grad[1] = A[2] - A[1]`
    - Last row/column: `grad[end] = A[end] - A[end-1]`
    - **Central differences** for interior points:
    - `grad[i] = (A[i+1] - A[i-1]) / 2`

    This behavior is identical to `numpy.gradient` for 2D arrays.    
    
    """
    function numpy_gradient_2d(A)
        grad_y = similar(A, Float64)
        grad_x = similar(A, Float64)
        
        ny, nx = size(A)
        
        # Gradient w kierunku Y (dim=1, wiersze)
        grad_y[1, :] = A[2, :] - A[1, :]
        grad_y[end, :] = A[end, :] - A[end-1, :]
        for i in 2:ny-1
            grad_y[i, :] = (A[i+1, :] - A[i-1, :]) / 2
        end
        
        # Gradient w kierunku X (dim=2, kolumny)
        grad_x[:, 1] = A[:, 2] - A[:, 1]
        grad_x[:, end] = A[:, end] - A[:, end-1]
        for j in 2:nx-1
            grad_x[:, j] = (A[:, j+1] - A[:, j-1]) / 2
        end
        
        return (grad_y, grad_x)
    end

    # Rodrigues rotation matrix to rotate vector a -> b (both 3-element)
    # If a and b are parallel, returns I.
    function rotation_matrix_from_to(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
        a_u = a / norm(a)
        b_u = b / norm(b)
        # handle near-parallel / anti-parallel
        c = dot(a_u, b_u)
        if isapprox(c, 1.0; atol=1e-12)
            return Matrix{Float64}(I, 3, 3)
        elseif isapprox(c, -1.0; atol=1e-12)
            # 180° rotation: pick an arbitrary orthogonal axis
            # find vector orthogonal to a_u
            v = abs(a_u[1]) < 0.9 ? [1.0,0.0,0.0] : [0.0,1.0,0.0]
            k = cross(a_u, v)
            k = k / norm(k)
            K = [  0.0  -k[3]  k[2];
                k[3]   0.0  -k[1];
                -k[2]  k[1]   0.0 ]
            return I - 2.0 * K * K
        else
            k = cross(a_u, b_u)
            K = [  0.0   -k[3]   k[2];
                k[3]   0.0   -k[1];
                -k[2]  k[1]    0.0 ]
            s = norm(k)
            return I + K + K*K * ((1 - c)/(s^2))
        end
    end

    # small helper: spherical (r, theta, phi) -> cartesian (x,y,z)
    # theta = polar angle from +z, phi = azimuth from +x
    sph_to_cartesian(r, θ, φ) = (r * sin(θ) * cos(φ), r * sin(θ) * sin(φ), r * cos(θ))

    # cartesian -> spherical (r, θ, φ)
    function cartesian_to_spherical(x::Real, y::Real, z::Real)
        r = sqrt(x^2 + y^2 + z^2)
        if r == 0.0
            return (0.0, 0.0, 0.0)
        end
        θ = acos(clamp(z / r, -1.0, 1.0))
        φ = atan(y, x)
        return (r, θ, φ)
    end

    # Transformations struct: precompute rotation that maps +z -> rotation_axis
    struct Transformations
        R::Matrix{Float64}   # rotation matrix: maps vector in local frame -> rotated (global) frame
    end

    function Transformations(rotation_axis::AbstractVector{<:Real})
        # We need rotation that maps z-axis (0,0,1) to rotation_axis
        z = [0.0, 0.0, 1.0]
        R = rotation_matrix_from_to(z, rotation_axis)
        return Transformations(R)
    end

    """
        beaming(t::Transformations; θ, φ)

    Return the direction (r, θ_rot, φ_rot) obtained by taking the
    unit vector with spherical angles (θ, φ) in the local frame (z-axis),
    rotating it by t.R, and converting to spherical coordinates in the rotated/global frame.
    """
    function beaming(t::Transformations; θ::Real, φ::Real)
        x, y, z = sph_to_cartesian(1.0, θ, φ)   # unit vector in local frame
        v = t.R * [x, y, z]                    # rotated vector in global frame
        return cartesian_to_spherical(v[1], v[2], v[3])  # (r, θ_rot, φ_rot)
    end


end # module end

