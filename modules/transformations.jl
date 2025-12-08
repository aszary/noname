module Transformations    

    include("functions.jl")

    """
        rotation_matrix_yx(α, β)

    Create rotation matrix for Euler angles: first rotate by α around X-axis, then by β around Y-axis.
    """
    function rotation_matrix_yx(α, β)
        return [
            cos(β)           sin(α)*sin(β)   sin(β)*cos(α)
            0.0              cos(α)          -sin(α)
            -sin(β)          cos(β)*sin(α)   cos(β)*cos(α)
        ]
    end

    """
        get_euler_angles(axis)

    Calculate Euler angles (α, β) needed to align z-axis with the given vector.
    """
    function get_euler_angles(axis)
        x, y, z = axis
        
        α = atan(y, z)
        
        if z != 0.0
            β = atan(-x / z * cos(α))
        elseif y != 0.0
            β = atan(-x / y * sin(α))
        else
            β = π / 2
        end
        
        return α, β
    end

    """
        rotate_to_vector(axis, vec)

    Transform vector `vec` to coordinate system where z-axis is aligned with `axis`.

    # Arguments
    - `axis`: vector in cartessian coordinates
    - `vec`: vector to be rotated (3-element array)

    # Returns
    - Rotated vector in the new coordinate system
    """
    function rotate_to_vector(axis, vec)
        α, β = get_euler_angles(axis)
        R = rotation_matrix_yx(α, β)
        return R * vec
    end

    """
        rotate_back(axis, vec_rot)

    Inverse transformation: rotate vector back to original coordinate system.

    # Arguments
    - `axis`: vector in cartesian coordinates
    - `vec_rot`: vector in rotated coordinate system

    # Returns
    - Vector in original coordinate system
    """
    function rotate_back(axis, vec_rot)
        α, β = get_euler_angles(axis)
        R = rotation_matrix_yx(α, β)
        return R' * vec_rot  # transpose = inverse for rotation matrix
    end


    """
        beaming(axis, θ, φ)

    Create a point on a cone around `axis` with opening angle θ and azimuth φ.

    This is the key method for pulsar beam geometry: given the magnetic axis direction,
    it generates points on the emission cone.

    # Arguments
    - `axis`: magnetic axis direction in cartesian coordinates
    - `θ`: opening angle of the cone (radians)
    - `φ`: azimuthal angle around the cone (radians)

    # Returns
    - Spherical coordinates (r, θ, φ) in the original frame
    """
    function beaming(axis, θ, φ)
        # Get spherical coords of vector in rotated frame (should be along z-axis)
        vec_rot = rotate_to_vector(axis, axis)
        sph = Functions.cartesian2spherical(vec_rot)
        
        # Add angular offsets to create point on cone
        sph[2] += θ  # polar angle offset
        sph[3] += φ  # azimuthal offset
        
        # Convert back to Cartesian, rotate back, then to spherical
        cart = Functions.spherical2cartesian(sph)
        cart_original = rotate_back(axis, cart)
        
        return Functions.cartesian2spherical(cart_original)
    end

    """
        beaming(t::Transformations, θ, φ::AbstractVector)

    Vectorized version: generate multiple points on the cone for an array of azimuthal angles.

    # Returns
    - Matrix where each column is (r, θ, φ) for one point
    """
    function beaming(vector, θ, φ::AbstractVector)
        return hcat([beaming(vector, θ, φ_i) for φ_i in φ]...)
    end

end # module