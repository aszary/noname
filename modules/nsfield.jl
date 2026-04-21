module NSField

    """
        Field

    Container for the magnetic field configuration.

    Fields:
    - `rmax`: maximum radius for field line tracing [m]
    - `size`: number of steps per field line
    - `anomalies`: list of `Anomaly` objects
    """
    mutable struct Field
        rmax # maximum radius for which we will calculate magnetic field lines
        size # number of points in a line
        #nlos # number of 
        anomalies # list of magnetic anomalies

        function Field()
            rmax = 500_000 # 500 km
            size = 100
            nlos = 100
            anomalies = []
            new(rmax, size, anomalies)
        end

        function Field(json)
            rmax = json.field.rmax
            size = json.field.size
            anomalies = [Anomaly(a) for a in json.anomalies]
            new(rmax, size, anomalies)
        end
    end


    """
        Anomaly

    A displaced magnetic dipole anomaly. All radial coordinates are in stellar radii
    (r = 1 at the pulsar surface). Angles are in radians.

    Fields:
    - `r`: radial position of the anomaly [stellar radii]
    - `theta_r`: polar angle of the anomaly position [rad]
    - `phi_r`: azimuthal angle of the anomaly position [rad]
    - `m`: dipole moment strength relative to the global dipole [dimensionless]
    - `theta_m`: polar angle of the dipole moment orientation [rad]
    - `phi_m`: azimuthal angle of the dipole moment orientation [rad]

    JSON input: `r` in stellar radii, `theta`/`phi`/`theta_m`/`phi_m` in degrees.
    """
    mutable struct Anomaly
        r # spherical coordiante r,  location of the anommaly [in stellar radius]
        theta_r # spherical coordiante theta, location of the anommaly [rad]
        phi_r # spherical coordiante phi, location of the anommaly [rad]
        m # strength of the anomally with respect to the global dipol
        theta_m # spherical coordiante theta, orientation of the anommaly [rad]
        phi_m # spherical coordiante phi, orientation of the anommaly [rad]

        function Anomaly(json)
            new(json.r, deg2rad(json.theta), deg2rad(json.phi), json.m, deg2rad(json.theta_m), deg2rad(json.phi_m))
        end

    end



    # Projections of the anomaly position vector onto the local spherical basis
    # at field point (theta, phi). Used in the displaced-dipole formula.
    # r [stellar radii], theta, phi [rad]
    r1(a::Anomaly, theta, phi) = a.r * (sin(a.theta_r)*sin(theta)*cos(phi - a.phi_r) + cos(a.theta_r)*cos(theta))
    r2(a::Anomaly, theta, phi) = a.r * (sin(a.theta_r)*cos(theta)*cos(phi - a.phi_r) - cos(a.theta_r)*sin(theta))
    r3(a::Anomaly, phi)        = -a.r * sin(a.theta_r) * sin(phi - a.phi_r)

    # Projections of the anomaly dipole moment onto the local spherical basis
    # at field point (theta, phi). m [dimensionless, relative to global dipole].
    m1(a::Anomaly, theta, phi) = a.m * (sin(a.theta_m)*sin(theta)*cos(phi - a.phi_m) + cos(a.theta_m)*cos(theta))
    m2(a::Anomaly, theta, phi) = a.m * (sin(a.theta_m)*cos(theta)*cos(phi - a.phi_m) - cos(a.theta_m)*sin(theta))
    m3(a::Anomaly, phi)        = -a.m * sin(a.theta_m) * sin(phi - a.phi_m)

    # m·r scalar product (numerator term in displaced-dipole formula)
    # r [stellar radii], theta, phi [rad]
    function _T(a::Anomaly, r, theta, phi)
        m1(a,theta,phi)*r - (m1(a,theta,phi)*r1(a,theta,phi) + m2(a,theta,phi)*r2(a,theta,phi) + m3(a,phi)*r3(a,phi))
    end

    # Squared distance between the field point and the anomaly position [stellar radii²]
    _D(a::Anomaly, r, theta, phi) = a.r^2 + r^2 - 2*a.r*r*(sin(a.theta_r)*sin(theta)*cos(phi - a.phi_r) + cos(a.theta_r)*cos(theta))

    """
        Bb(a, r, theta, phi) -> (b_r, b_θ, b_φ)

    Magnetic field contribution of a single anomaly `a` at point (r, theta, phi)
    in spherical basis (e_r, e_θ, e_φ). Based on the displaced dipole formula.

    - `r`: radial distance [stellar radii], r = 1 at the pulsar surface
    - `theta`: polar angle [rad]
    - `phi`: azimuthal angle [rad]

    Returns field normalized to the global dipole strength.
    """
    function Bb(a::Anomaly, r, theta, phi)
        d  = _D(a, r, theta, phi)
        dl = d^(-2.5)
        t  = _T(a, r, theta, phi)
        b1 = -dl * (3*t*r1(a,theta,phi) - 3*t*r + d*m1(a,theta,phi))
        b2 = -dl * (3*t*r2(a,theta,phi)           + d*m2(a,theta,phi))
        b3 = -dl * (3*t*r3(a,phi)                 + d*m3(a,phi))
        return b1, b2, b3
    end

    # Global dipole field components (e_r, e_θ); H_φ = 0 by symmetry.
    # Normalized so that |H| = 1 at the magnetic pole on the surface (r=1, theta=0).
    # r [stellar radii], theta [rad]
    H1(r, theta) = 2*cos(theta) / r^3
    H2(r, theta) =   sin(theta) / r^3

    """
        B(f, r, theta, phi) -> (b_r, b_θ, b_φ)

    Sum of magnetic field contributions from all anomalies in `f.anomalies`
    at point (r, theta, phi). Does not include the global dipole.

    - `r`: radial distance [stellar radii], r = 1 at the pulsar surface
    - `theta`: polar angle [rad]
    - `phi`: azimuthal angle [rad]
    """
    function B(f::Field, r, theta, phi)
        b1 = b2 = b3 = 0.0
        for a in f.anomalies
            db1, db2, db3 = Bb(a, r, theta, phi)
            b1 += db1; b2 += db2; b3 += db3
        end
        return b1, b2, b3
    end

    """
        BFull(f, r, theta, phi) -> scalar

    Magnitude of the total magnetic field (anomalies + global dipole) at (r, theta, phi).

    - `r`: radial distance [stellar radii], r = 1 at the pulsar surface
    - `theta`: polar angle [rad]
    - `phi`: azimuthal angle [rad]
    """
    function BFull(f::Field, r, theta, phi)
        b1, b2, b3 = B(f, r, theta, phi)
        sqrt((b1 + H1(r,theta))^2 + (b2 + H2(r,theta))^2 + b3^2)
    end

    """
        BSph(f, r, theta, phi) -> (b_r, b_θ, b_φ)

    Total magnetic field vector (anomalies + global dipole) in spherical coordinates.

    - `r`: radial distance [stellar radii], r = 1 at the pulsar surface
    - `theta`: polar angle [rad]
    - `phi`: azimuthal angle [rad]
    """
    function BSph(f::Field, r, theta, phi)
        br, bθ, bφ = B(f, r, theta, phi)
        br += H1(r, theta)
        bθ += H2(r, theta)
        return br, bθ, bφ
    end

    """
        BVec(f, r, theta, phi) -> (Bx, By, Bz)

    Total magnetic field vector (anomalies + global dipole) in Cartesian coordinates.

    - `r`: radial distance [stellar radii], r = 1 at the pulsar surface
    - `theta`: polar angle [rad]
    - `phi`: azimuthal angle [rad]
    """
    function BVec(f::Field, r, theta, phi)
        br, bθ, bφ = B(f, r, theta, phi)
        br += H1(r, theta)
        bθ += H2(r, theta)
        Bx = sin(theta)*cos(phi)*br + cos(theta)*cos(phi)*bθ - sin(phi)*bφ
        By = sin(theta)*sin(phi)*br + cos(theta)*sin(phi)*bθ + cos(phi)*bφ
        Bz = cos(theta)*br          - sin(theta)*bθ
        return Bx, By, Bz
    end

end # module end
