module NSField

    mutable struct Field
        rmax # maximum radius for which we will calculate magnetic field lines
        size # number of points in a line
        anomalies # list of magnetic anomalies

        function Field()
            rmax = 50_000 # 50 km 
            size = 100
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


    mutable struct Anomaly
        r # spherical coordiante r,  location of the anommaly [in stellar radius]
        theta_r # spherical coordiante theta, location of the anommaly [deg]
        phi_r # spherical coordiante phi, location of the anommaly [deg]
        m # strength of the anomally with respect to the global dipol
        theta_m # spherical coordiante theta, orientation of the anommaly [deg]
        phi_m # spherical coordiante phi, orientation of the anommaly [deg]
        
        function Anomaly(json)
            new(json.r, json.theta, json.phi, json.m, json.theta_m, json.phi_m)
        end

    end






end # module end