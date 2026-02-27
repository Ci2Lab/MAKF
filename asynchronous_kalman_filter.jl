using LinearAlgebra

mutable struct KalmanFilter
    x::Vector{Float64}          # state vector
    P::Matrix{Float64}          # state covariance matrix
    F::Matrix{Float64}          # state transition matrix
    Q::Matrix{Float64}          # process covariance matrix
    K::Matrix{Float64}          # Kalman filter gain
    H::Matrix{Float64}          # Measurement matrix
    R::Matrix{Float64}          # Measurement covariance matrix

    function KalmanFilter(x, P, F, Q, H, R, K)
        # default constructor, initialize variables
        new(x, P, F, Q, H, R, K)
    end
end

function predict(kf::KalmanFilter)
    kf.x = kf.F * kf.x
    kf.P = kf.F * kf.P * kf.F' + kf.Q
end

function update(kf::KalmanFilter, z::Vector{Float64})
    z_pred = kf.H * kf.x
    y = z - z_pred
    Ht = kf.H'
    S = kf.H * kf.P * Ht + kf.R
    kf.K = kf.P * Ht * inv(S)
    kf.x = kf.x + kf.K * y
    kf.P = (I - kf.K * kf.H) * kf.P
end


mutable struct Tracking

    kf::KalmanFilter
    is_initialized::Bool
    previous_timestamp::Int64

    # acceleration noise components
    noise_ax::Float64
    noise_ay::Float64

    function Tracking(x0, P0, F, Q, H, R, K)
        is_initialized = false
        previous_timestamp = 0

        # create a 6D state vector, we don't know yet the values of the x state
        kf = KalmanFilter(
            x0, P0, F, Q, H, R, K
        )

        # set the acceleration noise components
        noise_ax = 10.0
        noise_ay = 10.0

        new(kf, is_initialized, previous_timestamp, noise_ax, noise_ay)
    end
end


function process_measurement_acc(tracking::Tracking, measurement::Vector)

    if measurement[2] == "Total Station"
        # measurement is from total station
        tracking.kf.H = [1.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0] # Measurement matrix total station
        tracking.kf.R = diagm([10.0, 10.0, 10.0]) # Measurement covariance total station

        if !tracking.is_initialized && measurement[2] == "Total Station"
            println("Kalman Filter Initialization")
            # set the state with the initial location and zero velocity
            tracking.kf.x[1] = measurement[4][1]
            tracking.kf.x[2] = measurement[4][2]
            tracking.kf.x[3] = measurement[4][3]

            tracking.previous_timestamp = datetime2unix(measurement[1])
            tracking.is_initialized = true
            return
        end

    elseif measurement[2] == "SAT-INSAR"

        return

        # measurement is from satellite
        tracking.kf.H = [1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0] # Measurement matrix satellite
        tracking.kf.R = diagm([12.5^2, 10.0^2]) # Measurement covariance satellite

        if !tracking.is_initialized && measurement[2] == "SAT-INSAR"
            println("Kalman Filter Initialization")
            # set the state with the initial location and zero velocity
            tracking.kf.x[1] = measurement[4][1]
            tracking.kf.x[3] = measurement[4][2]

            tracking.previous_timestamp = datetime2unix(measurement[1])
            tracking.is_initialized = true
            return
        end

    elseif measurement[2] == "GB-INSAR"

        return

        # measurement is from ground based
        tracking.kf.H = [0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0] # Measurement matrix ground based
        tracking.kf.R = diagm([9.9869^2, 0.3061^2]) # Measurement covariance ground based

        if !tracking.is_initialized && measurement[2] == "GB-INSAR"
            println("Kalman Filter Initialization")
            # set the state with the initial location and zero velocity
            tracking.kf.x[2] = measurement[4][1]
            tracking.kf.x[3] = measurement[4][2]

            tracking.previous_timestamp = datetime2unix(measurement[1])
            tracking.is_initialized = true
            return
        end

    else
        println("Measurement type not recognized")
        return

    end

    # compute the time elapsed between the current and previous measurements
    dt = (datetime2unix(measurement[1]) - tracking.previous_timestamp)  # dt - expressed in seconds
    tracking.previous_timestamp = datetime2unix(measurement[1])

    dt_2 = dt * dt
    dt_3 = dt_2 * dt
    dt_4 = dt_3 * dt

    # Modify the F matrix so that the time is integrated
    tracking.kf.F[1, 4] = dt
    tracking.kf.F[2, 5] = dt
    tracking.kf.F[3, 6] = dt

    # set the process covariance matrix Q
    tracking.kf.Q = [(dt_4/4)*tracking.noise_ax 0.0 0.0 (dt_3/2)*tracking.noise_ax 0.0 0.0
        0.0 (dt_4/4)*tracking.noise_ay 0.0 0.0 dt_3/2*tracking.noise_ay 0.0
        0.0 0.0 (dt_4/4)*tracking.noise_ay 0.0 0.0 (dt_3/2)*tracking.noise_ay
        (dt_3/2)*tracking.noise_ax 0.0 0.0 dt_2*tracking.noise_ax 0.0 0.0
        0.0 (dt_3/2)*tracking.noise_ay 0.0 0.0 dt_2*tracking.noise_ay 0.0
        0.0 0.0 (dt_3/2)*tracking.noise_ay 0.0 0.0 dt_2*tracking.noise_ay]

    # predict
    predict(tracking.kf)

    # measurement update
    update(tracking.kf, vec(measurement[4]))

end


function process_measurement(tracking::Tracking, measurement::Vector, sensors::Vector)

    if measurement[2] == "Total Station" && "Total Station" in sensors
        # measurement is from total station
        tracking.kf.H = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0] # Measurement matrix total station
        tracking.kf.R = diagm([4.7384^2, 7.9214^2, 11.02575^2]) # Measurement covariance total station

        # if NaN in measurement, call only prediction step
        if isnan(measurement[4][1])
            predict(tracking.kf)
            return measurement[1]
        end

        if !tracking.is_initialized
            println("Kalman Filter Initialization")
            # set the state with the initial location and zero velocity
            tracking.kf.x[1] = measurement[4][1]
            tracking.kf.x[2] = measurement[4][2]
            tracking.kf.x[3] = measurement[4][3]

            tracking.previous_timestamp = datetime2unix(measurement[1])
            tracking.is_initialized = true
            return
        end

    elseif measurement[2] == "SAT-INSAR" && "SAT-INSAR" in sensors

        # measurement is from satellite
        tracking.kf.H = [1.0 0.0 0.0; 0.0 0.0 1.0] # Measurement matrix satellite
        tracking.kf.R = diagm([12.5^2, 10.0^2]) # Measurement covariance satellite

        if !tracking.is_initialized
            println("Kalman Filter Initialization")
            # set the state with the initial location and zero velocity
            tracking.kf.x[1] = measurement[4][1]
            tracking.kf.x[3] = measurement[4][2]

            tracking.previous_timestamp = datetime2unix(measurement[1])
            tracking.is_initialized = true
            return
        end

    elseif measurement[2] == "GB-INSAR" && "GB-INSAR" in sensors

        # measurement is from ground based
        tracking.kf.H = [0.0 1.0 0.0; 0.0 0.0 1.0] # Measurement matrix ground based
        tracking.kf.R = diagm([9.9869^2, 0.3061^2]) # Measurement covariance ground based

        if !tracking.is_initialized
            println("Kalman Filter Initialization")
            # set the state with the initial values
            tracking.kf.x[2] = measurement[4][1]
            tracking.kf.x[3] = measurement[4][2]

            tracking.previous_timestamp = datetime2unix(measurement[1])
            tracking.is_initialized = true
            return
        end

    elseif measurement[2] == "Extensometer" && "Extensometer" in sensors

        # measurement is from extensometer
        tracking.kf.H = [1.0 0.0 0.0; 0.0 0.0 1.0] # Measurement matrix ground based
        tracking.kf.R = diagm([100.0867^2, 100.3622^2]) # Measurement covariance ground based

        if !tracking.is_initialized
            println("Kalman Filter Initialization")
            # set the state with the initial location and zero velocity
            tracking.kf.x[2] = measurement[4][1]
            tracking.kf.x[3] = measurement[4][2]

            tracking.previous_timestamp = datetime2unix(measurement[1])
            tracking.is_initialized = true
            return
        end

    else

        #println("Measurement type not recognized")
        return

    end

    tracking.previous_timestamp = datetime2unix(measurement[1])

    # predict
    predict(tracking.kf)

    # measurement update
    update(tracking.kf, vec(measurement[4]))

    return measurement[1]

end



function process_measurement_vel(tracking::Tracking, measurement::Vector, sensors::Vector)

    if measurement[2] == "Total Station" && "Total Station" in sensors
        # measurement is from total station
        tracking.kf.H = [1.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0] # Measurement matrix total station
        tracking.kf.R = diagm([4.7384^2, 7.9214^2, 11.02575^2]) # Measurement covariance total station

        if !tracking.is_initialized
            println("Kalman Filter Initialization")
            # set the state with the initial location and zero velocity
            tracking.kf.x[1] = measurement[4][1]
            tracking.kf.x[2] = measurement[4][2]
            tracking.kf.x[3] = measurement[4][3]

            tracking.previous_timestamp = datetime2unix(measurement[1])
            tracking.is_initialized = true
            return
        end

    elseif measurement[2] == "SAT-INSAR" && "SAT-INSAR" in sensors

        # measurement is from satellite
        tracking.kf.H = [1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0] # Measurement matrix satellite
        tracking.kf.R = diagm([12.5^2, 10.0^2]) # Measurement covariance satellite

        if !tracking.is_initialized
            println("Kalman Filter Initialization")
            # set the state with the initial location and zero velocity
            tracking.kf.x[1] = measurement[4][1]
            tracking.kf.x[3] = measurement[4][2]

            tracking.previous_timestamp = datetime2unix(measurement[1])
            tracking.is_initialized = true
            return
        end

    elseif measurement[2] == "GB-INSAR" && "GB-INSAR" in sensors

        # measurement is from ground based
        tracking.kf.H = [0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0] # Measurement matrix ground based
        tracking.kf.R = diagm([9.9869^2, 0.3061^2]) # Measurement covariance ground based

        if !tracking.is_initialized
            println("Kalman Filter Initialization")
            # set the state with the initial location and zero velocity
            tracking.kf.x[2] = measurement[4][1]
            tracking.kf.x[3] = measurement[4][2]

            tracking.previous_timestamp = datetime2unix(measurement[1])
            tracking.is_initialized = true
            return
        end

    elseif measurement[2] == "Extensometer" && "Extensometer" in sensors

        # measurement is from extensometer
        tracking.kf.H = [1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0] # Measurement matrix ground based
        tracking.kf.R = diagm([1.0867^2, 0.3622^2]) # Measurement covariance ground based

        if !tracking.is_initialized
            println("Kalman Filter Initialization")
            # set the state with the initial location and zero velocity
            tracking.kf.x[2] = measurement[4][1]
            tracking.kf.x[3] = measurement[4][2]

            tracking.previous_timestamp = datetime2unix(measurement[1])
            tracking.is_initialized = true
            return
        end

    else

        #println("Measurement type not recognized")
        return

    end

    dt = (datetime2unix(measurement[1]) - tracking.previous_timestamp)
    tracking.previous_timestamp = datetime2unix(measurement[1])
    previous_state = tracking.kf.x

    # predict
    predict(tracking.kf)

    # measurement update
    update(tracking.kf, vec(measurement[4]))
    # tracking.kf.x[4:6] = (previous_state[1:3] - tracking.kf.x[1:3])/dt

    return measurement[1]
end



function process_measurement_global(tracking::Tracking, measurement::Vector)

    
    if !tracking.is_initialized
        println("Kalman Filter Initialization")
        tracking.previous_timestamp = datetime2unix(measurement[1])
        tracking.is_initialized = true
        return
    end



    dt = (datetime2unix(measurement[1]) - tracking.previous_timestamp)
    tracking.previous_timestamp = datetime2unix(measurement[1])
    previous_state = tracking.kf.x

    # predict
    predict(tracking.kf)

    # measurement update
    update(tracking.kf, measurement[2])
    # tracking.kf.x[4:6] = (previous_state[1:3] - tracking.kf.x[1:3])/dt

    return measurement[1]
end