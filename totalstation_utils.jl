using DataFrames
using XLSX
using FileIO
using ExcelFiles
using Dates
using Missings
using Rotations
using AlgebraOfGraphics
using LinearAlgebra
using GLMakie
using ColorTypes
using GLMakie, Random, Colors, LinearAlgebra
using GeometryBasics: Cylinder, Pyramid
using Makie
using Statistics
using GeometryBasics


"""
import_totalstation(path) -> (DataFrame, DataFrame, DataFrame, DataFrame)

Imports the totalstation data set from the given path

# Arguments
- `path`: Absolute Path to the XLSX file containing the data set.

# Returns
- `(DataFrame, DataFrame, DataFrame, DataFrame)`: Tuple of DataFrames with data for all prisms, Dates, East, North and Height
"""
function import_totalstation(path)

    COL_NAMES = ["t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10", "t11", "t12", "t13", "t14", "t15", "t16", "t17", "t18", "ref1", "ref2", "ref3", "ref4"];

    ts_df_east =  DataFrame(XLSX.readtable(path, "East")...);
    ts_df_north = DataFrame(XLSX.readtable(path, "North")...);
    ts_df_height = DataFrame(XLSX.readtable(path, "Height")...);
    ts_df_date = DataFrame(XLSX.readtable(path, "Date")...);

    rename!(ts_df_east, COL_NAMES);
    rename!(ts_df_north, COL_NAMES);
    rename!(ts_df_height, COL_NAMES);
    rename!(ts_df_date, ["Date"]);

    println("Dimensions:\n\tEast: ", size(ts_df_east), "\n\tNorth: ", size(ts_df_north), "\n\tHeight: ", size(ts_df_height), "\n\tDates: ", size(ts_df_date))

    return (ts_df_date, ts_df_east, ts_df_north, ts_df_height)

end


"""
    prism2df(name, size, replace_nan) -> DataFrame

Prepares data from a total station prism to be used in a DataFrame

# Arguments
- `name`: Name of the prism
- `size`: Size of the DataFrame (no of rows)
- `replace_nan`: whether or not to replace NaN values with missing
- `replace_nan`: whether or not to drop rows with NaN values
- `ts_df_date`: DataFrame with raw dates
- `ts_df_east`: DataFrame with raw east data
- `ts_df_north`: DataFrame with raw north data
- `ts_df_height`: DataFrame with raw height data

# Returns
- `DataFrame`: Data frame for the given prism with DateTime, East, North and Height columns
"""
function prism2df(name, size, replace_nan, drop_nan, ts_df_date, ts_df_east, ts_df_north, ts_df_height)

    df = DataFrame()
    df[:, "Date"] = ts_df_date[!, :Date];
    df[:, "East"] = ts_df_east[1:size, name];
    df[:, "North"] = ts_df_north[1:size, name];
    df[:, "Height"] = ts_df_height[1:size, name];

    if drop_nan
        df = df[(df.East .!= "NaN"), :]
    end

    for col in names(df[!, 2:end])        
        df[!, col] = parse.(Float64, df[!, col])                            # Parse to Float64
        # df[!, col] = (df[!, col] .- df[!, col][1]).*100                   # Rescale to start at 0
        df[!, col] = df[!, col].*1000                                       # scale to mm
        if replace_nan
            df[!, col] = [ isnan(val) ? missing : val for val in df[!, col] ]   # Replace Nan with missing
        end
    end

    # Parse Date Column
    date_format = DateFormat("\'dd.mm.yyyy HH:MM\'")
    df.Date = DateTime.(df.Date, date_format)

    return df
end



"""
    kalman_filter_1d(x0, P0, A, B, Q, H, R, y) -> Array{Float64}

Computes a series of filtered values for a one-dimensional system using a Kalman Filter with the following definition:
**Model**<br>
Model without control input (i.e. ``B = \\begin{bmatrix}0 & 0\\\\0&0\\end{bmatrix}``):<br>
``x_k = A x_{k-1} + \\omega_k``<br>
``y_k = H x_k + v_k`` <br>
where ``\\omega_k \\sim N(0, Q)`` and ``\\nu_k \\sim N(0,R)``
<br>
**Prediction Step**<br>
``\\hat{x}_k = A \\hat{x}_{k-1} + B u_k``<br>
``P_k = A p_{k-1} A^T``
<br>
**Update Step**<br>
``G_k = P_k C^T (C p_k C^T + r)^{-1}``<br>
``\\hat{x}_k \\leftarrow \\hat{x}_k + G_k (z_k - C \\hat{x}_k)``<br>
``P_k \\leftarrow (I - G_k C) P_k``<br>

# Arguments
- `x0`: initial system state
- `P0`: covariance matrix of estimation process
- `A`: System state transition matrix
- `B`: Control input matrix
- `Q`: Covariance matrix of system noise
- `R`: Covariance matrix of measurement noise
- `y`: measurements/observations

# Returns
- `Array{Float64}`: filtered values, same length as measurements/observations
"""
function kalman_filter_1d(x0::Float64, P0::Float64, Q::Float64, H::Float64, R::Float64, y::Vector{Union{Missing,Float64}})

    N = length(y);
    x_k = x0
    P_k = P0

    filtered_values = zeros(0)

    for i in 1:N

        y_k = y[i]

        # Prediction Step
        x_k = x_k 
        P_k = P_k + Q

        # Update Step
        G_k = (P_k*H) / (H*P_k*H + R)
        x_k = x_k + G_k*(y_k - H*x_k)
        P_k = (1 - G_k*H) * P_k

        append!(filtered_values, x_k)

    end

    println(length(filtered_values), " predictions.");

    return filtered_values
end


"""
 rotmat2rotang(R) -> DataFrame

Computes the rotation angles around the x y and z axis as defined by the rotation matrix R

# Arguments
- `R`: 3x3 Rotation Matrix

# Returns
- `DataFrame`: rotation angles in order about x, y and z, in degree
"""
function rotmat2rotang(R::RotMatrix3{Float64})

    phi = atand(R[3,2], R[3,3])
    psi = atand(R[2,1], R[1,1])

    if cosd(psi) == 0
        theta = atand(-R[3,1], (R[2,1]/sind(psi)))
    else
        theta = atand(-R[3,1], (R[1,1]/cosd(psi)))
    end

    angles = DataFrame(roll=phi, pitch=theta, yaw=psi)

    return angles
end



"""
 ts_rototrans(ts1,ts2,ts3) -> DataFrame

Computes the roto-translation of a surface described by three vertices of which the time series are given.
The input time series must have the following columns: East, North, Height, Date

# Arguments
- `ts1`: DataFrame with Time Series of Vertex 1
- `ts2`: DataFrame Time Series of Vertex 2
- `ts3`: DataFrame Time Series of Vertex 3

# Returns
- `DataFrame`: translation along East, North and Height, Normal vector of the computed surface with columns Date, East, North, Height and Orientation
"""
function ts_rototrans(ts1::DataFrame, ts2::DataFrame, ts3::DataFrame)

    # Filter incomplete rows as data for all three points is needed for any time step
    filter1 = completecases(ts1)
    filter2 = completecases(ts2)
    filter3 = completecases(ts3)
    filter_all = filter1 .& filter2 .& filter3
    
    ts1_filtered = ts1[filter_all , :]
    ts2_filtered = ts2[filter_all , :]
    ts3_filtered = ts3[filter_all , :]

    # Store position and normal vectors
    estimated_positions = DataFrame()
    normal_vectors = []

    # For every time step with full information compute translation and normal vector
    for i in 1:size(ts1_filtered, 1)

        # First Vertices
        v2 = [ts1_filtered.East[i], ts1_filtered.North[i], ts1_filtered.Height[i]]
        v1 = [ts2_filtered.East[i], ts2_filtered.North[i], ts2_filtered.Height[i]]
        v3 = [ts3_filtered.East[i], ts3_filtered.North[i], ts3_filtered.Height[i]]

        # Compute centroid
        centroid = (v1 + v2 + v3) / 3

        # Compute normal vector
        edge1 = v2 - v1
        edge2 = v3 - v1

        # Unit normal vector
        normal_vec = normalize(cross(edge1 , edge2 ))

        push!(estimated_positions, (Date=ts1_filtered.Date[i], East=centroid[1], North=centroid[2], Height=centroid[3]))
        push!(normal_vectors, normal_vec)

    end

    println(size(estimated_positions), " estimated positions\n", size(normal_vectors), " estimated normal vectors")

    # Compute cumulated differences for better visualisation
    cum_diffs = DataFrame()
    cum_diffs[:, "Date"] = estimated_positions[2:end, "Date"];
    cum_diffs[:, "East"] =  cumsum(diff(estimated_positions[!, "East"]));
    cum_diffs[:, "North"] =  cumsum(diff(estimated_positions[!, "North"]));
    cum_diffs[:, "Height"] =  cumsum(diff(estimated_positions[!, "Height"]));
    cum_diffs[:, "Orientation"] = normal_vectors[2:end];

    return (cum_diffs, estimated_positions)

end



"""
plot_rototranslation(df)

Visualises the translation and rotation of an object in 3 dimensions over time.
The DataFrame is expected to have the columns Date, East, North, Height and Orientation 
(the orientation being a 3 element vector).

# Arguments
- `df`: DataFrame containing the roto-translation data
"""
function plot_rototranslation(df::DataFrame)

    GLMakie.activate!()
    Makie.inline!(false);

    # Define the two vectors that indicate the orientation of the surface (replace with your own vectors)
    # Use mean values over 200 normal vectors to withstand the influence of apparent noise
    v1 = mean(df[1:400, "Orientation"])
    v2 = mean(df[end-400:end, "Orientation"])

    # Compute the rotation matrix that rotates from v1 to v2
    q = rotation_between(v1, v2)

    # Compute the rotation matrix
    aa = AngleAxis(q)
    R = RotMatrix(aa)

    # Define the rotation matrix and initial vectors
    v1 = [-100, 0, 0]  # Initial x-axis vector (inverted x-axis)
    v2 = [0, 100, 0]  # Initial y-axis vector
    v3 = [0, 0, 100]  # Initial z-axis vector

    # Apply the rotation matrix to the initial vectors
    v1_new = R * v1
    v2_new = R * v2
    v3_new = R * v3

    n = size(df.Date,1)

    arrow_size = Vec3f(20, 20, 50)
    tip_radius = 2
    tip_length = 10

    startpoint = Point3f(0,0,0)
    endpoint = Point3f(df.East[end],df.North[end],df.Height[end])

    for i in 1:100

        v1_new = R * v1_new
        v2_new = R * v2_new
        v3_new = R * v3_new

    end

    f = Figure(size = (800, 800))
    ax = Axis3(
        f[1, 1],
        title="System Translation and Rotation",
        xlabel = "East Displacement [mm]",
        ylabel = "North Displacement [mm]",
        zlabel = "Height Displacement [mm]",
        aspect=(1, 1, 1),
        perspectiveness=0.5,
    )

    hm = lines!(ax, df.East, df.North, df.Height, color=1:n, linewidth=3)
    ax2 = Colorbar(f[1, 2], hm; label="Date", vertical=true,
        flipaxis=true, ticksize=15, tickalign=1, width=20)
    ax2.tickformat = x -> Dates.format.(df.Date[Int.(x)], "mm-yyyy") 

        
    # Plot coordinate systems to visualise orientation
    arrows3d!([startpoint], [Vec3f(v1)], color=:red)
    arrows3d!([startpoint], [Vec3f(v2)], color=:green)
    arrows3d!([startpoint], [Vec3f(v3)], color=:blue)

    arrows3d!([endpoint], [Vec3f(v1_new)], color=:red)
    arrows3d!([endpoint], [Vec3f(v2_new)], color=:green)
    arrows3d!([endpoint], [Vec3f(v3_new)], color=:blue)

    wireframe!(ax, Rect3f(Point3f(df.East[1000]-2, df.North[1000]-20, df.Height[1000]-20), Vec3f(2,40,40)), color = :red)
    wireframe!(ax, Rect3f(Point3f(df.East[2000]-2, df.North[2000]-20, df.Height[2000]-20), Vec3f(2,40,40)), color = :red)
    wireframe!(ax, Rect3f(Point3f(df.East[3000]-2, df.North[3000]-20, df.Height[3000]-20), Vec3f(2,40,40)), color = :red)
    wireframe!(ax, Rect3f(Point3f(df.East[4000]-2, df.North[4000]-20, df.Height[4000]-20), Vec3f(2,40,40)), color = :red)

    

    # Set axis limits to have same scale
    xlims!(ax, (-750, 50))
    ylims!(ax, (-400, 400))
    zlims!(ax, (-750, 50))

    f
end



"""
rotmat2euler(rotmats::Array{Float64})

Computes the Euler angles from a rotation matrix.

# Arguments
- `rotmats`: Array of rotation matrices
"""
function rotmat2euler(rotmat)

    euler_angles = zeros(3)

    # Check for gimbal lock
    if abs(rotmat[3,1]) == 1
            
        euler_angles[1] = atan(-rotmat[3,1]*rotmat[1,2], -rotmat[3,1]*rotmat[1,3])
        euler_angles[2] = -rotmat[3,1]*pi/2
        euler_angles[3] = 0

    else

        # First solution
        euler_angles[2] = -asin(rotmat[3,1])
        euler_angles[1] = atan(rotmat[3,2] / cos(euler_angles[2]) , rotmat[3,3] / cos(euler_angles[2]))
        euler_angles[3] = atan(rotmat[2,1] / cos(euler_angles[2]) , rotmat[1,1] / cos(euler_angles[2]))

        # Second solution
        # euler_angles[2] = pi - euler_angles[2]
        # euler_angles[1] = atan(rotmat[3,2] / cos(euler_angles[2]) , rotmat[3,3] / cos(euler_angles[2]))
        # euler_angles[3] = atan(rotmat[2,1] / cos(euler_angles[2]) , rotmat[1,1] / cos(euler_angles[2]))
    
    end

    return euler_angles


end