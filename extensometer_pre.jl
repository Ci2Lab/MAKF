using CSV
using Missings
using Dates
using DataFrames
using Statistics

function import_extensometer_data(file_path::String)

    df = CSV.read(file_path, DataFrame)
    df.Date = map(parse_date, df.Date);

    return df

end


function parse_date(date_str::AbstractString)
    return DateTime(date_str, "dd.mm.yyyy HH:MM:SS")
end


function filter_extensometer_data(df::DataFrame, start_date::DateTime, end_date::DateTime)

    # Filter time series by date and align new start to zero
    df = df[(df.Date .> positions[1, "Date"]) .& (df.Date .< positions[end, "Date"]), :]

    # Drop data points with impossible values
    df = dropmissing!(ifelse.(df.Distance .< 0, missing, df))
    df = dropmissing!(ifelse.(df.Distance .> 1500, missing, df))

    # Compute acceleration for each point
    df.Acceleration = [0.0; diff(df.Distance) ./ diff(datetime2unix.(df.Date))] # unix timestamps in seconds

    # Drop large accelerations
    df = dropmissing!(ifelse.(abs.(df.Acceleration) .> 0.0005, missing, df))

    return df

end


function decompose_los_measurements(df::DataFrame, col_name_los::String, view_angle::Float64)

    df.Vertical = df[:, col_name_los] .* sin.(view_angle * pi / 180)
    df.Easting = df[:, col_name_los] .* cos.(view_angle * pi / 180)

    return df

end


function get_level_zero_extensometer(file_path::String, start_date::DateTime, end_date::DateTime)

    df = import_extensometer_data(file_path)
    df = filter_extensometer_data(df, start_date, end_date)
    df = decompose_los_measurements(df, "Distance", 37.0)

    return df

end