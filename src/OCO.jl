
function produce_mode_query(modes::Vector{String}) :: String

    mode_query = ""
    if length(modes) > 0
        mode_query = "AND ("
        for j in 1:length(modes)
            mode_query *= "mode like '%$(modes[j])%' "

            if j < length(modes)
                mode_query *= "OR "
            end

        end
        mode_query *= ")"
    end

    return mode_query

end

function extract_OCO_sounding_id(sounding_id::Int)
    return parse(Int, string(sounding_id)[end])
end

function convert_OCO_sounding_id_to_fp(sounding_id::Int)
    return(parse(Int, string(sounding_id)[end]))
end

function convert_OCO_sounding_id_to_date(sounding_id::Int)

    # Conert to String
    snid = string(sounding_id)

    # OCO-type sounding ID must be 16 digits long
    if length(snid) != 16
        @error "Sounding ID $(sonding_id) must be 16 characters/digits long!"
        return nothingx
    end

    id_year = parse(Int, snid[1:4])
    id_month = parse(Int, snid[5:6])
    id_day = parse(Int, snid[7:8])
    id_hour = parse(Int, snid[9:10])
    id_minute = parse(Int, snid[11:12])
    id_second = parse(Int, snid[13:14])
    id_milli = parse(Int, snid[15:16])

    return DateTime(id_year, id_month, id_day,
                    id_hour, id_minute, id_second, id_milli)
end



# These dictionaries contain the footprint-dependent coefficients
# to produce the SIF single-sounding uncertainty of the form
# sqrt(A + B * continuum_radiance)

OCO2_noise_dict = Dict(
    1 => (0.00243, 0.00192),
    2 => (0.00204, 0.00179),
    3 => (0.00365, 0.00174),
    4 => (0.00331, 0.0018) ,
    5 => (0.00329, 0.00188),
    6 => (0.00301, 0.00203),
    7 => (0.00241, 0.00211),
    8 => (0.00373, 0.00242)
)

OCO3_noise_dict = Dict(
    1 => (0.00016, 0.00161),
    2 => (0.00011, 0.00156),
    3 => (3.0e-5, 0.00144),
    4 => (3.0e-5, 0.00135),
    5 => (0.00011, 0.00131),
    6 => (0.00021, 0.0014),
    7 => (0.00037, 0.0015),
    8 => (0.00034, 0.00164)
)


function calculate_OCO_SIF_uncertainty(continuum::Real, instrument::String, fp::Int)
    if instrument == "OCO-3"
        noise_dict = OCO3_noise_dict
    elseif instrument == "OCO-2"
        noise_dict = OCO2_noise_dict
    else
        @error "Need to be either OCO-2 or OCO-3, but I got $(instrument)."
    end

    return sqrt(noise_dict[fp][1] + noise_dict[fp][2] * continuum)
end

function convert_OCO_df_to_scenes(df::DataFrame)

    # At this point, we expect all the instrument labels here
    # to be a unique

    unique_instrument = unique(df.instrument)

    if length(unique_instrument) != 1
        @error "We have multiple instruments via the SQLITE query: $(unique_instrument)"
        return nothing
    end

    instrument = unique_instrument[1]

    # Here we pre-load the surface data that we then sample
    # in the scene loop below. The strategy is clear: accessing
    # the surface data HDF file point-by-point is slow,
    # so we read in some intermediate objects here by only
    # loading a small subset that bounds the ROI. Once this is
    # in memory, we can very quickly sample and calculate reflectivity etc.

    lon_bound_min = minimum(df.lon)
    lon_bound_max = maximum(df.lon)
    lat_bound_min = minimum(df.lat)
    lat_bound_max = maximum(df.lat)

    vnp_sd = create_VNP_SD_from_locbounds(
        lon_bound_min,
        lat_bound_min,
        lon_bound_max,
        lat_bound_max
    )

    #=
    tropomi_data = create_TROPOMI_SIF_from_locbounds(
        lon_bound_min,
        lat_bound_min,
        lon_bound_max,
        lat_bound_max
    )
    =#

    locarray = Geolocation[]
    scenearray = Scene[]

    # Go through scenes and create scene and location objects
    for row in eachrow(df)

        this_footprint = convert_OCO_sounding_id_to_fp(row.sounding_id)

        this_loc = Geolocation(row.lon, row.lat)

        this_loctime = GeolocationTime(
            this_loc,
            convert_OCO_sounding_id_to_date(row.sounding_id)
        )

        # New unction for this particular footprint
        this_ucert_function(x) = calculate_OCO_SIF_uncertainty(
            x,
            instrument,
            this_footprint
        )

        this_scene = create_SIF_scene(
            instrument,
            row.mode,
            this_loctime,
            row.vza,
            this_ucert_function,
            vnp_sd
        )

        push!(locarray, this_loc)
        push!(scenearray, this_scene)
    end

    @info "Populated $(length(locarray)) scene objects."

    # We must time-order scenes (measurements)
    time_sort = sortperm(get_time.(scenearray))
    # Not sure what order we want for locations?

    return instrument, locarray[time_sort], scenearray[time_sort]

end



"""
    OCOSampling(target_lon, target_lat, radius, oco_location_file)

Constructor for an OCOSampling type object. Reads files from "datapath"
and subsets those scenes by only retain those locations which are
within the radius of the user-supplied target lon/lat location.

# Arguments
- 'target_lon::Number': Longitude of target location
- 'target_lat::Number': Latitude of target location
- 'radius::Number': Radius [km] within target location where scenes will be kept
- 'oco_location_file::String': Path to location of OCO-x locations file

# Returns
- 'OCOSampling': OCO sampling object

"""
function OCOSampling(target_lon::Real, target_lat::Real,
                     radius::Real, oco_db_file::String;
                     modes::Vector{String}=String[])


    # First - construct a sensible bounding box which will definitely contain the
    # locations in the circular search.

    wgs84 = Projection(Proj4.epsg[4326])
    xy = lonlat2xy([Float64(target_lon), Float64(target_lat)], wgs84)
    r = radius * 1000.0 * 1.02 # Add x percent as buffer, radius is given in [km] but we need [m]

    north = xy2lonlat(geod_destination(xy, 0.0, r, wgs84), wgs84)
    east = xy2lonlat(geod_destination(xy, 90.0, r, wgs84), wgs84)
    south = xy2lonlat(geod_destination(xy, 180.0, r, wgs84), wgs84)
    west = xy2lonlat(geod_destination(xy, 270.0, r, wgs84), wgs84)

    lon_min = minimum([west[1], south[1]])
    lat_min = minimum([west[2], south[2]])
    lon_max = maximum([east[1], north[1]])
    lat_max = maximum([east[2], north[2]])

    # Check for bounding box extent
    if (lon_min > lon_max)
        @error "Longitude maximum [$(lon_max)] is smaller than longitude minimum [$(lon_min)]"
        return nothing
    end

    if (lat_min > lat_max)
        @error "Latitude maximum [$(lon_max)] is smaller than latitude minimum [$(lat_min)]"
        return nothing
    end


    # If the user requests certain modes only, we inject the following additional query
    # into the main query.
    mode_query = produce_mode_query(modes)

    # Construct SQL query
    query = """
    SELECT sounding_id, vza, ST_X(geometry) AS lon, ST_Y(geometry) AS lat, mode
    FROM locations WHERE
    PtDistWithin (MakePoint($(target_lon), $(target_lat), 4326), geometry, $(radius * 1000.0)) = 1
    $(mode_query)
    AND rowid IN
    (
        SELECT rowid FROM SpatialIndex WHERE f_table_name = 'locations' AND
        search_frame = BuildMBR($(lon_min), $(lat_min), $(lon_max), $(lat_max))
    );
    """

    # Load up DB and execute query
    db = SQLite.DB(oco_db_file)
    SQLite.enable_load_extension(db, true)
    DBInterface.execute(db, "SELECT load_extension('/home/psomkuti/miniconda3/lib/mod_spatialite.so');")
    df = DBInterface.execute(db, query) |> DataFrame

    # Obtain instrument label from database and attach to dataframe
    this_instrument = SQLite.getvalue(DBInterface.execute(db, "SELECT * FROM instrument;"), 1, String)
    df[!, "instrument"] .= this_instrument

    # Turn DataFrame into location and scene arrays
    instrument, locarray, scenearray = convert_OCO_df_to_scenes(df)

    # Create OCOSampling object and return
    info = "OCO-type sampling with center point $(target_lon), $(target_lat) and radius $(radius) m."

    return OCOSampling(
        info,
        [instrument],
        locarray,
        scenearray
    )

end

"""
    OCOSampling(bbox, oco_location_file)

Constructor for an OCOSampling type object. Reads files from "datapath"
and subsets those scenes by only retain those locations which are
within the bounding box defined as [lon_min, lat_min, lon_max, lat_max]

# Arguments
- 'bbox::Array{Number, 1}': Bounding box [lon_min, lat_min, lon_max, lat_max]
- 'oco_db_file::String': Path to location of OCO-x locations SQLite file

# Returns
- 'OCOSampling': OCO sampling object

"""
function OCOSampling(bbox::Vector{<:Real}, oco_db_file::String)

    # Check for validity of bounding box length
    if length(bbox) != 4
        @error "Length of bounding box vector must be 4!"
        return nothing
    end

    lon_min, lat_min, lon_max, lat_max = bbox

    # Check for bounding box extent
    if (lon_min > lon_max)
        @error "Longitude maximum [$(lon_max)] is smaller than longitude minimum [$(lon_min)]"
        return nothing
    end

    if (lat_min > lat_max)
        @error "Latitude maximum [$(lon_max)] is smaller than latitude minimum [$(lat_min)]"
        return nothing
    end

    # Construct SQL query
    query = """
    SELECT sounding_id, vza, ST_X(geometry) AS lon, ST_Y(geometry) AS lat, mode
    FROM locations WHERE
    MBRWithin(geometry, BuildMBR($(lon_min), $(lat_min), $(lon_max), $(lat_max))) = 1
    AND rowid IN
    (
        SELECT rowid FROM SpatialIndex WHERE f_table_name = 'locations' AND
        search_frame = BuildMBR($(lon_min), $(lat_min), $(lon_max), $(lat_max))
    );
    """

    # Load up DB and execute query
    db = SQLite.DB(oco_db_file)
    SQLite.enable_load_extension(db, true)
    DBInterface.execute(db, "SELECT load_extension('/home/psomkuti/miniconda3/lib/mod_spatialite.so');")
    df = DBInterface.execute(db, query) |> DataFrame

    # Turn DataFrame into location and scene arrays

    instrument, locarray, scenearray = convert_OCO_df_to_scenes(df)

    # Create OCOSampling object and return
    info = "OCO-type sampling with bounding box $(bbox) with N=$(size(df)[1])."

    return OCOSampling(
        info,
        [instrument],
        locarray,
        scenearray
    )

end
