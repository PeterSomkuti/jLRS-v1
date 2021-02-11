
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

function convert_date_to_sounding_id(date::DateTime)
    
end


function extract_OCO_sounding_id(sounding_id::Int)
    return parse(Int, string(sounding_id)[end])
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
    locarray = Geolocation[]
    scenearray = Scene[]

    # Go through scenes and create scene and location objects
    for row in eachrow(df)

        this_loc = Geolocation(row.lon, row.lat)

        this_loctime = GeolocationTime(
            this_loc,
            convert_OCO_sounding_id_to_date(row.sounding_id)
        )

        # solar azimuth (saa) at this point unused!
        this_sza, this_saa = calculate_solar_angles(this_loctime)

        # Obtain irradiance information (downwelling radiace at surface)
        this_irradiance = calculate_BOA_irradiance(this_sza, 757.0)

        # Calculate PPFD for this scene
        this_PPFD = calculate_PPFD(this_sza)

        # Replace this with samplers
        this_sif = rand()
        this_sif_ucert = rand()

        this_nir = calculate_reflectance(this_loctime, this_sza, "NIR")
        this_vis = calculate_reflectance(this_loctime, this_sza, "VIS")
        this_ndvi = (this_nir - this_vis) / (this_nir + this_vis)
        this_nirv = this_ndvi * this_nir

        this_reflectance = calculate_reflectance(this_loctime, this_sza, "755")

        # At the moment no cloud or aerosol data
        this_od = 0.0

        this_scene = Scene(
            instrument,
            row.mode, # sampling mode comes from the instrument
            this_loctime, # location time comes from the instrument
            this_sif,
            this_sif_ucert,
            this_sza, # SZA comes from calculcations (via loctime)
            row.vza, # viewing zenith comes from the instrument
            this_ndvi, # NDVI from VIIRS
            this_nirv, # NIRv calculated from NDVI * NIR
            this_irradiance, # Irradiance at the surface and some ref. wl
            this_PPFD, # Integrated irradiance at PAR wavelengths
            this_reflectance, # Reflectance comes from BRDF sampling and SZA
            this_od # optical depth may come from ISCCP one day
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
