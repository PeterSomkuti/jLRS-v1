"""
    read_OCO_file(fname)

Reads the relevant contents of an OCO-2/3 collection file and
returns them as a dict.

# Arguments
- 'fname::String': Location of L1bSc or GeoSc (HDF5) file

# Returns
- 'Dict' : Dictionary containing the contents of the L1bSc or GeoSc file.
"""
function read_OCO_file(fname::String)

    out = Dict()

    h5open(fname, "r") do h5

        out["lon"] = h5["lon"][:]
        out["lat"] = h5["lat"][:]
        out["sza"] = h5["sza"][:]
        out["vza"] = h5["vza"][:]
        out["mode"] = h5["mode"][:]

        # Convert tai93 to DateTime objects
        out["datetime"] = DateTime(1993,1,1) .+ (Dates.Microsecond).(h5["tai93"][:] .* 1e6)
    end

    # In case we need to log whether this is OCO-2 or OCO-3
    out["instrument"] = h5read(fname, "/instrument")

    return out
end


function covert_OCO_dict_to_scenes(dict)

    locarray = Geolocation[]
    scenearray = Scene[]

    # Go through scenes and create scene and location objects
    for i in 1:length(dict["lon"])
        this_loc = Geolocation(
            dict["lon"][i],
            dict["lat"][i]
        )

        this_loctime = GeolocationTime(
            this_loc,
            dict["datetime"][i]
        )

        # Replace this with samplers
        this_sif = rand()
        this_sif_ucert = rand()
        this_nirv = rand()
        this_albedo = rand()

        this_scene = Scene(
            this_loctime,
            this_sif,
            this_sif_ucert,
            dict["sza"][i],
            dict["vza"][i],
            this_nirv,
            this_albedo
        )
        push!(locarray, this_loc)
        push!(scenearray, this_scene)
    end

    # We must time-order scenes (measurements)
    time_sort = sortperm((p -> p.loctime.time).(scenearray))
    # Not sure what order we want for locations?
    

    return dict["instrument"], locarray, scenearray[time_sort]

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
function OCOSampling(target_lon::Number, target_lat::Number,
                     radius::Number, oco_location_file::String)

    tmp_dict = read_OCO_file(oco_location_file)
    instrument, locarray, scenearray = covert_OCO_dict_to_scenes(tmp_dict)

    # Select only locations which lie within radius
    mask = mask_locations_within_radius(
        target_lon,
        target_lat,
        radius,
        locarray
    )

    # Create OCOSampling object and return
    info = "OCO-type sampling with radius $(radius)km centered at ($(target_lon), $(target_lat)) with N=$(sum(mask))"

    return OCOSampling(
        info,
        [instrument],
        locarray[mask],
        scenearray[mask]
    )

end

"""
    OCOSampling(bbox, oco_location_file)

Constructor for an OCOSampling type object. Reads files from "datapath"
and subsets those scenes by only retain those locations which are
within the bounding box defined as [lon_min, lat_min, lon_max, lat_max]

# Arguments
- 'bbox::Array{Number, 1}': Bounding box [lon_min, lat_min, lon_max, lat_max]
- 'oco_location_file::String': Path to location of OCO-x locations file

# Returns
- 'OCOSampling': OCO sampling object

"""
function OCOSampling(bbox::Array{<:Number, 1}, oco_location_file::String)

    tmp_dict = read_OCO_file(oco_location_file)
    instrument, locarray, scenearray = covert_OCO_dict_to_scenes(tmp_dict)

    # Select only locations which lie within bounding box
    lons = (p -> p.lon).(locarray)
    lats = (p -> p.lat).(locarray)

    mask = (
        (lons .>= bbox[1]) .&
        (lons .<= bbox[3]) .&
        (lats .>= bbox[2]) .&
        (lats .<= bbox[4])
    )

    # Create OCOSampling object and return
    info = "OCO-type sampling with bounding box $(bbox) with N=$(sum(mask))."

    return OCOSampling(
        info,
        [instrument],
        locarray[mask],
        scenearray[mask]
    )

end
