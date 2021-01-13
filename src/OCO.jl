"""
    read_OCO_file(fname)

Reads the relevant contents of an OCO-2/3 L1bSc or GeoSc file and
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

    return out
end

"""
    list_OCO_files(datapath)

Looks for L1bSc and GeoSc files within either a directory, or from a
list of filenames.

# Arguments
- 'datapath::String': Folder possibly containing L1bSc and GeoSc files

# Returns
- 'String[]' : Array of strings containing the file paths
"""
function list_OCO_files(datapath::String)

    # "datapath" can either be a directory, so every *.h5 file will be
    # read in, as long as it has either L1bSc or GeoSc in its name.

    filelist = String[]

    if isfile(datapath)
        error("Not implemented yet!")
    elseif ispath(datapath)
        # We have been given a path, so scan it for
        # GeoSc and L1bSc files.
        for pattern in ["*GeoSc*.h5", "*L1bSc*.h5"]
            for tmp_file in glob(pattern, datapath)
                push!(filelist, tmp_file)
            end
        end
    end

    return filelist

end



"""
    OCOSampling(target_lon, target_lat, radius, datapath)

Constructor for an OCOSampling type object. Reads files from "datapath"
and subsets those scenes by only retain those locations which are
within the radius of the user-supplied target lon/lat location.

# Arguments
- 'target_lon::Number': Longitude of target location
- 'target_lat::Number': Latitude of target location
- 'radius::Number': Radius [km] within target location where scenes will be kept
- 'datapath::String': Path to folder with L1bSc/GeoSc files, or path to text file with a list of those files

# Returns
- 'OCOSampling': OCO sampling object

"""
function OCOSampling(target_lon::Number, target_lat::Number,
                     radius::Number, oco_location_file::String)

    locarray = Geolocation[]
    scenearray = Scene[]

    tmp_dict = read_OCO_file(oco_location_file)

    # Go through scenes and create scene and location objects
    for i in 1:length(tmp_dict["lon"])
        this_loc = Geolocation(
            tmp_dict["lon"][i],
            tmp_dict["lat"][i]
        )

        this_loctime = GeolocationTime(
            this_loc,
            tmp_dict["datetime"][i]
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
            tmp_dict["sza"][i],
            tmp_dict["vza"][i],
            this_nirv,
            this_albedo
        )
        push!(locarray, this_loc)
        push!(scenearray, this_scene)
    end

    # Select only locations which lie within radius
    mask = mask_locations_within_radius(
        target_lon,
        target_lat,
        radius,
        locarray
    )

    # Create OCOSampling object and return

    info = "OCO-type sampling with radius $(radius)km centered at ($(target_lon), $(target_lat))."

    return OCOSampling(
        info,
        locarray[mask],
        scenearray[mask]
    )

end
