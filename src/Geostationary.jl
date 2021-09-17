function create_geostationary_granule(
    target_lon::Real,
    target_lat::Real,
    central_longitude::Real,
    N_step=35,
    N_pixel=1016,
    x_space=5400, # footprint size E-W in m
    y_space=2700 # footprint size N-S in m
)

    N = N_step * N_pixel

    # lon/lat projection object
    wgs84 = Projection("+proj=longlat +datum=WGS84 +no_defs")
    # geostationary projection
    geo_proj = Projection("+proj=geos +lon_0=$(central_longitude) +h=35785831.0 +sweep=y")

    # Produce a rectangular grid given the parameters
    c_x = 0.0
    c_y = 0.0
    try
        c_x, c_y = Proj4.transform(wgs84, geo_proj, [target_lon, target_lat])
    catch
        return nothing
    end


    # Move up to the "top left" scene of the grid
    topleft_x = c_x + (N_step - 1) / 2 * x_space
    topleft_y = c_y + (N_pixel - 1) / 2 * y_space

    x_coords = Float64[]
    y_coords = Float64[]
    dts = DateTime[]

    for xstep in 0:N_step-1
        for ystep in 0:N_pixel-1

            this_x = topleft_x - x_space * xstep
            this_y = topleft_y - y_space * ystep

            push!(x_coords, this_x)
            push!(y_coords, this_y)

        end
    end

    lons = Float64[]
    lats = Float64[]

    for i in 1:N

        this_lon = NaN
        this_lat = NaN

        try
            this_lon, this_lat =
                Proj4.transform(geo_proj, wgs84, [x_coords[i], y_coords[i]])

        catch

        end

        push!(lons, this_lon)
        push!(lats, this_lat)

    end

    locarray = Geolocation[]

    # Produce the location array
    # These are the full granule!
    for i in 1:N
        this_loc = Geolocation(lons[i] ,lats[i])
        push!(locarray, this_loc)
    end

    return geo_proj, locarray

end


function calculate_geocarb_uncertainty(continuum::Real)

    #a = 833661.0
    #b = 4109.0

    # THIS IS RELATIVE UNCERTAINTY!
    # To get absolute uncertainty for 757nm you must multiply
    # by the reflected continuum level radiance!
    #return 1.0 / sqrt(a * albedo * cos(deg2rad(sza)) - b)

    a = 0.0004
    b = 0.0016

    return sqrt(a + continuum * b)

end


function GeostationaryIntensiveSampling(
    radius::Real,
    target_lon::Real,
    target_lat::Real,
    central_longitude::Real,
    start_time::DateTime,
    end_time::DateTime;
    N_step=35,
    N_pixel=1016,
    x_space=5400, # footprint size E-W in m
    y_space=2700, # footprint size N-S in m
    stare_time=Dates.Second(9) + Dates.Millisecond(600) # How long does one frame take?
)

    N = N_step * N_pixel
    instrument = "Geostationary $(central_longitude)"

    # This would make conversions from compound periods possible
    stare_time_period = Dates.Millisecond(Dates.toms(stare_time))


    geo_proj, locarray = create_geostationary_granule(
        target_lon,
        target_lat,
        central_longitude,
        N_step,
        N_pixel,
        x_space,
        y_space
    )

    lons = (p -> p.lon).(locarray)
    lats = (p -> p.lat).(locarray)

    lon_bound_min = minimum(lons)
    lon_bound_max = maximum(lons)
    lat_bound_min = minimum(lats)
    lat_bound_max = maximum(lats)

    vnp_sd = create_VNP_SD_from_locbounds(
        lon_bound_min,
        lat_bound_min,
        lon_bound_max,
        lat_bound_max
    )


    scenearray = Scene[]

    # repeat the loop over all locations as long as the measurement
    # time stays within the stated limit..

    current_frame_time = start_time
    current_granule = 1
    current_frame = 1

    while current_frame_time <= end_time

        #println(current_frame, "/", current_granule)

        for i in 1:N_pixel

            # Grab the location corresponding to pixel/step
            this_loc = locarray[i + (current_frame - 1) * N_pixel]

            this_loctime = GeolocationTime(
                this_loc,
                current_frame_time
            )

            # #############################
            # Calculate viewing zenith here
            # #############################
            #
            # Position of scene in ECEF
            r_location = geodetic_to_ecef(this_loc.lon, this_loc.lat, 0.0)
            # Position of satellite in ECEF
            r_satellite = geodetic_to_ecef(central_longitude, 0.0, 0.0)
            # Normalized location
            r_norm = normalize(r_location)
            r_loc_to_sat_norm = normalize(r_satellite - r_location)
            _tmp = dot(r_loc_to_sat_norm, r_norm)
            if (_tmp > 1) & (_tmp < 1 + 1e-6)
                _tmp = 1.0
            end

            this_vza = rad2deg(acos(_tmp))

            this_scene = create_SIF_scene(
                instrument,
                "N/A",
                this_loctime,
                this_vza,
                calculate_geocarb_uncertainty,
                vnp_sd
            )

            push!(scenearray, this_scene)

        end

        if current_frame < N_step
            current_frame += 1
        else
            current_frame = 1
            current_granule += 1
        end

        current_frame_time += stare_time_period
    end


    N_scene = length(scenearray)
    # Once the full scenes have been done, we need to subset
    # to the radius given by the user. The reason why we can
    # do this only AFTER all the scenes have been calculated,
    # is that we must have scenes corresponding to some real
    # scanning operation, and can only filter afterwards.

    radius_mask = Int[]
    for i in 1:N_scene
        if check_location_within_radius(target_lon, target_lat,
                                        radius, scenearray[i].loctime)

            push!(radius_mask, i)

        end
    end


    # Match the location array with the locations found in the
    # scenes. Drop any locations that don't appear in scenes.
    new_locarray = unique((p -> p.loctime.loc).(scenearray))

    info = "Geostationary intensive sampling at $(target_lon), $(target_lat)"

    return GeostationaryIntensiveSampling(
        info,
        [instrument],
        locarray,
        scenearray[radius_mask]
    )

end
