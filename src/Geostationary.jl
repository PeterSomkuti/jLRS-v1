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

            # solar azimuth (saa) at this point unused!
            this_sza, this_saa = calculate_solar_angles(this_loctime)

            # Obtain irradiance information (downwelling radiance at surface),
            # already SZA-corrected!
            # (unlike in L2 algorithms, where irradiance needs to be multiplied by
            #  mu0 to account for normal component)
            this_irradiance = calculate_BOA_irradiance(this_sza, 757.0)

            # Calculate PPFD for this scene
            this_PPFD = calculate_PPFD(this_sza)

            # black sky albedos from BRDFs?
            this_nir = calculate_reflectance(this_loctime, this_sza, "M7", vnp_sd)
            this_vis = calculate_reflectance(this_loctime, this_sza, "M5", vnp_sd)

            # These wavelengths are for VIIRS M5 and M7,
            # factor of 1000 takes us from W/m2/nm/sr to W/m2/um/sr
            this_nir_radiance = this_nir * calculate_BOA_irradiance(this_sza, 865.0) * 1000 / pi
            this_vis_radiance = this_vis * calculate_BOA_irradiance(this_sza, 672.0) * 1000 / pi

            this_ndvi = (this_nir - this_vis) / (this_nir + this_vis)

            this_nirv = this_ndvi * this_nir
            this_nirv_radiance = this_nirv * calculate_BOA_irradiance(this_sza, 865.0) * 1000 / pi 

            # Reflectance at ~757 nm is roughly between VIIRS bands M5 and M7
            refl_M5 = calculate_reflectance(this_loctime, this_sza, "M5", vnp_sd)
            refl_M7 = calculate_reflectance(this_loctime, this_sza, "M7", vnp_sd)

            this_reflectance = 0.5 * (refl_M5 + refl_M7)

            # irradiance is in /nm, but we want /um
            this_TOA_radiance = this_irradiance * this_reflectance * 1000 / pi

            # This is in W/m2/sr/um
            this_sif = model_fluorescence(
                this_PPFD,
                25.0,
                200.0,
                209.0,
                757.0, # wavelength in nm
                this_ndvi
            )

            # Estimate sigma based on empirical GeoCarb model
            this_sif_ucert = calculate_geocarb_uncertainty(this_TOA_radiance)

            # At the moment no cloud or aerosol data,
            # could add ISCCP sampler in here
            this_od = 0.0

            this_scene = Scene(
                instrument,
                "N/A", # sampling mode comes from the instrument
                this_loctime, # location time comes from the instrument
                this_sif,
                this_sif_ucert,
                this_sza, # SZA comes from calculcations (via loctime)
                0.0, # viewing zenith comes from the instrument,
                this_nir, # NIR reflectance
                this_vis, # VIS reflectance
                this_ndvi, # NDVI from VIIRS
                this_nirv, # NIRv calculated from NDVI * NIR
                this_nirv_radiance, # NIRv * L0(868nm)
                this_irradiance, # Irradiance at the surface and some ref. wl
                this_PPFD, # Integrated irradiance at PAR wavelengths 400nm to 700nm
                this_reflectance, # Reflectance comes from BRDF sampling and SZA
                this_od # optical depth may come from ISCCP one day
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
