tropomi_tle_string = """
       S5P
       1 42969U 17064A   21193.45435625 -.00000014  00000-0  14036-4 0  9992
       2 42969  98.7381 133.0236 0001130  90.2950 269.8356 14.19543373194087
       """



function intersect_with_wgs84!(
    origin,
    look_vector,
    output,
    a=6378137.0, # m
    b=6356752.314245 # m
    )

    x, y, z = origin
    u, v, w = look_vector

    t = -(1.0/(b^2 * (u^2 + v^2) + (a^2 * w^2))) * (b^2 * (u*x + v*y) +
        (a^2 * w*z) + 0.5*sqrt(4 * (b^2 * (u*x + v*y) + (a^2 * w*z))^2 -
        4 * (b^2 * (u^2 + v^2) + (a^2 * w^2)) * (b^2 * (-a^2 + x^2 + y^2) + (a^2 * z^2))))

    output[1] = x + t * u
    output[2] = y + t * v
    output[3] = z + t * w

end

function NadirSwathOrbiterSampling(
    target_lon::Real,
    target_lat::Real,
    radius::Real,
    start_time::DateTime,
    end_time::DateTime,
    tle_string_OR_kepler,
    field_of_view::Real,
    N_pixel::Integer,
    exposure_time::Real,
    uncertainty_function::Function;
    instrument::String="instrument"
    )

    locarray = Geolocation[]
    scenearray = Scene[]

    #wgs84_a = 6378137.0 # m
    #wgs84_b = 6356752.314245 # m

    #wgs84_a2 = wgs84_a^2
    #wgs84_b2 = wgs84_b^2

    # This is needed for whatever reason?
    eop = get_iers_eop()

    # Overall strategy:
    # First produce orbit, propagate orbit and generate
    # sounding locations, after which we can kick out
    # any scenes that are nighttime etc.

    sat_lons = Float32[]
    sat_lats = Float32[]

    sounding_lons = Float32[]
    sounding_lats = Float32[]

    # Check the type of 
    if tle_string_OR_kepler isa String

        tle_string = tle_string_OR_kepler
        tle = read_tle_from_string(tle_string)
        orbp = init_orbit_propagator(Val(:sgp4), tle[1])

    elseif tle_string_OR_kepler isa KeplerianElements

        k = tle_string_OR_kepler
        orbp = init_orbit_propagator(Val(:J4),
                                     k.t, # epoch
                                     k.a, # semi-major axis [m]
                                     k.e, # eccentricity
                                     k.i, # inclination [rad]
                                     k.Ω, # right ascension of ascending node [rad]
                                     k.ω, # argument of perigee / periapsis [rad],
                                     k.f  # true anomaly [rad] 
                                     )

    end

    # Array of angles (swath angles away from full nadir)
    theta_array_deg = collect(-N_pixel//2:1:N_pixel//2) .* (field_of_view / N_pixel) # TROPOMI: 216 pixels with 0.5 deg each
    theta_array = deg2rad.(theta_array_deg)
    cos_theta_array = cos.(theta_array)
    sin_theta_array = sin.(theta_array)

    # Rotation around x, y and z
    Rx = Matrix(1.0I, 3, 3)
    Ry = Matrix(1.0I, 3, 3)
    Rz = Matrix(1.0I, 3, 3)

    # temp matrices
    tmp1 = Matrix(1.0I, 3, 3)
    tmp2 = Matrix(1.0I, 3, 3)
    tmp3 = ones(3)
    tmp4 = ones(3)
    look_vector = ones(3)

    r_intersection = ones(3)

    # First - construct a sensible bounding box to load the surface data
    wgs84 = Projection(Proj4.epsg[4326])
    xy = lonlat2xy([Float64(target_lon), Float64(target_lat)], wgs84)
    r = radius * 1000.0 * 1.10 # Add x percent as buffer, radius is given in [km] but we need [m]

    north = xy2lonlat(geod_destination(xy, 0.0, r, wgs84), wgs84)
    east = xy2lonlat(geod_destination(xy, 90.0, r, wgs84), wgs84)
    south = xy2lonlat(geod_destination(xy, 180.0, r, wgs84), wgs84)
    west = xy2lonlat(geod_destination(xy, 270.0, r, wgs84), wgs84)

    lon_min = minimum([west[1], south[1]])
    lat_min = minimum([west[2], south[2]])
    lon_max = maximum([east[1], north[1]])
    lat_max = maximum([east[2], north[2]])

    vnp_sd = create_VNP_SD_from_locbounds(
        lon_min,
        lat_min,
        lon_max,
        lat_max
    )


    # this does not create a copy, however if you add a time period
    # it *will* create a new object, such that "start_time" will not
    # be modified.
    current_time = start_time

    pasted_swath_width = false

    while current_time <= end_time

        # Increment by exposure time
        current_time += Dates.Millisecond(exposure_time * 1000.0)

        # Turn the current time into a julian date
        current_epoch = date_to_jd(
            Dates.year(current_time),
            Dates.month(current_time),
            Dates.day(current_time),
            Dates.hour(current_time),
            Dates.minute(current_time),
            Dates.second(current_time) + Dates.millisecond(current_time) / 1000.0
        )

        # Propagate the orbit to the requested date
        r_teme, v_teme = propagate_to_epoch!(orbp, current_epoch)
        D_ITRF_TEME = r_eci_to_ecef(TEME(), ITRF(), current_epoch, eop)

        r_itrf = D_ITRF_TEME * r_teme
        x, y, z = r_itrf
        x2 = x*x
        y2 = y*y
        z2 = z*z

        v_itrf = D_ITRF_TEME * v_teme

        # Calculate the off-nadir footprints from the current
        # satellite position!
        this_lat, this_lon, this_alt = ecef_to_geodetic(r_itrf)
        sat_lon = rad2deg(this_lon)
        sat_lat = rad2deg(this_lat)

        # Normalized vector pointing from the satellite down
        # to center of the Earth
        r_sat_to_origin = normalize(-r_itrf)

        # Normalized velocity vector gives us flight direction
        v_norm = normalize(v_itrf)

        # Compute the rotation matrices for this point in time
        u, v, w = r_sat_to_origin # The vector to be rotated around rotation axis
        a, b, c = v_norm # The rotation axis
        d = sqrt(b^2 + c^2)

        # Rotations around x and y are independent of the swath angle
        Rx[2,2] = c / d
        Rx[3,3] = c / d
        Rx[2,3] = -b / d
        Rx[3,2] = b / d

        Ry[1,1] = d
        Ry[3,3] = d
        Ry[3,1] = a
        Ry[1,3] = -a

        @views tmp1[:,:] = Rx' * Ry'
        #=
        tmp1[1,1] = d
        tmp1[1,3] = a
        tmp1[2,1] = -a*b/d
        tmp1[2,2] = c/d
        tmp1[2,3] = b
        tmp1[3,1] = -a*c/d
        tmp1[3,2] = -b/d
        tmp1[3,3] = c
        =#


        @views tmp2[:,:] = Ry * Rx
        #=
        @views tmp2[:,:] = tmp1[:,:]
        tmp2[3,1] *= -1
        tmp2[3,2] *= -1
        tmp2[1,3] *= -1
        tmp2[2,3] *= -1
        =#

        # tmp3 = Ry * Rx * r_sat_to_origin
        @views tmp3[:] = tmp2 * r_sat_to_origin
        #=
        tmp3[1] = u*d - a*b*v/d - a*c*w/d #u*d - a*w
        tmp3[2] = c*v/d - b*w/d #-a*b*u/d + c*v/d - b*w
        tmp3[3] = a*u + b*v + c*w #a*c*u/d + b*v/d + c*w
        =#

        # Some smart guessing method here, to determine whether we actually
        # need to calculate the sounding locations of the entire swath - or
        # whether we can just move on to the next frame.

        # If the distance between ROI center and subsattelite point is larger than twice the
        # swath extremal point distance, skip this!

        swath_end_lons = zeros(2)
        swath_end_lats = zeros(2)

        for i in 1:2

            if i == 1
                idx = 1
            else
                idx = length(theta_array)
            end

            Rz[1,1] = cos_theta_array[idx]
            Rz[2,2] = cos_theta_array[idx]
            Rz[1,2] = -sin_theta_array[idx]
            Rz[2,1] = sin_theta_array[idx]

            #tmp4 = Rz * tmp3
            tmp4[1] = Rz[1,1] * tmp3[1] + Rz[1,2] * tmp3[2]
            tmp4[2] = Rz[2,1] * tmp3[1] + Rz[2,2] * tmp3[2]
            tmp4[3] = tmp3[3]

            look_vector[1] = tmp1[1,1] * tmp4[1] + tmp1[1,2] * tmp4[2] + tmp1[1,3] * tmp4[3]
            look_vector[2] = tmp1[2,1] * tmp4[1] + tmp1[2,2] * tmp4[2] + tmp1[2,3] * tmp4[3]
            look_vector[3] = tmp1[3,1] * tmp4[1] + tmp1[3,2] * tmp4[2] + tmp1[3,3] * tmp4[3]

            # Intersect look vector with WGS84 geoid
            intersect_with_wgs84!(r_itrf, look_vector, r_intersection)
            this_lat, this_lon, this_altitude = ecef_to_geodetic(r_intersection)

            swath_end_lons[i] = rad2deg(this_lon)
            swath_end_lats[i] = rad2deg(this_lat)

        end

        # Distance between first and last footprint in a swath
        distance_swath = calculate_distance(swath_end_lons[1], swath_end_lats[1],
                                            swath_end_lons[2], swath_end_lats[2])

        if ~pasted_swath_width
            println("Swath width: $(distance_swath)")
            pasted_swath_width = true
        end

        # Distance between subsatellite point and target ROI center
        distance_ss_to_ROI_center = calculate_distance(
            sat_lon, sat_lat,
            target_lon, target_lat
        )

        if distance_ss_to_ROI_center - radius > distance_swath
            continue
        end

        # Process the full swath
        for i in 1:length(theta_array)

            Rz[1,1] = cos_theta_array[i]
            Rz[2,2] = cos_theta_array[i]
            Rz[1,2] = -sin_theta_array[i]
            Rz[2,1] = sin_theta_array[i]

            #tmp4 = Rz * tmp3
            tmp4[1] = Rz[1,1] * tmp3[1] + Rz[1,2] * tmp3[2]
            tmp4[2] = Rz[2,1] * tmp3[1] + Rz[2,2] * tmp3[2]
            tmp4[3] = tmp3[3]

            look_vector[1] = tmp1[1,1] * tmp4[1] + tmp1[1,2] * tmp4[2] + tmp1[1,3] * tmp4[3]
            look_vector[2] = tmp1[2,1] * tmp4[1] + tmp1[2,2] * tmp4[2] + tmp1[2,3] * tmp4[3]
            look_vector[3] = tmp1[3,1] * tmp4[1] + tmp1[3,2] * tmp4[2] + tmp1[3,3] * tmp4[3]

            # @views look_vector[:] = Rx' * Ry' * Rz * Ry * Rx * r_sat_to_origin
            # @views look_vector[:] = tmp1 * tmp4

            # Intersect look vector with WGS84 geoid
            intersect_with_wgs84!(r_itrf, look_vector, r_intersection)
            this_lat, this_lon, this_altitude = ecef_to_geodetic(r_intersection)

            this_loc = Geolocation(rad2deg(this_lon), rad2deg(this_lat))

            # This footprint is indeed within our requested target
            if check_location_within_radius(target_lon, target_lat, radius, this_loc)

                this_loctime = GeolocationTime(this_loc, current_time)

                # solar azimuth (saa) at this point unused!
                this_sza, this_saa = calculate_solar_angles(this_loctime)

                # If we are looking at nighttime - skip
                if isnan(this_sza) | (this_sza > 90.0)
                    continue
                end


                # Calculate viewing zenith angles
                # This should be the angle between the vector pointing
                # from measurement location (r_intersections) and the
                # satellite position (r_itrf) in ECEF, and the measurement location
                # and the normal (going from earth center to measurement location)

                r_norm = normalize(r_intersection)
                r_loc_to_sat_norm = normalize(r_itrf - r_intersection)

                _tmp = dot(r_loc_to_sat_norm, r_norm)
                if (_tmp > 1) & (_tmp < 1 + 1e-6)
                    _tmp = 1.0
                end

                this_vza = rad2deg(acos(_tmp))

                this_scene = create_SIF_scene(
                    instrument,
                    "nadir",
                    this_loctime,
                    this_vza,
                    uncertainty_function,
                    vnp_sd
                )

                push!(locarray, this_loc)
                push!(scenearray, this_scene)

            end

        end

        # Compute the new look vectors going from center of spacecraft
        # in ECEF frame.
        #new_look_vectors = rotate_vec_arbitrary.(
        #    Ref(r_sat_to_origin), # point to be rotated
        #    Ref(r_itrf), # position of rotation origin
        #    Ref(v_norm), # rotation axis relative to origin
        #    theta_array) # rotation angle

    end

    return NadirSwathOrbiterSampling(
        "",
        ["instrument"],
        locarray,
        scenearray
    )


end
