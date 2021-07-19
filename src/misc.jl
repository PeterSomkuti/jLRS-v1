"""
    calculate_distance (lon1, lat1, lon2, lat2)

Returns the distance between coordinates (lon1, lat1) and (lon2, lat2) in
units of [km].

# Arguments
- 'lon1::Number': Longitude of point 1
- 'lat1::Number': Latitude of point 1
- 'lon2::Number': Longitude of point 2
- 'lat2::Number': Latitude of point 2

# Returns
- 'Number': Distance between (lon1, lat1) and (lon2, lat2) in kilometers.
"""
function calculate_distance(lon1, lat1, lon2, lat2)

    # Perform in-range checks
    if (lon1 < -180) | (lon1 > 180)
        throw(DomainError(lon1, "lon1 must be [-180 180]"))
    end

    if (lon2 < -180) | (lon2 > 180)
        throw(DomainError(lon2, "lon2 must be [-180 180]"))
    end

    if (lat1 < -90) | (lat1 > 90)
        throw(DomainError(lat1, "lat1 must be [-90 90]"))
    end

    if (lat2 < -90) | (lat2 > 90)
        throw(DomainError(lat2, "lat2 must be [-90 90]"))
    end

    deltalon = deg2rad(lon2 - lon1)

    lat1_rad = deg2rad(lat1)
    lat2_rad = deg2rad(lat2)
    deltalat = lat2_rad - lat1_rad

    a = sin(deltalat / 2.0)^2 + cos(lat1_rad) * cos(lat2_rad) * sin(deltalon / 2.0)^2

    return 2.0 * 6378.0 * asin(sqrt(a))

end


"""
    mask_locations_within_radius (lon, lat, radius, locations)

Takes an array of locations {Geolocation} along with a center lon and
center lat coordinate, a radius (in [km]) and returns a Bool array which
represents whether a location is within the radius or not.

# Arguments
- 'lon::Number': Longitude of central coordinate
- 'lat::Number': Latitude of central coordinate
- 'radius::Number': Radius [km] within which points will be kept
- 'locations::Array{Geolocation, 1}': Array of Geolocations to be masked

# Returns
- 'Array{Bool, 1}': Boolean mask array

"""
function check_location_within_radius(lon, lat, radius, location::Geolocation)

    dist = calculate_distance(
            lon, lat,
            location.lon,
            location.lat
        )

    return dist <= radius

end


function check_location_within_radius(lon, lat, radius, loctime::GeolocationTime)

    dist = calculate_distance(
        lon, lat,
        loctime.loc.lon,
        loctime.loc.lat
    )

    return dist <= radius

end


"""


"""
function calculate_solar_angles(loctime::GeolocationTime)

    # Convert to Julian date
    jdate = jdcnv(loctime.time)
    # Obtain solar position
    sun_radec = sunpos(jdate)
    # Get the local position
    loc_radec = zenpos(jdate, loctime.loc.lat, loctime.loc.lon)
    # Compute apparent ALT and AZ of the sun from (lon, lat) position
    sun_alt, sun_az, sun_hour = eq2hor(sun_radec[1], sun_radec[2], jdate, loctime.loc.lat, loctime.loc.lon)

    if sun_alt < 0
        return NaN, sun_az
    else
        # Return SZA and SAA in degrees
        return 90.0 - sun_alt, sun_az
    end

end


"""
    Function to return bin indices in one dimension, taken from:
https://stackoverflow.com/questions/54879412/get-the-mapping-from-each-element-of-input-to-the-bin-of-the-histogram-in-julia

"""
binindices(edges, data) = searchsortedlast.(Ref(edges), data)

"""

Function to produce a 2D histogram which also returns the indices
of the 2D grid in which certain coordinates fall into (stdlib Julia functions
don't return the indices sadly, so we built our own)

"""
function calculate_regular_2d_histogram(lons::Array{<:Number, 1},
                                        lats::Array{<:Number, 1},
                                        delta_lon,
                                        delta_lat)

    # Construct the lon and lat boundaries
    # ------------------------------------

    # The grids are constructed such that the space between first
    # bin edge and the smallest data point is the same as the
    # last bin edge and the largest data point.

    Nlon = convert(Int, ceil((maximum(lons) - minimum(lons)) / delta_lon))
    Nlat = convert(Int, ceil((maximum(lats) - minimum(lats)) / delta_lat))

    lon_bin0 = (minimum(lons) + maximum(lons) - (Nlon) * delta_lon) / 2
    lat_bin0 = (minimum(lats) + maximum(lats) - (Nlat) * delta_lat) / 2

    lon_grid = collect(range(lon_bin0, step=delta_lon, length=Nlon+1))
    lat_grid = collect(range(lat_bin0, step=delta_lat, length=Nlat+1))

    # Assign to bin indices
    lon_idx = binindices(lon_grid, lons)
    lat_idx = binindices(lat_grid, lats)

    cell_idx = collect(zip(lon_idx, lat_idx))

    return (lon_grid, lat_grid), cell_idx

end

function subset_scene_time(IS::InstrumentSampling,
                           date_min::DateTime,
                           date_max::DateTime)

    T = typeof(IS)

    # Understand which scene time fall into the time window
    idx_times = findall(
        (get_time(IS) .>= date_min) .&
        (get_time(IS) .<= date_max)
    )

    # Which of those resulting locations drop out of the locations list?
    newloc = intersect(
        get_loc(IS),
        (p -> p.loctime.loc).(IS.scenes[idx_times])
    )

    return T(
        IS.info,
        IS.instrument,
        newloc,
        IS.scenes[idx_times]
    )

end
