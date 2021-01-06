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
function mask_locations_within_radius(lon, lat, radius, locations::Array{Geolocation, 1})

    mask = zeros(Bool, length(locations))

    for i in 1:length(locations)

        dist = calculate_distance(
            lon, lat,
            locations[i].lon,
            locations[i].lat
        )

        mask[i] = dist <= radius
    end

    return mask

end
