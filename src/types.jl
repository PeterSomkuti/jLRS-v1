# Standard library
using Dates


"""
    Geolocation (lon, lat)

Type to hold a lon/lat pair and internally converts
to -180 to 180 lon range. Checks if supplied coordinates
are within bounds [-180 180] and [-90 90] and throws
an error if not.

# Arguments
- 'lon::Number': Longitude of location
- 'lat::Number': Latitude of location

# Returns
- 'Geolocation': Geolocation object
"""
struct Geolocation

    lon::Number
    lat::Number

    function Geolocation(lon::Number, lat::Number)

        # Convert to -180 to 180
        lon = mod((lon + 180), 360) - 180

        # Check for range
        if (-180 <= lon <= 180) & (-90 <= lat <= 90)
            new(lon, lat)
        else
            error("Location out of range [-180 180], [-90 90].")
        end

    end
end

"""
    GeolocationTime (location, time)

Type to hold both a location {Geolocation} and a {DateTime}
time which holds the measurement time.
"""
struct GeolocationTime
    location::Geolocation
    time::DateTime
end


"""
    Scene

This type represents a measurement, has a location and time, as
well as the associated solar and viewing geometries, the surface
SIF value and a NIRv value.

"""
struct Scene
    # Location and time of scene
    loctime::GeolocationTime
    # SIF value
    SIF::Number
    # SIF single-sounding uncertainty
    SIF_ucert::Number
    # solar zenith angle
    SZA::Number
    # viewing zenith angle
    VZA::Number
    # NIRv
    NIRv::Number
    # Surface albedo (at the SIF wavelength)
    albedo::Number
end


"""
    Some design thoughts on the various sampling types

The main goal of this tool is to provide some OSSE-type insights into
how light response curves for some region of interest are sampled given
some either existing or fictional space-based instrument. Since we
want to be able to "play" with various configurations, we are keeping
the various data of a *sampling* somewhat redundant.

For example, scene locations are separately stored in several places. This
allows to separate the sampling times and repetitions from the location
data, and e.g. allows us to use a function on a sampling object which
adds a repetition.


"""
abstract type Sampling end




# GeoCarb-type sampling for a full day
struct GeostationaryFullDaySampling <: Sampling

end

# Think "repeated Granule"
struct GeostationaryIntensiveSampling <: Sampling

end

# OCO-2/3 type sampling pattern (this uses real data)
struct OCOSampling <: Sampling
    # Scene locations obtained from files
    locations::Array{Geolocation, 1}

    # This is merely an array of scenes which are
    # obtained through files. Can be of different
    # dimensions than scenelocs
    scenes::Array{Scene, 1}
end
