# Standard library
using Dates
import Base: length

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
    Aggregate

Represents a SIF aggregate
"""
struct Aggregate
    # Number of measurements in aggregate
    N::Int
    # Scene indices, in case you want to refer back to other data
    scene_idx::Array{Int, 1}
    # Start time
    start_time::DateTime
    # End time
    end_time::DateTime
    # SIF
    SIF::Array{<:Number, 1}
    # SIF single-sounding uncertainty
    SIF_ucert::Array{<:Number, 1}
    # solar zenith angle
    SZA::Array{<:Number, 1}
    # viewing zenith angle
    VZA::Array{<:Number, 1}
    # NIRv
    NIRv::Array{<:Number, 1}
    # Surface albedo (at the SIF wavelength)
    albedo::Array{<:Number, 1}
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
abstract type InstrumentSampling end

"""
    length (IS)

# Arguments
- 'IS::InstrumentSampling': InstrumentSampling object

# Returns
- 'Int': Number of scenes (not scene locations) within the sampling
"""
function length(IS::InstrumentSampling)
    return length(IS.scenes)
end


# GeoCarb-type sampling for a full day
struct GeostationaryFullDaySampling<:InstrumentSampling

end

# Think "repeated Granule"
struct GeostationaryIntensiveSampling<:InstrumentSampling

end

# OCO-2/3 type sampling pattern (this uses real data)
struct OCOSampling<:InstrumentSampling
    info::String
    # Scene locations obtained from files
    locations::Array{Geolocation, 1}

    # This is merely an array of scenes which are
    # obtained through files. Can be of different
    # dimensions than scenelocs
    scenes::Array{Scene, 1}
end



abstract type TemporalSampling end

# Some regular temporal sampling, think "X days/weeks/months"
struct RegularTemporalSampling<:TemporalSampling
    # Time in [s] between two aggregation boundaries
    period::Number
end

function WeeklySampling(interval::Number)
    # There are 60 * 60 * 24 * 7 = 604800 seconds in a week
    return RegularTemporalSampling(604800 * interval)
end


# In case we want to investigate some more complex irregular
# temporal samplings.
# (e.g. correlations with periods of stable precipitation etc.)
struct IrregularTemporalSampling<:TemporalSampling
    time_boundaries::Array{DateTime, 1}
end



abstract type SpatialSampling end

# No further subsampling in spatial dimension
# treat whole ROI of InstrumentSampling object
# as one.
struct FullROI<:SpatialSampling
end

# Regular X by Y grid cells in lon / lat
struct RegularGridCells<:SpatialSampling
    delta_lon::Number
    delta_lat::Number
end


