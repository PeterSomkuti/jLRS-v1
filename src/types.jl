# Standard library
using Dates
import Base: length,+

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
    loc::Geolocation
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
    # Scenes
    scenes::Array{Scene, 1}
    # Start time
    start_time::DateTime
    # End time
    end_time::DateTime
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
    # Information about how this object was created
    info::String

    # Which instrument?
    instrument::Array{String, 1}

    # Scene locations obtained from files
    locations::Array{Geolocation, 1}

    # This is merely an array of scenes which are
    # obtained through files. Can be of different
    # dimensions than scenelocs
    scenes::Array{Scene, 1}
end

# Some new combined sampling
# (fill this in when you have a better idea of what to do with mixed samplings)
struct CombinedSampling<:InstrumentSampling
    info::String
    instruments::Array{String, 1}
    # Scene locations obtained from files
    locations::Array{Geolocation, 1}

    # This is merely an array of scenes which are
    # obtained through files. Can be of different
    # dimensions than scenelocs
    scenes::Array{Scene, 1}
end


############################################################
# Various getter functions to obtain data from scenes etc. #
############################################################

get_instrument(IS::OCOSampling) = IS.instrument

get_locations(IS::InstrumentSampling) = IS.locations
get_locations_lons(IS::InstrumentSampling) = (p -> p.lon).(get_locations(IS))
get_locations_lats(IS::InstrumentSampling) = (p -> p.lat).(get_locations(IS))

get_scene_times(IS::InstrumentSampling) = (p -> p.loctime.time).(IS.scenes)
get_scene_locations(IS::InstrumentSampling) = (p -> p.loctime.loc).(IS.scenes)
get_scene_lons(IS::InstrumentSampling) = (p -> p.lon).(get_scene_locations(IS))
get_scene_lats(IS::InstrumentSampling) = (p -> p.lat).(get_scene_locations(IS))

# Define an addition (merger) of two InstrumentSamplings
# (difficult? this should work for all instrument types?)
function +(IS1::T1, IS2::T2) where {T1<:InstrumentSampling, T2<:InstrumentSampling}

    newinfo = "$(IS1.info) + $(IS2.info)"

    # We do not want any duplicates in the merged IS object,
    # so that we don't overcount measurements
    # ------------------------------------------------------

    # Simple set union works just fine, and is fast
    newlocations = union(IS1.locations, IS2.locations)

    # Scenes are a little more tricky:
    # It could happen that the same exact location is sampled by
    # multiple instruments at the same time. So we consider a scene
    # a duplicate if the location, the time, and the VZA is the same.

    # We are therefore performing a negative search - find scenes for
    # with the time is the same (IS.scenes[].loctime.time is time-ordered),
    # and then check for those duplicate scenes whether the other
    # quantities are the same. If they are, those scenes of IS2 are
    # added to the exclusion list which are not taken up into the
    # new merged IS object.

    IS1_instrument = get_instrument(IS1)
    IS2_instrument = get_instrument(IS2)

    IS1_times = get_scene_times(IS1)
    IS2_times = get_scene_times(IS2)

    intersect_times = intersect(IS1_times, IS2_times)
    srt_IS1 = searchsortedfirst.(Ref(IS1_times), intersect_times)
    srt_IS2 = searchsortedfirst.(Ref(IS2_times), intersect_times)

    exclude_list = Int[]

    for i in 1:length(srt_IS1)

        s1 = IS1.scenes[srt_IS1[i]]
        s2 = IS2.scenes[srt_IS2[i]]

        check = (
            (s1.loctime == s2.loctime) &
            (s1.VZA == s2.VZA)
        )

        if check

            push!(exclude_list, srt_IS2[i])

        end
    end

    # From the exclude list, build an include list
    include_list = setdiff(
        collect(1:length(IS2.scenes)),
        exclude_list
    )


    # Depending on what type combinations we have, create a new object
    if T1 == T2

        if IS1_instrument == IS2_instrument
            # If the two instrument descriptors are the same,
            # just use the same instrument descriptor
            newinstrument = IS1_instrument
        else
            # Otherwise, build a new array with unique entries only
            newinstrument = unique(
                vcat(IS1_instrument, IS1_instrument)
            )
        end

        return T1(
            newinfo,
            newinstrument,
            newlocations,
            vcat(IS1.scenes, IS2.scenes[include_list])
        )

    else

        return nothing

    end


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


