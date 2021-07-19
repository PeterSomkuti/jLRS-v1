# Standard library
using Printf
using Dates
import Base: length, +, show, zero


abstract type AbstractLocation end

"""
    Geolocation (lon, lat)

Type to hold a lon/lat pair and internally converts
to -180 to 180 lon range. Checks if supplied coordinates
are within bounds [-180 180] and [-90 90] and throws
an error if not.

# Arguments
- 'lon<:Real': Longitude of location
- 'lat<:Real': Latitude of location

# Returns
- 'Geolocation': Geolocation object
"""
struct Geolocation <: AbstractLocation

    lon::Real
    lat::Real

    function Geolocation(lon, lat)

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
struct GeolocationTime <: AbstractLocation
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
    # What instrument?
    instrument::String
    # Potential sampling mode (nadir, glint, SAM, ..)
    mode::String
    # Location and time of scene
    loctime::GeolocationTime
    # SIF value
    SIF::Real
    # SIF single-sounding uncertainty
    SIF_ucert::Real
    # solar zenith angle
    SZA::Real
    # viewing zenith angle
    VZA::Real
    # NIR radiance
    NIR::Real
    # VIS radiance
    VIS::Real
    # NDVI
    NDVI::Real
    # NIRv (reflectance)
    NIRv::Real
    # NIRv (radiance)
    NIRv_radiance::Real
    # BOA irradiance in SIF units W/m2/nm/sr
    irradiance::Real
    # BOA PPFD
    PPFD::Real
    # Surface reflectance (ideally computed from SZA, VZA via BRDF)
    # such that L_out = reflectance * L_in  * 1000 / pi
    reflectance::Real
    # Total column optical depth (e.g. clouds, aerosols)
    OD::Real
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
    # Information about how this object was created
    info::String

    # Which instrument(s)?
    instrument::Vector{String}

    # Scene locations obtained from files
    locations::Vector{Geolocation}

    # This is merely an array of scenes which are
    # obtained through files. Can be of different
    # dimensions than scenelocs
    scenes::Vector{Scene}
end

# OCO-2/3 type sampling pattern (this uses real data)
struct OCOSampling<:InstrumentSampling
    # Information about how this object was created
    info::String

    # Which instrument(s)?
    instrument::Vector{String}

    # Scene locations obtained from files
    locations::Vector{Geolocation}

    # This is merely an array of scenes which are
    # obtained through files. Can be of different
    # dimensions than scenelocs
    scenes::Vector{Scene}
end

# Orbiter with a grating-type nadir-looking swath instrument (e.g. TROPOMI, CO2M)
struct NadirSwathOrbiterSampling<:InstrumentSampling
    # Information about how this object was created
    info::String

    # Which instrument(s)?
    instrument::Vector{String}

    # Scene locations obtained from files
    locations::Vector{Geolocation}

    # This is merely an array of scenes which are
    # obtained through files. Can be of different
    # dimensions than scenelocs
    scenes::Vector{Scene}
end


# Some new combined sampling
# (fill this in when you have a better idea of what to do with mixed samplings)
struct CombinedSampling<:InstrumentSampling
    # Information about how this object was created
    info::String

    # Which instrument(s)?
    instruments::Vector{String}

    # Scene locations obtained from files
    locations::Vector{Geolocation}

    # This is merely an array of scenes which are
    # obtained through files. Can be of different
    # dimensions than scenelocs
    scenes::Vector{Scene}
end


############################################################
# Various getter functions to obtain data from scenes etc. #
############################################################

get_loc(S::Scene) = S.loctime.loc
get_loctime(S::Scene) = S.loctime
get_lon(S::Scene) = S.loctime.loc.lon
get_lat(S::Scene) = S.loctime.loc.lat
get_time(S::Scene) = S.loctime.time
get_mode(S::Scene) = S.mode
get_sza(S::Scene) = S.SZA
get_mu0(S::Scene) = cos(deg2rad(S.SZA))
get_vza(S::Scene) = S.VZA
get_nir(S::Scene) = S.NIR
get_vis(S::Scene) = S.VIS
get_sif(S::Scene) = S.SIF
get_ppfd(S::Scene) = S.PPFD
get_reflectance(S::Scene) = S.reflectance
get_irradiance(S::Scene) = S.irradiance
get_ndvi(S::Scene) = S.NDVI
get_nirv(S::Scene) = S.NIRv
get_nirv_radiance(S::Scene) = S.NIRv_radiance
get_instrument(S::Scene) = S.instrument
get_sif_ucert(S::Scene) = S.SIF_ucert

get_instrument(IS::OCOSampling) = IS.instrument
get_instrument(IS::GeostationaryIntensiveSampling) = IS.instrument

get_loc(IS::InstrumentSampling) = IS.locations
get_loc_lons(IS::InstrumentSampling) = (p -> p.lon).(IS.locations)
get_loc_lats(IS::InstrumentSampling) = (p -> p.lat).(IS.locations)

get_time(IS::InstrumentSampling) = get_time.(IS.scenes)

get_lon(IS::InstrumentSampling) = get_lon.(IS.scenes)
get_lat(IS::InstrumentSampling) = get_lat.(IS.scenes)
get_sza(IS::InstrumentSampling) = get_sza.(IS.scenes)
get_vza(IS::InstrumentSampling) = get_vza.(IS.scenes)
get_nir(IS::InstrumentSampling) = get_nir.(IS.scenes)
get_vis(IS::InstrumentSampling) = get_vis.(IS.scenes)
get_sif(IS::InstrumentSampling) = get_sif.(IS.scenes)
get_mode(IS::InstrumentSampling) = get_mode.(IS.scenes)
get_ppfd(IS::InstrumentSampling) = get_ppfd.(IS.scenes)
get_reflectance(IS::InstrumentSampling) = get_reflectance.(IS.scenes)
get_irradiance(IS::InstrumentSampling) = get_irradiance.(IS.scenes)
get_ndvi(IS::InstrumentSampling) = get_ndvi.(IS.scenes)
get_nirv(IS::InstrumentSampling) = get_nirv.(IS.scenes)
get_nirv_radiance(IS::InstrumentSampling) = get_nirv_radiance.(IS.scenes)
get_sif_ucert(IS::InstrumentSampling) = get_sif_ucert.(IS.scenes)



# Define some pretty printing for our types

function show(io::IO, S::Scene)
    print(io, [get_lon(S), get_lat(S)])
end

function show(io::IO, ::MIME"text/plain", S::Scene)
    println(io, @sprintf("%s scene at [%0.4f, %0.4f]", get_instrument(S), get_lon(S), get_lat(S)))
    println(io, Dates.format(get_time(S), "yyyy u dd, HH:MM:SS"))
    println(io, @sprintf("VZA: %0.3f deg", get_vza(S)))
    println(io, @sprintf("SZA: %0.3f deg", get_sza(S)))
    println(io, @sprintf("Reflectance@755: %0.3f", get_reflectance(S)))
    println(io, @sprintf("NDVI: %0.3f", get_ndvi(S)))
    println(io, @sprintf("NIRv: %0.3f", get_nirv(S)))
end


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

    IS1_times = get_time.(IS1.scenes)
    IS2_times = get_time.(IS2.scenes)

    intersect_times = intersect(IS1_times, IS2_times)
    srt_IS1 = searchsortedfirst.(Ref(IS1_times), intersect_times)
    srt_IS2 = searchsortedfirst.(Ref(IS2_times), intersect_times)

    exclude_list = Int[]

    # Loop through scenes, determing if location+time
    # and viewing zenith are the same. If so, push them into
    # exclusion list
    for i in 1:length(srt_IS1)

        s1 = IS1.scenes[srt_IS1[i]]
        s2 = IS2.scenes[srt_IS2[i]]

        check::Bool = (
            (get_loctime(s1) == get_loctime(s2)) &
            (get_vza(s1) == get_vza(s2))
        )

        if check
            push!(exclude_list, srt_IS2[i])
        end
    end

    # From the exclude list, build an include list via
    # set operations.
    include_list = setdiff(
        range(1, length=length(IS2)),
        exclude_list
    )

    # Depending on what type combinations we have, create a new object
    if T1 == T2

        # If both InstrumentSamplings are of the same type
        # (i.e. from the satellite kind), we create a new
        # object of the same type.

        if IS1_instrument == IS2_instrument
            # If the two instrument descriptors are the same,
            # just use the same instrument descriptor
            newinstrument = IS1_instrument
        else
            # Otherwise, build a new array with unique entries only.
            # For example ["oco2"] and ["oco3"] will become
            # ["oco2", "oco3"].
            newinstrument = sort(unique(
                vcat(IS1_instrument, IS2_instrument)
            ))
        end

        # New scenes will be also sorted by time, even though it will
        # shuffle instruments.
        newscenes = vcat(IS1.scenes, IS2.scenes[include_list])
        newsort = sortperm(get_time.(newscenes))

        newsampling = T1(
            newinfo,
            newinstrument,
            newlocations,
            newscenes[newsort] # Pass time-sorted scenes
        )

        # Last check if scenes are truly sorted by time
        if issorted(get_time(newsampling))
            return newsampling
        else
            @error "Resulting InstrumentSampling is not time-sorted. Debug."
            return nothing
        end

    else

        @error "Not implemented yet!"
        return nothing

    end

end


abstract type TemporalSampling end

# Some regular temporal sampling, think "X days/weeks/months"
struct RegularTemporalSampling<:TemporalSampling
    # Time in [s] between two aggregation boundaries
    period::Real
end

function HourlySampling(interval::Number)
    # There are 60 * 60 = 3600 seconds in an hour
    return RegularTemporalSampling(3600 * interval)
end

function WeeklySampling(interval::Number)
    # There are 60 * 60 * 24 * 7 = 604800 seconds in a week
    return RegularTemporalSampling(604800 * interval)
end

function DailySampling(interval::Number)
    # There are 60 * 60 * 24 = 86400 seconds in a day
    return RegularTemporalSampling(86400 * interval)
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
    delta_lon::Real
    delta_lat::Real
end


"""
    Aggregate

Represents a SIF aggregate
"""
struct Aggregate
    # Number of measurements in aggregate
    N::Int
    # Scene indices, in case you want to refer back to other data
    scene_idx::Vector{Int}
    # Scenes
    scenes::Vector{Scene}
    # Start time
    start_time::DateTime
    # End time
    end_time::DateTime
end


# Struct for holding surface data which is then accessed
# by sampling functions

abstract type SurfaceData end

struct VNPData <: SurfaceData
    data::Array{Int, 4} ## Lat, Lon, DOY, variable
    lon_bounds::Vector{Real}
    lat_bounds::Vector{Real}
    time_bounds::Vector{DateTime}
    data_year::Int
    # How many entries per coordinate?
    names::Vector{String}
    fill_values::Vector{Real}
end


# Struct to hold SIF data given some lon/lat boundaries
# that can then be sampled by the sampling functions

abstract type SIFData end

struct TROPOMISIFData <: SIFData
    data::Array{Float32, 3} ### Lat, Lon, Month
    lon_bounds::Vector{Real}
    lat_bounds::Vector{Real}
    time_bounds::Vector{DateTime}
    fill_values::Real
end
