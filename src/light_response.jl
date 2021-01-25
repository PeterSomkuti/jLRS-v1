
"""
    compute_2d_boundaries_indices

Computes the boundaries and indices for a regular spatial sampling,
think e.g. 2x2 degree grid boxes.

"""

function compute_2d_boundaries_indices(lons::Vector{<:Real},
                                       lats::Array{<:Real},
                                       SS::RegularGridCells)

    return calculate_regular_2d_histogram(
        lons, lats,
        SS.delta_lon,
        SS.delta_lat
    )
end


"""
    compute_2d_boundaries_indices

Computes the boundaries and indices for a full ROI, i.e. no further
aggregation and only one (maybe big) grid cell and all scenes fall
into the same [1] grid cell.

"""
function compute_2d_boundaries_indices(lons::Vector{<:Real},
                                       lats::Vector{<:Real},
                                       SS::FullROI)

    lon_grid = [minimum(lons), maximum(lons)]
    lat_grid = [minimum(lats), maximum(lats)]
    idx = repeat([1], length(lons))

    return (lon_grid, lat_grid), idx
end




"""
    aggregate_scenes(IS, SS, TS)

This function will aggregate an array of scenes according to a
regular time sampler. Regular time sampling means "stick all scenes
in a given regular time interval into a separate aggregation container".
Classic example: collect in daily/weekly/monthly intervals.

# Arguments
- 'IS::InstrumentSampling': An instrument sampling object (contains the scenes)
- 'SS::SpatialSampling': A spatial sampling object
- 'TS::RegularTemporalSampling': Regular temporal sampling object

# Returns
- 'Array{Aggregate, 1}': aggregate objects

"""
function aggregate_scenes(IS::InstrumentSampling,
                          SS::SpatialSampling,
                          TS::RegularTemporalSampling; minimum_number=1)

    # For empty instrument samplings, just return an empty container
    if length(IS) == 0
        return Aggregate[]
    end

    if minimum_number < 1
        @error "Minimum number must be >= 1 (you wanted $(minimum_number))"
        return Aggregate[]
    end

    # Just in case the scenes are not ordered in time, established ordering index
    # Perform the temporal aggregation according to TS
    time_idx = sortperm(get_time(IS))

    # Extract first and last times
    start_date = get_time(IS)[time_idx[1]]
    end_date = get_time(IS)[time_idx[end]]

    # This is the time boundary grid
    time_boundary = collect(start_date:Dates.Second(TS.period):end_date)

    # Compute spatial aggregation here
    spatial_bound, spatial_idx = compute_2d_boundaries_indices(
        get_lon(IS),
        get_lat(IS),
        SS
    )

    # This is an empty list of Aggregate objects, which
    # we fill up in the next section of the code.
    aggregates = Aggregate[]

    # General idea:
    # -------------
    # Go through every unique spatial index, and perform a time aggregate
    # on the subset of all scenes which fall into this spatial index. We
    # then end up with a spatio-temporal aggregation according to SS and TS


    for ss_idx in unique(spatial_idx)

        # Remember - spatial index can be of any arbitrary dimension, but
        # will most likely be 1- or 2-dim tuple

        # Run through all time bins
        last_start_index = 1
        for i in 1:length(time_boundary) - 1

            agg_idx = Int[]

            for j in last_start_index:length(IS)

                if (IS.scenes[time_idx[j]].loctime.time >= time_boundary[i]) &
                    (IS.scenes[time_idx[j]].loctime.time < time_boundary[i+1]) &
                    (ss_idx == spatial_idx[time_idx[j]])

                    # time_idx[j] belongs to aggregate bin bounded by
                    # time_boundary[i] and time_boundary[i+1]

                    # Thus, scene with index "j" belongs to spatial subset "ss_idx"
                    # and temporal bin between time_boundary[i] and time_boundary[i+1]

                    # Save current index, so that next iteration of the scene loop
                    # can start from here
                    last_start_index = j

                    # Push this index into index list
                    push!(agg_idx, time_idx[j])

                end
            end

            # Continue to next bin if below threshold
            if length(agg_idx) < minimum_number
                continue
            end

            # Otherwise, create a new aggregate and
            # push it into list
            this_aggregate = Aggregate(
                length(agg_idx),
                agg_idx,
                IS.scenes[agg_idx],
                time_boundary[i],
                time_boundary[i+1],
            )

            push!(aggregates, this_aggregate)

        end # This ends the temporal index loop
    end # This ends the spatial index loop

    return aggregates

end

