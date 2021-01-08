
"""
    aggregate_scenes(scenes, TS)

This function will aggregate an array of scenes according to a
regular time sampler.

# Arguments
- 'scenes::Array{Scene, 1}': Array of scenes
- 'TS::RegularTemporalSampling': Regular temporal sampling object

# Returns
- 'Array{Aggregate, 1}': aggregate objects

"""
function aggregate_scenes(scenes::Array{Scene, 1}, TS::RegularTemporalSampling)

    # Just in case the scenes are not ordered in time, established ordering index

    # Perform the temporal aggregation according to TS
    time_idx = sortperm([scenes[i].loctime.time for i in 1:length(scenes)])

    # Extract first and last times
    start_date = scenes[time_idx[1]].loctime.time
    end_date = scenes[time_idx[end]].loctime.time

    # This is the time boundary grid
    time_boundary = collect(start_date:Dates.Second(TS.period):end_date)

    aggregates = Aggregate[]

    last_start_index = 1
    for i in 1:length(time_boundary) - 1

        agg_idx = Int[]

        for j in last_start_index:length(scenes)

            if (scenes[time_idx[j]].loctime.time >= time_boundary[i]) &
                (scenes[time_idx[j]].loctime.time < time_boundary[i+1])

                # time_idx[j] belongs to aggregate bin bounded by
                # time_boundary[i] and time_boundary[i+1]

                # Save current index, so that next loop can start from here
                last_start_index = j

                # Push this index into index list
                push!(agg_idx, time_idx[j])

            end

        end

        # Continue to next bin if empty
        if length(agg_idx) == 0
            continue
        end

        # Otherwise, create a new aggregate and
        # push it into list
        this_aggregate = Aggregate(
            length(agg_idx),
            time_boundary[i],
            time_boundary[i+1],
            (p -> p.SIF).(scenes[agg_idx]),
            (p -> p.SIF_ucert).(scenes[agg_idx]),
            (p -> p.SZA).(scenes[agg_idx]),
            (p -> p.VZA).(scenes[agg_idx]),
            (p -> p.NIRv).(scenes[agg_idx]),
            (p -> p.albedo).(scenes[agg_idx])
        )

        push!(aggregates, this_aggregate)

    end

    return aggregates

end


function calculate_light_response(SS::SpatialSampling, TS::TemporalSampling)

    println("$(typeof(SS)), $(typeof(TS))")

    println("Aggregating scenes")
    aggregate_scenes(SS.scenes, TS)





end
