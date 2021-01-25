using VegaLite

function build_plot_dataframe(agg::T) where {T}

    if T == Vector{Aggregate}
        # all fine
    elseif T == Aggregate
        agg = Aggregate[agg]
    else
        @error "Unsupported type - needs to be Array{Aggregate, 1}"
        return nothing
    end

    # Figure out the total number of entries required here:
    N = sum((p -> p.N).(agg))

    # Return empty DF if aggregate is empty
    if N == 0
        return DataFrame()
    end

    # Not sure if this is the best option:
    # create a dict with some keys and then cast it into
    # a DataFrame object. Pre-allocating a df and filling
    # it later on doesn't seem to be in the spirit of the
    # DataFrames package..

    df = DataFrame(
        index = fill(-1, N),
        instrument = fill("none", N),
        mode = fill("unknown", N),
        SZA = fill(-1.0, N),
        VZA = fill(-1.0, N),
        lon = fill(-999.99, N),
        lat = fill(-999.99, N)
        )

    last_index = 1
    for i in 1:length(agg)
        df[last_index:last_index + agg[i].N - 1, :index] .= i
        df[last_index:last_index + agg[i].N - 1, :instrument] = (p -> p.instrument).(agg[i].scenes)
        df[last_index:last_index + agg[i].N - 1, :mode] = get_mode.(agg[i].scenes)
        df[last_index:last_index + agg[i].N - 1, :SZA] = get_SZA.(agg[i].scenes)
        df[last_index:last_index + agg[i].N - 1, :VZA] = get_VZA.(agg[i].scenes)
        df[last_index:last_index + agg[i].N - 1, :lon] = get_lon.(agg[i].scenes)
        df[last_index:last_index + agg[i].N - 1, :lat] = get_lon.(agg[i].scenes)

        last_index = last_index + agg[i].N
    end

    return df
end


function aggregate_overview_vl(agg::Vector{Aggregate})

    # Construct DataFrame that holds most needed into
    df = build_plot_dataframe(agg)

    # From that, let's construct a grouped dataframe with some aggregate (of aggregates)
    # values ..

    dfg = groupby(df, :index)

    # Containing the number of scenes per instrument mode
    mode_symbols = [Symbol(y) for y in unique(df.mode)]
    df_modes = combine(dfg, ([:mode => (x -> sum(x .== y)) => Symbol(y) for y in unique(df.mode)])...)



    #################
    # Produce plots #
    #################


    # A set of bar plots showing the occurrence of different instrument
    # modes for the array of aggregates
    p1 = df_modes |> @vlplot(repeat=mode_symbols, columns=3) +
        @vlplot(mark={:bar, "tooltip" = "content" => "data"},
                y={field={repeat=:repeat}, bin=false, type=:quantitative},
                x={"index", title="Aggregate #"})
    display(p1)


    # A bar plot showing the number of scenes per aggregate
    p2 = combine(dfg, :index => length => :count) |> @vlplot() + @vlplot(
        mark={:bar, "tooltip" = "content" => "data"},
        x={"index:o", title="Aggregate #"}, y={"count:q", title="Number of scenes"}
    )
    display(p2)


end
