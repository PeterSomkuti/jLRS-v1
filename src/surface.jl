# 1.) Build BRDF data grid from files
# 2.) Define access functions that operate on those grids

# Populate arrays
vnp_lons = collect(-180:0.05:180)
vnp_lats = collect(-90:0.05:90)
vnp_flist = readdir("/data10/psomkuti/VNP43C1/aggregate_18day")
# VNP dataset is every six days in the year 2020 (366 days)
vnp_year = 2020
vnp_chunk = 18
# Create an array of DateTimes corresponding to the time slots according
# to the chunk size (in days). Maxes out at 366 days, so the last chunk
# will always be shorter if chunk size is not divisible by 366.
vnp_timebounds = DateTime(2020, 1, 1) .+ (p -> Dates.Day(min(p, 366))).(0:vnp_chunk:366+vnp_chunk)

vnp_variables = [
    ("BRDF_Albedo_Parameter1_nir", 32767),
    ("BRDF_Albedo_Parameter2_nir", 32767),
    ("BRDF_Albedo_Parameter3_nir", 32767),
    ("BRDF_Albedo_Parameter1_vis", 32767),
    ("BRDF_Albedo_Parameter2_vis", 32767),
    ("BRDF_Albedo_Parameter3_vis", 32767),
    ("BRDF_Albedo_Parameter1_M5M7", 32767),
    ("BRDF_Albedo_Parameter2_M5M7", 32767),
    ("BRDF_Albedo_Parameter3_M5M7", 32767),
    ("Percent_Snow", 255),
]

vnp_assignments = Dict(
    "BRDF_Albedo_NIR" => ["BRDF_Albedo_Parameter1_nir",
                          "BRDF_Albedo_Parameter2_nir",
                          "BRDF_Albedo_Parameter3_nir"],
    "BRDF_Albedo_VIS" => ["BRDF_Albedo_Parameter1_vis",
                          "BRDF_Albedo_Parameter2_vis",
                          "BRDF_Albedo_Parameter3_vis"],
    "BRDF_Albedo_755" => ["BRDF_Albedo_Parameter1_M5M7",
                          "BRDF_Albedo_Parameter2_M5M7",
                          "BRDF_Albedo_Parameter3_M5M7"],
    "Percent_Snow" => ["Percent_Snow"]
)

vnp_fill_values = Dict(
    "BRDF_Albedo_Parameter1_nir" => 32767,
    "BRDF_Albedo_Parameter2_nir" => 32767,
    "BRDF_Albedo_Parameter3_nir" => 32767,
    "BRDF_Albedo_Parameter1_vis" => 32767,
    "BRDF_Albedo_Parameter2_vis" => 32767,
    "BRDF_Albedo_Parameter3_vis" => 32767,
    "BRDF_Albedo_Parameter1_M5M7" => 32767,
    "BRDF_Albedo_Parameter2_M5M7" => 32767,
    "BRDF_Albedo_Parameter3_M5M7" => 32767,
    "Percent_Snow" => 255
)

vnp_sparse_dict = Dict()
for var in keys(vnp_assignments)
    this_N = length(vnp_assignments[var])
    # This dictionary will be of type (Int32 index -> Tuple Int16s)
    vnp_sparse_dict[var]= SparseMatrixCSC{NTuple{this_N, Int16}, Int32}[]
end

@info "Read-IN of VIIRS BRDF data."
for fname in vnp_flist

    h5open("/data10/psomkuti/VNP43C1/aggregate_18day/" * fname, "r") do h5

        @info "Processing $(fname)"
        for vgroup in keys(vnp_assignments)

            x_coords = []
            y_coords = []
            z_coords = []

            for var in vnp_assignments[vgroup]

                # Remember! in latitude dimension, arrays need to be flipped
                # so that lower indices are at lower latitudes
                tmp_data = h5[var][:,:][:,end:-1:1]
                tmp_idx = findall(tmp_data .!= vnp_fill_values[var])

                push!(x_coords, (p -> p[1]).(tmp_idx))
                push!(y_coords, (p -> p[2]).(tmp_idx))
                push!(z_coords, tmp_data[tmp_idx])

            end

            # Results will be all tuples!
            z_out = collect(zip(z_coords...))

            tmp_sparse = sparse(
                convert.(Int32, x_coords[1]),
                convert.(Int32, y_coords[1]),
                z_out,
                7200, 3600
            )

            push!(vnp_sparse_dict[vgroup], tmp_sparse)

        end
    end
end

vnp_surface_data = Dict()

for var in keys(vnp_assignments)

    vnp_surface_data[var] = VNPSparseData(
        vnp_sparse_dict[var],
        vnp_lons,
        vnp_lats,
        vnp_timebounds,
        vnp_year,
        length(vnp_assignments[var])
    )

    vnp_sparse_dict[var] = Nothing

end

vnp_sparse_dict = Nothing

function sample_surface_data(loctime::GeolocationTime, SD::SurfaceData)

    #Smoothing out makes things rather slow!

    # Find spatial indices
    idx_x = searchsortedfirst(SD.lon_bounds, loctime.loc.lon) - 1
    #idx_xplus = min(idx_x + 1, length(SD.lon_bounds))
    #idx_xminus = max(idx_x - 1, 1)
    #x_fac = (loctime.loc.lon - SD.lon_bounds[idx_x]) / (SD.lon_bounds[idx_xplus] - SD.lon_bounds[idx_x])

    idx_y = searchsortedfirst(SD.lat_bounds, loctime.loc.lat) - 1
    #idx_yplus = min(idx_y + 1, length(SD.lat_bounds))
    #idx_yminus = max(idx_y - 1, 1)
    #y_fac = (loctime.loc.lat - SD.lat_bounds[idx_y]) / (SD.lat_bounds[idx_yplus] - SD.lat_bounds[idx_y])

    #center_weight = 0.25
    #w = (1.0 - center_weight) / 2.0
    # Spatial weights using x_fac and y_fac
    #sp_w = [0       (1-y_fac) * w      0;
    #        (1-x_fac) * w    center_weight      (x_fac) * w;
    #        0        (y_fac) * w       0]
    # Create new time based on the year of the surface dataset
    newtime = DateTime(
        (SD.data_year)...,
        Dates.yearmonthday(loctime.time)[2:end]...,
        (Dates.hour(loctime.time))...,
        (Dates.minute(loctime.time))...,
        (Dates.second(loctime.time))...
    )

    idx_t = searchsortedfirst(SD.time_bounds, newtime) - 1
    idx_tplus = min(idx_t + 1, length(SD.time_bounds) - 1)
    t_fac = (newtime - SD.time_bounds[idx_t]) / (SD.time_bounds[idx_tplus] - SD.time_bounds[idx_t])

    #println(idx_tplus, " / ", length(SD.time_bounds))
    #println(loctime)
    #println(idx_t, ", ", idx_tplus, ", factor: ", t_fac)
    #println(idx_xminus, ", ", idx_x, ", ", idx_xplus, ", factor: ", x_fac)
    #println(idx_yminus, ", ", idx_y, ", ", idx_yplus, ", factor: ", y_fac)

    #sampled = SD.data[idx_t][idx_xminus:idx_xplus, idx_yminus:idx_yplus]

    #output = tuple(
    #    [round(sum(((p -> p[i]).(sampled)) .* sp_w)) for i in 1:SD.length]...
    #)

    #return (1 - t_fac) .* SD.data[idx_t][idx_x, idx_y] .+ (t_fac .* SD.data[idx_tplus][idx_x, idx_y])
    return SD.data[idx_t][idx_x, idx_y]
end

function sample_surface_data_old(loctime::GeolocationTime, SD::SurfaceData)

    # Find spatial indices
    idx_x = searchsortedfirst(SD.lon_bounds, loctime.loc.lon) - 1
    idx_y = searchsortedfirst(SD.lat_bounds, loctime.loc.lat) - 1

    # Create new time based on the year of the surface dataset
    newtime = DateTime(
        (SD.data_year)...,
        Dates.yearmonthday(loctime.time)[2:end]...,
        (Dates.hour(loctime.time))...,
        (Dates.minute(loctime.time))...,
        (Dates.second(loctime.time))...
    )

    idx_t = searchsortedfirst(SD.time_bounds, newtime) - 1

    return SD.data[idx_t][idx_x, idx_y]
end


function calculate_reflectance(loctime::GeolocationTime, sza::Real, band::String)

    g0_iso = 1.0
    g1_iso = 0.0
    g2_iso = 0.0

    g0_vol = -0.007574
    g1_vol = -0.070987
    g2_vol =  0.307588

    g0_geo = -1.284909
    g1_geo = -0.166314
    g2_geo =  0.041840

    f_iso, f_vol, f_geo = sample_surface_data(loctime, vnp_surface_data["BRDF_Albedo_$(band)"])

    # Compute black-sky albedo
    sza2 = deg2rad(sza)
    sza3 = deg2rad(sza)

    # 1 / 1000.0 factor is needed to convert Int-valued VIIRS data into
    # weights.

    #return (f_iso * 1.0 + f_vol * 0.189184 + f_geo * (-1.377622)) * 0.001

    return f_iso * 0.001 *  (g0_iso + g1_iso * sza2 + g2_iso * sza3) +
        f_vol * 0.001 * (g0_vol + g1_vol * sza2 + g2_vol * sza3) +
        f_geo * 0.001 * (g0_geo + g1_geo * sza2 + g2_geo * sza3)
end


function calculate_NDVI(loctime::GeolocationTime, sza::Real)

    R_nir = calculate_reflectance(loctime, sza, "NIR")
    R_vis = calculate_reflectance(loctime, sza, "VIS")

    if R_nir == 0.0
        return 0.0
    elseif R_vis == 0.0
        return 0.0
    else
        return (R_nir - R_vis) / (R_nir + R_vis)
    end

end

function calculate_NIRv(loctime::GeolocationTime, sza::Real)

    R_nir = calculate_reflectance(loctime, sza, "NIR")
    NDVI = calculate_NDVI(loctime, sza)

    return R_nir * NDVI
end
