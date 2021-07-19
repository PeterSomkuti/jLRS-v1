vnp_year = 2020
vnp_fname = "/home/psomkuti/VNP43C1_2020_full_compressed.h5"

function create_VNP_SD_from_locbounds(lon_min::Real, lat_min::Real,
                                      lon_max::Real, lat_max::Real;
                                      vnp_fname::String=vnp_fname)

    vnp_h5 = h5open(vnp_fname, "r")

    vnp_lons = collect(-180:0.05:180)
    vnp_lats = collect(-60:0.05:90)
    vnp_variables = vnp_h5["names"][:]
    vnp_fill_values = vnp_h5["fill_values"][:]
    vnp_year = 2020
    vnp_times = DateTime(vnp_year, 1, 1) .+ (p -> Dates.Day(min(p, 366))).(0:size(vnp_h5["VNP43C1"])[3] - 1)

    # Trim down to our Antarctica limit
    if lat_min < -60
        lat_min = -59.99999
    end


    idx_lon_min = searchsortedfirst(vnp_lons, lon_min) - 1
    idx_lon_max = searchsortedfirst(vnp_lons, lon_max) - 1

    idx_lat_min = 3001 - searchsortedfirst(vnp_lats, lat_max) + 1
    idx_lat_max = 3001 - searchsortedfirst(vnp_lats, lat_min) + 1

    #println(idx_lon_min, ": ", vnp_lons[idx_lon_min])
    #println(idx_lon_max, ": ", vnp_lons[idx_lon_max])
    #println(idx_lat_min, ": ", vnp_lats[end:-1:1][idx_lat_max])
    #println(idx_lat_max, ": ", vnp_lats[end:-1:1][idx_lat_min])

    # Load data and flip along lat dimension (original data is in decreasing lats)
    data = vnp_h5["VNP43C1"][idx_lat_min:idx_lat_max, idx_lon_min:idx_lon_max,:,:][end:-1:1,:,:,:]

    close(vnp_h5)

    return VNPData(
        data,
        vnp_lons[idx_lon_min:idx_lon_max+1],
        vnp_lats[end:-1:1][idx_lat_min:idx_lat_max+1][end:-1:1],
        vnp_times,
        vnp_year,
        vnp_variables,
        vnp_fill_values
    )

end


function sample_vnp_data(loctime::GeolocationTime, var_list, vnp_sd::VNPData)

    if !(typeof(var_list) <: Vector)
         var_list = [var_list]
    end 

    for var in var_list
        if !(var in vnp_sd.names)
            @error "Sorry! Variable \"$(var)\" not preseng in VNPData object."
            return nothing
        end
        @debug "Variable check OK - all supplied names are found in VNP SD."
    end

    if ((loctime.loc.lon < vnp_sd.lon_bounds[1]) |
        (loctime.loc.lon > vnp_sd.lon_bounds[end]))
        @error "Sorry! Location $(loctime) is outside of longitude bounds!"
        return nothing
    end

    if ((loctime.loc.lat < vnp_sd.lat_bounds[1]) |
        (loctime.loc.lat > vnp_sd.lat_bounds[end]))
        @error "Sorry! Location $(loctime) is outside of latitude bounds!"
        return nothing
    end



    # Find spatial indices
    idx_x = searchsortedfirst(vnp_sd.lon_bounds, loctime.loc.lon) - 1
    idx_y = searchsortedfirst(vnp_sd.lat_bounds, loctime.loc.lat) - 1

    #idx_xplus = min(idx_x + 1, length(SD.lon_bounds))
    #idx_xminus = max(idx_x - 1, 1)
    #x_fac = (loctime.loc.lon - SD.lon_bounds[idx_x]) / (SD.lon_bounds[idx_xplus] - SD.lon_bounds[idx_x])

    #idx_y = searchsortedfirst(SD.lat_bounds, loctime.loc.lat) - 1
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
    # (just replacing the year according to the dataset year)
    newtime = DateTime(
        vnp_sd.data_year,
        Dates.month(loctime.time),
        Dates.day(loctime.time),
        Dates.hour(loctime.time),
        Dates.minute(loctime.time),
        Dates.second(loctime.time)
    )

    idx_t = searchsortedfirst(vnp_sd.time_bounds, newtime) - 1

    idx_var = searchsortedfirst.(Ref(vnp_sd.names), var_list)

    return vnp_sd.data[idx_y, idx_x, idx_t, idx_var]

    #idx_tplus = min(idx_t + 1, length(SD.time_bounds) - 1)
    #t_fac = (newtime - SD.time_bounds[idx_t]) / (SD.time_bounds[idx_tplus] - SD.time_bounds[idx_t])

    #println(idx_tplus, " / ", length(SD.time_bounds))
    #println(loctime)
    #println(idx_t, ", ", idx_tplus, ", factor: ", t_fac)
    #println(idx_xminus, ", ", idx_x, ", ", idx_xplus, ", factor: ", x_fac)
    #println(idx_yminus, ", ", idx_y, ", ", idx_yplus, ", factor: ", y_fac)

    #sampled = SD.data[idx_t][idx_xminus:idx_xplus, idx_yminus:idx_yplus]

    #output = tuple(
    #    [round(sum(((p -> p[i]).(sampled)) .* sp_w)) for i in 1:SD.length]...
    #)

end

function calculate_reflectance(loctime::GeolocationTime, sza::Real, band::String, vnp_sd::VNPData)

    # Stitch together variable names
    var_list = ["BRDF_Albedo_Parameter$(i)_$(band)" for i in 1:3]

    g0_iso = 1.0
    g1_iso = 0.0
    g2_iso = 0.0

    g0_vol = -0.007574
    g1_vol = -0.070987
    g2_vol =  0.307588

    g0_geo = -1.284909
    g1_geo = -0.166314
    g2_geo =  0.041840

    f_iso, f_vol, f_geo = sample_vnp_data(loctime, var_list, vnp_sd)

    # Replace fill values by NaNs

    if f_iso == 32767
        f_iso = NaN
    end

    if f_geo == 32767
        f_geo = NaN
    end

    if f_vol == 32767
        f_vol = NaN
    end

    # Compute black-sky albedo
    sza2 = deg2rad(sza)
    sza3 = deg2rad(sza)

    # 1 / 1000.0 factor is needed to convert Int-valued VIIRS data into
    # weights.

    # White sky albedo
    #return (f_iso * 1.0 + f_vol * 0.189184 + f_geo * (-1.377622)) * 0.001

    # Black sky albedo
    return f_iso * 0.001 *  (g0_iso + g1_iso * sza2 + g2_iso * sza3) +
        f_vol * 0.001 * (g0_vol + g1_vol * sza2 + g2_vol * sza3) +
        f_geo * 0.001 * (g0_geo + g1_geo * sza2 + g2_geo * sza3)
end


function calculate_NDVI(loctime::GeolocationTime, sza::Real)

    R_nir, R_vis = calculate_reflectance(loctime, sza, ["M7", "M5"])

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
