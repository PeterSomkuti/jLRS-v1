# SIF portion

sif_fname = "/home/psomkuti/s5p_sif_2018-2019_gap_filled.h5"
spectra_fname = "/home/psomkuti/jLRS/data/photosystem_spectra.txt"

ps_spectra_wl = Float32[]
ps_spectra_phiI = Float32[]
ps_spectra_phiII = Float32[]
ps_spectra_raw = readlines(spectra_fname)

for line in ps_spectra_raw
    a, b, c = (p -> parse(Float32, p)).(split(line, ","))
    push!(ps_spectra_wl, a)
    push!(ps_spectra_phiI, b)
    push!(ps_spectra_phiII, c)
end


function create_TROPOMI_SIF_from_locbounds(lon_min,
                                           lat_min,
                                           lon_max,
                                           lat_max)

    # SIF data array: month, lon, lat
    sif_h5 = h5open(sif_fname, "r")

    lon_bounds = sif_h5["lon"][:]
    lat_bounds = sif_h5["lat"][:]

    idx_lon_min = searchsortedfirst(lon_bounds, lon_min) - 1
    idx_lon_max = searchsortedfirst(lon_bounds, lon_max)

    idx_lat_min = searchsortedfirst(lat_bounds, lat_min) - 1
    idx_lat_max = searchsortedfirst(lat_bounds, lat_max)

    # Load SIF subset and switch dimensions so they're
    # time, lon, lat ->
    # lat, lon, time
    sif_subset = permutedims(
        sif_h5["sif"][:,idx_lon_min:idx_lon_max,idx_lat_min:idx_lat_max],
        [3,2,1]
    )

    dates = [DateTime(1999, i, 15) for i in 1:12]

    return TROPOMISIFData(
        sif_subset,
        lon_bounds[idx_lon_min:idx_lon_max],
        lat_bounds[idx_lat_min:idx_lat_max],
        dates,
        NaN
    )

end


function calculate_tropomi_sif(loctime::GeolocationTime,
                               sza::Real,
                               sif_data::TROPOMISIFData)

    idx_x = searchsortedfirst(sif_data.lon_bounds, loctime.loc.lon) - 1
    idx_y = searchsortedfirst(sif_data.lat_bounds, loctime.loc.lat) - 1

    newtime = DateTime(
        Dates.year(sif_data.time_bounds[1]),
        Dates.month(loctime.time),
        Dates.day(loctime.time),
        Dates.hour(loctime.time),
        Dates.minute(loctime.time),
        Dates.second(loctime.time)
    )

    idx_t_left = searchsortedfirst(sif_data.time_bounds, newtime) - 1

    if idx_t_left == 0
        idx_t_left = 12
    end

    if idx_t_left == 12
        idx_t_right = 1
    else
        idx_t_right = idx_t_left + 1
    end

    t_fac = (newtime - sif_data.time_bounds[idx_t_left]) /
        (sif_data.time_bounds[idx_t_right] - sif_data.time_bounds[idx_t_left])

    sif_val = (1.0 - t_fac) * sif_data.data[idx_y, idx_x, idx_t_left] +
        t_fac * sif_data.data[idx_y, idx_x, idx_t_right]

    # Keep in mind: this is SIF at 755nm
    return sif_val * cos(deg2rad(sza))


end



function calculate_fAPAR(NDVI::Real)

    linear_fAPAR = 1.35 * NDVI - 0.32

    if linear_fAPAR < 0.0
        linear_fAPAR = 0.0
    end

    if linear_fAPAR > 1.0
        linear_fAPAR = 1.0
    end

    return linear_fAPAR

end



######
###### J. Johnson's SIF model
###### https://github.com/jenjohnson/johnson-berry-2021-pres
######


function photosynthetic_yields(
    Qin, # input PAR, umol PPFD / m2 / s - you want this to be an array
    Tin, # input leaf temp, C
    Cin, # input mesophyll CO2, ubar
    Oin;  # input atmospheric O2, mbar
    Abs=0.85, # this is fAPAR essentially (?), mol PPFD absorbed by the leaf, PPFD incident / mol absorbed
    alpha_opt="static",
    beta = 0.52, # mol PPFD absorbed by PSII, PPFD / mol absorbed
    CB6F = (350.0 / 300.0) / 1e6, # Cyt b6f density, mol sites per m2
    RUB = (100.0 / 3.6) / 1e6, # Rubisco density, mol sites per m2
    Rds = 0.01, # dark respiration scalar
    Kf = 0.05e9, # Rate constant for fluorescence at PSII and PSI, 1/s
    Kd = 0.55e9, # Rate constant for constitutive heat loss at PSII and PSI, 1/s
    Kp1 = 14.5e9, # Rate constant for photochemistry at PSI, 1/s
    Kp2 = 4.5e9,  # Rate constant for photochemistry at PSII, 1/s
    Kn1 = 14.5e9, # Rate constant for regulated heat loss at PSI, 1/s
    Ku2 = 0.0e9,  # Rate constant for exciton sharing PSII, 1/s 
    kq = 300.0, # Cyt b6f kcat for PQH2, mol mol / e- / sites / s
    nl = 0.75,  # ATP per e- in linear flow, ATP/e-
    nc = 1.00,  # ATP per e- in cyclic flow, ATP/e-
    kc = 3.6,   # Rubisco kcat for CO2, mol CO2 mol / sites / s
    ko = 3.6 * 0.27, # Rubisco kcat for O2, mol O2 mol / sites / s
    Kc = 260.0 / 1e6,  # Rubisco Km for CO2, bar
    Ko = 179000.0 / 1e6, # Rubisco Km for O2, bar
    theta1 = 1., # curvature parameter for Aj/Ac transition
    eps1 = 0., # PSI transfer function
    eps2 = 1. # PSII transfer function
)

    @assert (alpha_opt == "static") | (alpha_opt == "dynamic") "Bad parameter!"

    # Inputs
    Q = Qin ./ 1e6 # PAR in mol PPFD / m2 /s
    T = Tin
    C = Cin / 1e6 # mesophyll CO2 partial pressure, bar
    O = Oin / 1e3 # atmospheric O2 partial pressure, O2

    # Calculations

    Vqmax = CB6F * kq           # Maximum Cyt b6f activity, mol e-1 m-2 s-1
    Vcmax = RUB * kc            # Maximum Rubisco activity, mol CO2 m-2 s-1
    Rd = Vcmax * Rds            # Mitochondrial respiration, mol CO2 m-2 s-1
    S = (kc / Kc) * (Ko ./ ko)   # Rubisco specificity for CO2/O2, dimensionless
    gammas = O / (2 * S)        # CO2 compensation point in the absence of Rd, bar
    eta = (1.0 - (nl / nc) + (3.0 + 7.0 * gammas / C) / ((4.0 + 8.0 * gammas / C) * nc)) # PS I/II ETR
    phi1P_max = Kp1 / (Kp1 + Kd + Kf) # Maximum photochemical yield PS I

    if alpha_opt == "static"
        a2 = Abs .* beta
        a1 = Abs .- a2
    elseif alpha_opt == "dynamic"
        solve_xcs(Abs, CB6F, Kd, Kf, Kp2, Ku2, Q, eta, kq, phi1P_max) = 
            @. (-sqrt((Kd + Kf + Ku2) * (CB6F^2 * Kd^3 * kq^2 * phi1P_max^2 + CB6F^2 * Kf^3 * kq^2 * phi1P_max^2 + CB6F^2 * Kd * Kp2^2 * eta^2 * kq^2 + CB6F^2 * Kf * Kp2^2 * eta^2 * kq^2 + CB6F^2 * Kp2^2 * Ku2 * eta^2 * kq^2 + CB6F^2 * Kd * Kf^2 * kq^2 * phi1P_max^2 * 3.0 + CB6F^2 * Kd^2 * Kf * kq^2 * phi1P_max^2 * 3.0 + CB6F^2 * Kd * Kp2^2 * kq^2 * phi1P_max^2 + CB6F^2 * Kd^2 * Kp2 * kq^2 * phi1P_max^2 * 2.0 + CB6F^2 * Kf * Kp2^2 * kq^2 * phi1P_max^2 + CB6F^2 * Kf^2 * Kp2 * kq^2 * phi1P_max^2 * 2.0 + CB6F^2 * Kd^2 * Ku2 * kq^2 * phi1P_max^2 + CB6F^2 * Kf^2 * Ku2 * kq^2 * phi1P_max^2 + CB6F^2 * Kp2^2 * Ku2 * kq^2 * phi1P_max^2 + Abs^2 * Kd * Kp2^2 * Q^2 * eta^2 * phi1P_max^2 + Abs^2 * Kf * Kp2^2 * Q^2 * eta^2 * phi1P_max^2 + Abs^2 * Kp2^2 * Ku2 * Q^2 * eta^2 * phi1P_max^2 + CB6F^2 * Kd * Kf * Kp2 * kq^2 * phi1P_max^2 * 4.0 + CB6F^2 * Kd * Kf * Ku2 * kq^2 * phi1P_max^2 * 2.0 + CB6F^2 * Kd * Kp2 * Ku2 * kq^2 * phi1P_max^2 * 2.0 + CB6F^2 * Kf * Kp2 * Ku2 * kq^2 * phi1P_max^2 * 2.0 + CB6F^2 * Kd * Kp2^2 * eta * kq^2 * phi1P_max * 2.0 + CB6F^2 * Kd^2 * Kp2 * eta * kq^2 * phi1P_max * 2.0 + CB6F^2 * Kf * Kp2^2 * eta * kq^2 * phi1P_max * 2.0 + CB6F^2 * Kf^2 * Kp2 * eta * kq^2 * phi1P_max * 2.0 + CB6F^2 * Kp2^2 * Ku2 * eta * kq^2 * phi1P_max * 2.0 + CB6F^2 * Kd * Kf * Kp2 * eta * kq^2 * phi1P_max * 4.0 + CB6F^2 * Kd * Kp2 * Ku2 * eta * kq^2 * phi1P_max * 2.0 + CB6F^2 * Kf * Kp2 * Ku2 * eta * kq^2 * phi1P_max * 2.0 + Abs * CB6F * Kd * Kp2^2 * Q * eta * kq * phi1P_max^2 * 2.0 + Abs * CB6F * Kd * Kp2^2 * Q * eta^2 * kq * phi1P_max * 2.0 + Abs * CB6F * Kd^2 * Kp2 * Q * eta * kq * phi1P_max^2 * 2.0 + Abs * CB6F * Kf * Kp2^2 * Q * eta * kq * phi1P_max^2 * 2.0 + Abs * CB6F * Kf * Kp2^2 * Q * eta^2 * kq * phi1P_max * 2.0 + Abs * CB6F * Kf^2 * Kp2 * Q * eta * kq * phi1P_max^2 * 2.0 - Abs * CB6F * Kp2^2 * Ku2 * Q * eta * kq * phi1P_max^2 * 2.0 + Abs * CB6F * Kp2^2 * Ku2 * Q * eta^2 * kq * phi1P_max * 2.0 + Abs * CB6F * Kd * Kf * Kp2 * Q * eta * kq * phi1P_max^2 * 4.0 + Abs * CB6F * Kd * Kp2 * Ku2 * Q * eta * kq * phi1P_max^2 * 2.0 + Abs * CB6F * Kf * Kp2 * Ku2 * Q * eta * kq * phi1P_max^2 * 2.0)) + CB6F * Kd^2 * kq * phi1P_max + CB6F * Kf^2 * kq * phi1P_max + Abs * Kd^2 * Q * phi1P_max^2 * 2.0 + Abs * Kf^2 * Q * phi1P_max^2  *  2.0 + CB6F * Kd * Kp2 * eta * kq + CB6F * Kf * Kp2 * eta * kq + CB6F * Kp2 * Ku2 * eta * kq + CB6F * Kd * Kf * kq * phi1P_max * 2.0 + CB6F * Kd * Kp2 * kq * phi1P_max + CB6F * Kf * Kp2 * kq * phi1P_max + CB6F * Kd * Ku2 * kq * phi1P_max + CB6F * Kf * Ku2 * kq * phi1P_max + CB6F * Kp2 * Ku2 * kq * phi1P_max + Abs * Kd * Kf * Q * phi1P_max^2 * 4.0 + Abs * Kd * Kp2 * Q * phi1P_max^2 * 2.0 + Abs * Kf * Kp2 * Q * phi1P_max^2  *  2.0 + Abs * Kd * Ku2 * Q * phi1P_max^2 * 2.0 + Abs * Kf * Ku2 * Q * phi1P_max^2 * 2.0 + Abs * Kd * Kp2 * Q * eta * phi1P_max + Abs * Kf * Kp2 * Q * eta * phi1P_max + Abs * Kp2 * Ku2 * Q * eta * phi1P_max) ./ (Q * phi1P_max * (Kd^2 * phi1P_max + Kf^2 * phi1P_max + Kd * Kp2 * eta + Kf * Kp2 * eta + Kp2 * Ku2 * eta + Kd * Kf * phi1P_max * 2.0 + Kd * Kp2 * phi1P_max + Kf * Kp2 * phi1P_max + Kd * Ku2 * phi1P_max + Kf * Ku2 * phi1P_max) * 2.0)
            #(-sqrt((Kd+Kf+Ku2)*(CB6F^2*Kd^3*kq^2*phi1P_max^2+CB6F^2*Kf^3*kq^2*phi1P_max^2+CB6F^2*Kd*Kp2^2*eta^2*kq^2+CB6F^2*Kf*Kp2^2*eta^2*kq^2+CB6F^2*Kp2^2*Ku2*eta^2*kq^2+CB6F^2*Kd*Kf^2*kq^2*phi1P_max^2*3.0+CB6F^2*Kd^2*Kf*kq^2*phi1P_max^2*3.0+CB6F^2*Kd*Kp2^2*kq^2*phi1P_max^2+CB6F^2*Kd^2*Kp2*kq^2*phi1P_max^2*2.0+CB6F^2*Kf*Kp2^2*kq^2*phi1P_max^2+CB6F^2*Kf^2*Kp2*kq^2*phi1P_max^2*2.0+CB6F^2*Kd^2*Ku2*kq^2*phi1P_max^2+CB6F^2*Kf^2*Ku2*kq^2*phi1P_max^2+CB6F^2*Kp2^2*Ku2*kq^2*phi1P_max^2+Abs^2*Kd*Kp2^2*Q^2*eta^2*phi1P_max^2+Abs^2*Kf*Kp2^2*Q^2*eta^2*phi1P_max^2+Abs^2*Kp2^2*Ku2*Q^2*eta^2*phi1P_max^2+CB6F^2*Kd*Kf*Kp2*kq^2*phi1P_max^2*4.0+CB6F^2*Kd*Kf*Ku2*kq^2*phi1P_max^2*2.0+CB6F^2*Kd*Kp2*Ku2*kq^2*phi1P_max^2*2.0+CB6F^2*Kf*Kp2*Ku2*kq^2*phi1P_max^2*2.0+CB6F^2*Kd*Kp2^2*eta*kq^2*phi1P_max*2.0+CB6F^2*Kd^2*Kp2*eta*kq^2*phi1P_max*2.0+CB6F^2*Kf*Kp2^2*eta*kq^2*phi1P_max*2.0+CB6F^2*Kf^2*Kp2*eta*kq^2*phi1P_max*2.0+CB6F^2*Kp2^2*Ku2*eta*kq^2*phi1P_max*2.0+CB6F^2*Kd*Kf*Kp2*eta*kq^2*phi1P_max*4.0+CB6F^2*Kd*Kp2*Ku2*eta*kq^2*phi1P_max*2.0+CB6F^2*Kf*Kp2*Ku2*eta*kq^2*phi1P_max*2.0+Abs*CB6F*Kd*Kp2^2*Q*eta*kq*phi1P_max^2*2.0+Abs*CB6F*Kd*Kp2^2*Q*eta^2*kq*phi1P_max*2.0+Abs*CB6F*Kd^2*Kp2*Q*eta*kq*phi1P_max^2*2.0+Abs*CB6F*Kf*Kp2^2*Q*eta*kq*phi1P_max^2*2.0+Abs*CB6F*Kf*Kp2^2*Q*eta^2*kq*phi1P_max*2.0+Abs*CB6F*Kf^2*Kp2*Q*eta*kq*phi1P_max^2*2.0-Abs*CB6F*Kp2^2*Ku2*Q*eta*kq*phi1P_max^2*2.0+Abs*CB6F*Kp2^2*Ku2*Q*eta^2*kq*phi1P_max*2.0+Abs*CB6F*Kd*Kf*Kp2*Q*eta*kq*phi1P_max^2*4.0+Abs*CB6F*Kd*Kp2*Ku2*Q*eta*kq*phi1P_max^2*2.0+Abs*CB6F*Kf*Kp2*Ku2*Q*eta*kq*phi1P_max^2*2.0))+CB6F*Kd^2*kq*phi1P_max+CB6F*Kf^2*kq*phi1P_max+Abs*Kd^2*Q*phi1P_max^2*2.0+Abs*Kf^2*Q*phi1P_max^2*2.0+CB6F*Kd*Kp2*eta*kq+CB6F*Kf*Kp2*eta*kq+CB6F*Kp2*Ku2*eta*kq+CB6F*Kd*Kf*kq*phi1P_max*2.0+CB6F*Kd*Kp2*kq*phi1P_max+CB6F*Kf*Kp2*kq*phi1P_max+CB6F*Kd*Ku2*kq*phi1P_max+CB6F*Kf*Ku2*kq*phi1P_max+CB6F*Kp2*Ku2*kq*phi1P_max+Abs*Kd*Kf*Q*phi1P_max^2*4.0+Abs*Kd*Kp2*Q*phi1P_max^2*2.0+Abs*Kf*Kp2*Q*phi1P_max^2*2.0+Abs*Kd*Ku2*Q*phi1P_max^2*2.0+Abs*Kf*Ku2*Q*phi1P_max^2*2.0+Abs*Kd*Kp2*Q*eta*phi1P_max+Abs*Kf*Kp2*Q*eta*phi1P_max+Abs*Kp2*Ku2*Q*eta*phi1P_max)/(Q*phi1P_max*(Kd^2*phi1P_max+Kf^2*phi1P_max+Kd*Kp2*eta+Kf*Kp2*eta+Kp2*Ku2*eta+Kd*Kf*phi1P_max*2.0+Kd*Kp2*phi1P_max+Kf*Kp2*phi1P_max+Kd*Ku2*phi1P_max+Kf*Ku2*phi1P_max)*2.0)

        a2 = solve_xcs(Abs, CB6F, Kd, Kf, Kp2, Ku2, Q, eta, kq, phi1P_max)
        a1 = Abs .- a2
    end

    # Calculate limiting rates for gas-exchange and electron transport

    # Expressions for potential Rubisco-limited rates (_c)
    #   N.B., see Eqns. 32-33
    JP700_j = @. (Q * Vqmax) / (Q + Vqmax / (a1 * phi1P_max))
    JP680_j = @. JP700_j / eta
    Vc_j = @. JP680_j / (4.0 * (1.0 + 2.0 * gammas / C))
    Vo_j = @. Vc_j * 2.0 * gammas / C
    Ag_j = @. Vc_j - Vo_j / 2.0

    # Expressions for potential Rubisco-limited rates (_c)
    Vc_c = @. C * Vcmax / (C + Kc * (1.0 + O / Ko))
    Vo_c = @. Vc_c * 2.0 * gammas / C
    Ag_c = @. Vc_c - Vo_c / 2.0
    JP680_c = @. Ag_c * 4.0 * (1.0 + 2.0 * gammas / C) / (1.0 - gammas / C)
    # This quantity here is a scalar, but needs to be repeated to be of the
    # same size as JP680_c.
    JP700_c = repeat([JP680_c * eta], length(JP700_j))

    # Define anonymous function for quadratic to smooth transitions
    #   N.B., this returns an array with the two roots of the quadratic

    tr(l1, l2, th) = [
        ((l1 + l2) + sqrt((l1 + l2)^2 - 4.0 * th * l1 * l2)) / (2.0 * th),
        ((l1 + l2) - sqrt((l1 + l2)^2 - 4.0 * th * l1 * l2)) / (2.0 * th)
    ]

    # Select minimum PS1 ETR
    JP700_a = minimum.(tr.(JP700_j, JP700_c, Ref(theta1)))
    JP680_a = minimum.(tr.(JP680_j, JP680_c, Ref(theta1)))

    # Select minimum Ag_a
    Ag_a  = (
        minimum.(tr.(Ag_j, Ag_c, Ref(theta1))) .* (C > gammas) +
        maximum.(tr.(Ag_j, Ag_c, Ref(theta1))) .* (C <= gammas)
    )
    An_a = Ag_a .- Rd

    if alpha_opt == "dynamic"
        # Derive a2/a1 at light saturation point and update
        # N.B., this represents dynamic optimization of a2/a1 under 
        # limiting light, and then under saturating light holds a2/a1 at 
        # the values attained at the light saturation point.

        #I = findmin(hcat(JP700_j, JP700_c), dims=2)[2]
        I = (p -> p.I[2]).(findmin(hcat(JP700_j, JP700_c), dims=2)[2])
        which = findall(diff(I, dims=1) .== 1)[end]

        a2_new = repeat([a2[which.I[1]]], length(a2))
        a2_old = copy(a2)

        a2_update = zeros(length(a2))
        a2_update[a2_old .> a2_new] = a2_old[a2_old .> a2_new]
        a2_update[a2_old .<= a2_new] = a2_new[a2_old .<= a2_new]

        #println(a2_old .> a2_new)
        #println(a2_old .<= a2_new)

        a2 = copy(a2_update)
        a1 = Abs .- a2
    end

    CB6F_a = @. JP700_j / kq         # Eqns. 21, 30a, 34
    phi1P_a = @. JP700_a / (Q * a1)  # Eqn. 20
    q1_a = @. phi1P_a / phi1P_max    # Eqn. 19a
    phi2P_a = @. JP680_a / (Q * a2)  # Eqn. 26
    q2_a = @. 1.0 - CB6F_a / CB6F    # Eqns. 28 and 34

    # N.B., rearrange Eqn. 25a to solve for Kn2_a

    Kn2_a = @. (
        sqrt(Kp2^2 * phi2P_a^2 - 2.0 * Kp2^2 * phi2P_a * q2_a +
             Kp2^2 * q2_a.^2 - 4.0 * Kp2 * Ku2 * phi2P_a^2 * q2_a .+
             2.0 * Kp2 * Ku2 * phi2P_a^2 + 2.0 * Kp2 * Ku2 * phi2P_a * q2_a .+
             Ku2.^2 * phi2P_a^2)
        - Kp2 * phi2P_a + Ku2 * phi2P_a + Kp2 * q2_a) / (2.0 * phi2P_a) - Kf - Ku2 - Kd


    # Photosystem II (Eqns. 23a-23e and 25a-25d)
    phi2p_a = @. (q2_a) * Kp2 / (Kp2 + Kn2_a + Kd + Kf + Ku2)
    phi2n_a = @. (q2_a) * Kn2_a / (Kp2 + Kn2_a + Kd + Kf + Ku2) +
        (1 - q2_a) * Kn2_a / (Kn2_a + Kd + Kf + Ku2)
    phi2d_a = @. (q2_a) * Kd / (Kp2 + Kn2_a + Kd + Kf + Ku2) +
        (1 - q2_a) * Kd / (Kn2_a + Kd + Kf + Ku2)
    phi2f_a = @. (q2_a) * Kf / (Kp2 + Kn2_a + Kd + Kf + Ku2) +
        (1 - q2_a) * Kf ./ (Kn2_a + Kd + Kf + Ku2)
    phi2u_a = @. (q2_a) * Ku2 / (Kp2 + Kn2_a + Kd + Kf + Ku2) +
        (1 - q2_a) * Ku2 / (Kn2_a + Kd + Kf + Ku2)

    phi2P_a = @. phi2p_a / (1.0 - phi2u_a)
    phi2N_a = @. phi2n_a / (1.0 - phi2u_a)
    phi2D_a = @. phi2d_a / (1.0 - phi2u_a)
    phi2F_a = @. phi2f_a / (1.0 - phi2u_a)

    # For Photosystem I (Eqns. 19a-19d)
    phi1P_a = @. q1_a * Kp1 / (Kp1 + Kd + Kf)
    phi1N_a = @. (1.0 - q1_a) * Kn1 / (Kn1 + Kd + Kf)
    phi1D_a = @. q1_a * Kd / (Kp1 + Kd + Kf) + (1.0 - q1_a) * Kd / (Kn1 + Kd + Kf)
    phi1F_a = @. q1_a * Kf / (Kp1 + Kd + Kf) + (1.0 - q1_a) * Kf / (Kn1 + Kd + Kf)

    # PAM measured fluorescence levels (Eqns. 38-42)
    #   N.B., hardcoding of a2(1) for dark-adapted value
    ### this means that we MUST supply a wavelength array with a low radiance
    ### as its first value - for the case of the dynamic model!!
    Fm_a = @. a2[1] * Kf / (Kd + Kf) * eps2 + a1[1] * Kf / (Kn1 + Kd + Kf) * eps1
    Fo_a = @. a2[1] * Kf / (Kp2 + Kd + Kf) * eps2 + a1[1] * Kf / (Kp1 + Kd + Kf) * eps1
    Fmp_a = @. a2 * Kf / (Kn2_a + Kd + Kf) * eps2 + a1 * Kf / (Kn1 + Kd + Kf) * eps1
    Fop_a = @. a2 * Kf / (Kp2 + Kn2_a + Kd + Kf) * eps2 + a1 * Kf / (Kp1 + Kd + Kf) * eps1
    Fs_a = @. a2 * phi2F_a * eps2 + a1 * phi1F_a * eps1

    # PAM indices
    PAM1_a = @. 1.0 - Fs_a / Fmp_a # PhiP
    PAM2_a = @. Fs_a * (1.0 / Fmp_a - 1.0 / Fm_a) # PhiN
    PAM3_a = @. Fs_a / Fm_a # PhiD + PhiF
    PAM4_a = @. Q * 0.85 / 2 * PAM1_a # ETR
    PAM5_a = @. (Fmp_a - Fs_a) / (Fmp_a - Fop_a) # qP
    PAM6_a = @. (Fmp_a - Fs_a) * Fop_a / ((Fmp_a - Fop_a) * Fs_a) # qL
    PAM7_a = @. PAM4_a / (1.0 - PAM5_a) # kPuddle
    PAM8_a = @. PAM4_a / (1.0 - PAM6_a) # kLake
    PAM9_a = @. Fm_a / Fmp_a - 1.0 # NPQ

    return PAM1_a, PAM2_a, PAM3_a, PAM4_a, PAM5_a, PAM6_a, PAM7_a, PAM8_a, PAM9_a, Fs_a, An_a, a1 .* phi1F_a, a2 .* phi2F_a

end

function model_fluorescence(
    Q::Number, # We want this to be a number, not a vector
    T,
    C,
    O,
    wavelength, # this should be in nm
    NDVI
)

    fAPAR_chl = calculate_fAPAR(NDVI)

    if fAPAR_chl <= 1e-3
        return 0.0
    end

    fAPAR_abs = 0.85
    fAPAR_int = fAPAR_chl / fAPAR_abs
    fESC = 1.0

    p1, p2 = photosynthetic_yields(
        Q, # Q in umol PPFD m^-2 s^-1, INCIDENT PAR
        T, # leaf temp in C
        C, # CO2 pp in ubar
        O, # O2 pp in mbar
        Abs=fAPAR_abs * fAPAR_chl, # Total leaf absorbance to PAR [mol / mol]
        alpha_opt="static"
    )[[12, 13]]

    if p1 isa Vector
        p1 = p1[1]
    end

    if p2 isa Vector
        p2 = p2[1]
    end


    # Find where we need to sample the spectra
    wl_idx = argmin(abs.(ps_spectra_wl .- wavelength))
    hcNA = 1.0 / 8.359347229111778 # kg m^3 / mol / s2

    # Q: incoming radiance in PPFD mol / m2 / s
    # outgoing radiance at wavelength [PPFD umol/m2/s/nm]
    # (this is radiating out in to the full space)
    # (phiI, phiII is in units of 1/nm, p1 and p2 is in [1])
    SIF_PPFD = fAPAR_chl * Q * (p1 * ps_spectra_phiI[wl_idx] + p2 * ps_spectra_phiII[wl_idx])

    # outgoing radiance at wavelength [PPFD mol/m2/s/nm]
    SIF_PPFD *= 1e-3
    # outgoing radiance at wavelength [W/m2]
    SIF_PPFD *= hcNA * 1e9
    # outgoing radiance at wavelength [W/m2/um]
    SIF_PPFD /= 1e0 * wavelength
    # outgoing radiance at wavelength [W/m2/um/sr]
    SIF_PPFD /= 4 * pi
    # this final quantity should be roughly peaking around 1-3 W/m2/sr/um

end


function create_SIF_scene(
    instrument::String,
    mode::String,
    loctime::GeolocationTime,
    vza::Real,
    uncertainty_function::Function,
    vnp_sd::VNPData
    )

    # solar azimuth (saa) at this point unused!
    this_sza, this_saa = calculate_solar_angles(loctime)

    # Obtain irradiance information (downwelling radiance at surface), already SZA-corrected
    # (unlike in L2 algorithms, where irradiance needs to be multiplied by
    #  mu0 to account for normal component)
    this_irradiance = calculate_BOA_irradiance(this_sza, 757.0)

    # Calculate PPFD for this scene
    this_PPFD = calculate_PPFD(this_sza)

    # black sky albedos from BRDFs?
    this_nir = calculate_reflectance(loctime, this_sza, "M7", vnp_sd)
    this_vis = calculate_reflectance(loctime, this_sza, "M5", vnp_sd)

    # These wavelengths are for VIIRS M5 and M7,
    # factor of 1000 takes us from W/m2/nm/sr to W/m2/um/sr
    this_nir_radiance = this_nir * calculate_BOA_irradiance(this_sza, 865.0) * 1000 / pi
    this_vis_radiance = this_vis * calculate_BOA_irradiance(this_sza, 672.0) * 1000 / pi

    this_ndvi = (this_nir - this_vis) / (this_nir + this_vis)

    this_nirv = this_ndvi * this_nir
    this_nirv_radiance = this_nirv * calculate_BOA_irradiance(this_sza, 865.0) * 1000 / pi

    # Reflectance at ~757 nm is roughly between VIIRS bands M5 and M7
    refl_M5 = calculate_reflectance(loctime, this_sza, "M5", vnp_sd)
    refl_M7 = calculate_reflectance(loctime, this_sza, "M7", vnp_sd)

    this_reflectance = 0.5 * (refl_M5 + refl_M7)

    # irradiance is in /nm, but we want /um
    this_TOA_radiance = this_irradiance * this_reflectance * 1000 / pi

    # This is in W/m2/sr/um
    this_sif = model_fluorescence(
        this_PPFD,
        25.0,
        200.0,
        209.0,
        757.0, # wavelength in nm
        this_ndvi
    )

    # The function needed to calculate the uncertainty must be supplied
    this_sif_ucert = uncertainty_function(this_TOA_radiance)

    # Scale SIF and uncertainty values up for 740nm
    # Note that the native calculation is still done for ~755nm,
    # however the resultws are scaled up to be comparable to e.g. TROPOMI

    this_sif *= 1.56
    this_sif_ucert *= 1.56


    # Noisified SIF
    if this_sif_ucert > 0.
        # Apply global seed variable to sampling
        # Random.seed!(jLRS_seed)
        this_measured_sif = rand(Normal(this_sif, this_sif_ucert))
    else
        this_measured_sif = this_sif
    end

    # At the moment no cloud or aerosol data,
    # could add ISCCP sampler in here
    this_od = 0.0

    this_scene = Scene(
        instrument,
        mode, # sampling mode comes from the instrument
        loctime, # location time comes from the instrument
        this_measured_sif,
        this_sif,
        this_sif_ucert,
        this_sza, # SZA comes from calculcations (via loctime)
        vza, # viewing zenith comes from the instrument,
        this_nir, # NIR reflectance
        this_vis, # VIS reflectance
        this_ndvi, # NDVI from VIIRS
        this_nirv, # NIRv calculated from NDVI * NIR
        this_nirv_radiance, # NIRv * L0(868nm)
        this_irradiance, # Irradiance at the surface and some ref. wl
        this_PPFD, # Integrated irradiance at PAR wavelengths 400nm to 700nm
        this_reflectance, # Reflectance comes from BRDF sampling and SZA
        this_od # optical depth may come from ISCCP one day
    )

    return this_scene

end
