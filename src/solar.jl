using Interpolations

# Loads the data table
solar_h5 = h5open(String(@__DIR__) * "/../data/solar_irradiance_table.h5", "r")

# This structure will hold the necessary data which other
# functions then operate on to return e.g. integrated PPDF
# or the irradiance at 757nm etc.

struct SolarData
    sza::Vector{Real}
    wavelength::Vector{Real}
    solar::Vector{Real}

    GHI::ScaledInterpolation
    tgas::ScaledInterpolation
    trayleigh::ScaledInterpolation
    PPFD::ScaledInterpolation
end

# Construct interpolation objects for the three
# main variables: "direct", "diffuse", "transmittance"

function create_solar_data(h5) :: SolarData

    # Interpolations uses ranges rather than vectors/arrays
    # to scale the grid axes. So here we check if the SZA and
    # wavelength variables inside the HDF5 file can be properly
    # represented as a range object.

    sza_vec = h5["sza"][:]
    wl_vec = h5["wavelength"][:]

    sza_inc = sza_vec[2] - sza_vec[1]
    sza_range = minimum(sza_vec):sza_inc:maximum(sza_vec)

    if !(sza_range == sza_vec)
        @error "Sorry - SZA vector cannot be represented as range."
    end

    wl_inc = wl_vec[2] - wl_vec[1]
    wl_range = minimum(wl_vec):wl_inc:maximum(wl_vec)

    if !(wl_range == wl_vec)
        @error "Sorry - wavelength vector cannot be represented as range."
    end


    sdict = Dict{String, ScaledInterpolation}()

    # Produce interpolation objects for the tabulated data
    for var in ["GHI", "tgas", "trayleigh"]

        tmp_itp = interpolate(solar_h5[var][:,:], BSpline(Linear()))
        tmp_sitp = Interpolations.scale(tmp_itp, wl_range, sza_range)

        sdict[var] = tmp_sitp
    end

    # Produce an interpolation object for the newly
    # pre-calculated PPFD quantity

    PPFD_table = zeros(Float64, length(sza_vec))
    hcNA_inv = 8.359347229111778 # 1 / (h * c * NA)

    # Use a wavelength iterator
    # (which can be different from the wavelengths supplied in the file..)
    # PPFD is integrated along the visible spectrum, so we say 400 to 700
    this_wl = 400:1:700
    for (i, this_sza) in enumerate(sza_vec)
        rad = sdict["GHI"].(this_wl, Ref(this_sza)) # This is in W/m2
        f = rad .* this_wl # why do we have to do this?

        PPFD = 0.0
        # Integrate spectrum from 400nm to 700nm
        for j in 2:length(this_wl)
            PPFD += 0.5 * (f[j-1] + f[j]) * (this_wl[j] - this_wl[j-1])
        end

        # Factor of 1000 needed to get to umol/s/m2
        PPFD *= hcNA_inv / 1000
        PPFD_table[i] = PPFD
    end

    PPFD_itp = interpolate(PPFD_table, BSpline(Linear()))
    sdict["PPFD"] = Interpolations.scale(PPFD_itp, sza_range)

    return SolarData(
        sza_vec,
        wl_vec,
        h5["solar"][:],
        sdict["GHI"],
        sdict["tgas"],
        sdict["trayleigh"],
        sdict["PPFD"]
    )
end

# This creates the solar data object
solardata = create_solar_data(solar_h5)

function calculate_BOA_irradiance(sza::Real, wavelength::Real; solardata=solardata)

    if isnan(sza)
        return NaN
    else
        return solardata.GHI(wavelength, sza)
    end
end

function calculate_PPFD(sza::Real; solardata=solardata)
    if isnan(sza)
        return NaN
    else
        return solardata.PPFD(sza)
    end
end
