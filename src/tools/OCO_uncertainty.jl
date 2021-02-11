using HDF5
using LsqFit
using Glob

function derive_footprint_fits(fname::String)

    nc = h5open(fname, "r")

    # Find out what footprints we have
    fp = nc["Metadata/FootprintId"][:]
    un_fp = unique(nc["Metadata/FootprintId"][:])
    sort!(un_fp)

    SIF_ucert = nc["SIF_Uncertainty_740nm"][:]
    cont = nc["Science/continuum_radiance_757nm"][:]
    qual = nc["Quality_Flag"][:]

    output = Dict()

    @. model(x, p) = p[1] + p[2] * x

    for this_fp in un_fp
        idx = findall((qual .== 0) .& (fp .== this_fp))
        fit = curve_fit(model, cont[idx], (SIF_ucert[idx]).^2,
                        [0.0, 1.0],
                        lower=[0.0, 0.0], upper=[Inf, Inf])

        output[this_fp] = coef(fit)
    end

    return output

end
