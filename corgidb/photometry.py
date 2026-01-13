import os
import pickle
import numpy as np
from scipy.interpolate import RectBivariateSpline, interp1d
from sqlalchemy import create_engine
import pandas as pd


def load_point_cloud(
    fname,
    broadcast_arrays=True,
):

    # Load pickle
    with open(fname, "rb") as f:
        point_cloud = pickle.load(f)

    # recast data to consistent array size
    if broadcast_arrays:
        arr_shape = point_cloud["sep_mas"].shape
        for param, arr in point_cloud.items():
            if arr.shape == (arr_shape[0],):
                point_cloud[param] = np.full(arr_shape, arr[:, np.newaxis])
            elif arr.shape == (arr_shape[1],):
                point_cloud[param] = np.full(arr_shape, arr)
            elif arr.shape != arr_shape:
                raise ValueError(
                    f"param {param} has unexpected array shape {arr.shape}, "
                    "when sep_mas has shape {arr_shape}."
                )

    return point_cloud


def packagePhotometryData(dbfile="AlbedoModels.db", outname="allphotdata.npz"):
    """
    Read photometry data from database of grid models and repackage into single ndarray
    of data

    Download db from https://zenodo.org/records/2003949
    """

    # grab photometry data
    enginel = create_engine("sqlite:///" + dbfile)

    # getting values
    meta_alb = pd.read_sql_table("header", enginel)
    metallicities = meta_alb.metallicity.unique()
    metallicities.sort()
    betas = meta_alb.phase.unique()
    betas.sort()
    dists = meta_alb.distance.unique()
    dists.sort()
    clouds = meta_alb.cloud.unique()
    clouds.sort()
    cloudstr = clouds.astype(str)
    for j in range(len(cloudstr)):
        cloudstr[j] = "f" + cloudstr[j]
    cloudstr[cloudstr == "f0.0"] = "NC"
    cloudstr[cloudstr == "f1.0"] = "f1"
    cloudstr[cloudstr == "f3.0"] = "f3"
    cloudstr[cloudstr == "f6.0"] = "f6"

    tmp = pd.read_sql_table("g25_t150_m0.0_d0.5_NC_phang000", enginel)
    wavelns = tmp.WAVELN.values

    allphotdata = np.zeros(
        (metallicities.size, dists.size, clouds.size, betas.size, wavelns.size)
    )
    for i, fe in enumerate(metallicities):
        basename = "g25_t150_m" + str(fe) + "_d"
        for j, d in enumerate(dists):
            basename2 = basename + str(d) + "_"
            for k, cloud in enumerate(clouds):
                basename3 = basename2 + cloudstr[k] + "_phang"
                print(basename3)
                for l, beta in enumerate(betas):
                    name = basename3 + "%03d" % beta
                    try:
                        tmp = pd.read_sql_table(name, enginel)
                    except:  # noqa
                        print("Missing: %s" % name)
                        allphotdata[i, j, k, l, :] = np.nan
                        continue
                    pvals = tmp["GEOMALB"].values
                    if len(tmp) != len(wavelns):
                        missing = list(set(wavelns) - set(tmp.WAVELN.values))
                        inds = np.searchsorted(tmp["WAVELN"].values, missing)
                        pvals = np.insert(pvals, inds, np.nan)
                        assert np.isnan(pvals[wavelns == missing[0]])
                        print("Filled value: %s, %s" % (name, missing))
                    allphotdata[i, j, k, l, :] = pvals

    # patch individual nans
    for i, fe in enumerate(metallicities):
        for j, d in enumerate(dists):
            for k, cloud in enumerate(clouds):
                for l, beta in enumerate(betas):
                    nans = np.isnan(allphotdata[i, j, k, l, :])
                    if np.any(nans) & ~np.all(nans):
                        tmp = interp1d(
                            wavelns[~nans], allphotdata[i, j, k, l, ~nans], kind="cubic"
                        )
                        allphotdata[i, j, k, l, nans] = tmp(wavelns[nans])

    np.savez(
        outname,
        metallicities=metallicities,
        dists=dists,
        clouds=clouds,
        cloudstr=cloudstr,
        betas=betas,
        wavelns=wavelns,
        allphotdata=allphotdata,
    )


def loadPhotometryData(infile="allphotdata.npz"):
    """
    Read stored photometry data from disk and generate interpolants over data
    """

    tmp = np.load(infile)
    allphotdata = tmp["allphotdata"]
    clouds = tmp["clouds"]
    cloudstr = tmp["cloudstr"].astype(str)
    wavelns = tmp["wavelns"]
    betas = tmp["betas"]
    dists = tmp["dists"]
    metallicities = tmp["metallicities"]

    def makeninterp(vals):
        ii = interp1d(
            vals,
            vals,
            kind="nearest",
            bounds_error=False,
            fill_value=(vals.min(), vals.max()),
        )
        return ii

    distinterp = makeninterp(dists)
    betainterp = makeninterp(betas)
    feinterp = makeninterp(metallicities)
    cloudinterp = makeninterp(clouds)

    photinterps2 = {}
    quadinterps = {}
    for i, fe in enumerate(metallicities):
        photinterps2[fe] = {}
        quadinterps[fe] = {}
        for j, d in enumerate(dists):
            photinterps2[fe][d] = {}
            quadinterps[fe][d] = {}
            for k, cloud in enumerate(clouds):
                if np.any(np.isnan(allphotdata[i, j, k, :, :])):
                    # remove whole rows of betas
                    goodbetas = np.array(
                        list(
                            set(range(len(betas)))
                            - set(
                                np.unique(
                                    np.where(np.isnan(allphotdata[i, j, k, :, :]))[0]
                                )
                            )
                        )
                    )
                    photinterps2[fe][d][cloud] = RectBivariateSpline(
                        betas[goodbetas], wavelns, allphotdata[i, j, k, goodbetas, :]
                    )
                else:
                    photinterps2[fe][d][cloud] = RectBivariateSpline(
                        betas, wavelns, allphotdata[i, j, k, :, :]
                    )
                quadinterps[fe][d][cloud] = interp1d(
                    wavelns, allphotdata[i, j, k, 9, :].flatten()
                )

    return {
        "allphotdata": allphotdata,
        "clouds": clouds,
        "cloudstr": cloudstr,
        "wavelns": wavelns,
        "betas": betas,
        "dists": dists,
        "metallicities": metallicities,
        "distinterp": distinterp,
        "betainterp": betainterp,
        "feinterp": feinterp,
        "cloudinterp": cloudinterp,
        "photinterps": photinterps2,
        "quadinterps": quadinterps,
    }


# Generates fsed based on a random number: 0 <= num < 1
def get_fsed(num):
    """Generate random value of f_sed based on distribution provided by Mark Marley:
    f_sed             Frequency
    0.000000         0.099
    0.010000         0.001
    0.030000         0.005
    0.100000         0.010
    0.300000         0.025
    1.000000         0.280
    3.000000         0.300
    6.000000         0.280

    Input num is a uniform random value between 0 and 1
    """

    if num < 0.099:
        r = 0
    elif num < 0.1:
        r = 0.01
    elif num < 0.105:
        r = 0.03
    elif num < 0.115:
        r = 0.1
    elif num < 0.14:
        r = 0.3
    elif num < 0.42:
        r = 1
    elif num < 0.72:
        r = 3
    else:
        r = 6
    return float(r)
