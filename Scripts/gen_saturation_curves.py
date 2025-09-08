import corgietc  # noqa
import os
import json
import EXOSIMS.Prototypes.TargetList
import EXOSIMS.Prototypes.TimeKeeping
import copy
import astropy.units as u
import numpy as np
import pandas
from synphot import units


def gen_saturation_curves():
    """Create a dataframe containing the data for the SaturationCurves table

    Returns:
        pandas.DataFrame:
            The final table

    """

    # define units
    PHOTLAM_sr_decomposed_val = (1 * units.PHOTLAM / u.sr).to_value(
        u.ph / u.s / u.arcsec**2 / u.cm**2 / u.nm
    )
    F0_unit = u.ph / u.s / u.cm**2
    deltaLam_unit = u.nm
    inv_arcsec2 = 1 / u.arcsec**2

    # set up objects
    scriptfile = os.path.join(
        os.environ["CORGIETC_DATA_DIR"], "scripts", "CGI_Noise.json"
    )
    with open(scriptfile, "r") as f:
        specs = json.loads(f.read())

    TL = EXOSIMS.Prototypes.TargetList.TargetList(**copy.deepcopy(specs))
    OS = TL.OpticalSystem

    # set fixed inputs
    sInds = 0
    Izod = TL.ZodiacalLight.zodi_intensity_at_location(135 * u.deg, 30 * u.deg)
    npoints = 100

    # loop through all modes and populate arrays
    scenario_name = []
    r_lamD = []
    r_as = []
    contrast = []
    dMag = []
    t_int_hr_99percent_V5 = []

    for jj, mode in enumerate(OS.observingModes):
        # compute local zodi and exozodi
        Izod_color = Izod * TL.ZodiacalLight.zodi_color_correction_factor(mode["lam"])

        factor = (
            units.convert_flux(mode["lam"], Izod_color * u.sr, units.PHOTLAM).value
            / Izod_color.value
        )

        Izod_photons = Izod_color.value * factor * PHOTLAM_sr_decomposed_val

        fZ = (
            Izod_photons
            / (mode["F0"].to_value(F0_unit) / mode["deltaLam"].to_value(deltaLam_unit))
        ) * inv_arcsec2

        JEZ = TL.JEZ0[mode["hex"]][sInds]

        # populate scenario name
        scenario_name += [mode["Scenario"]] * npoints

        # populate working angles
        WAs = (
            np.linspace(mode["IWA"].value, mode["OWA"].value, npoints)
            * mode["IWA"].unit
        )

        r_as += WAs.to_value(u.arcsec).tolist()
        r_lamD += (WAs / mode["syst"]["input_angle_unit_value"]).value.tolist()

        # compute saturation dMag and contrast
        sat_dMags = OS.calc_saturation_dMag(
            TL,
            [sInds] * len(WAs),
            np.repeat(fZ, len(WAs)),
            np.repeat(JEZ, len(WAs)),
            WAs,
            mode,
        )
        dMag += sat_dMags.tolist()
        contrast += (10 ** (-0.4 * sat_dMags)).tolist()

        # compute integration time for 99% of sat dMag
        itimes = OS.calc_intTime(
            TL,
            [sInds] * len(WAs),
            np.repeat(fZ, len(WAs)),
            np.repeat(JEZ, len(WAs)),
            sat_dMags * 0.99,
            WAs,
            mode,
        )
        t_int_hr_99percent_V5 += itimes.to_value(u.hour).tolist()

        # end mode loop

    out = pandas.DataFrame(
        {
            "scenario_name": scenario_name,
            "r_lamD": r_lamD,
            "r_as": r_as,
            "contrast": contrast,
            "dMag": dMag,
            "t_int_hr_99percent_V5": t_int_hr_99percent_V5,
        }
    )

    return out
