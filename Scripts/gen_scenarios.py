import corgietc  # noqa
import os
import json
import EXOSIMS.Prototypes.TargetList
import EXOSIMS.Prototypes.TimeKeeping
import copy
import astropy.units as u
import pandas


def gen_scenarios():
    """Create a dataframe containing the data for the Scenarios table

    Returns:
        pandas.DataFrame:
            The final table

    """

    # set up objects
    scriptfile = os.path.join(
        os.environ["CORGIETC_DATA_DIR"], "scripts", "CGI_Noise.json"
    )
    with open(scriptfile, "r") as f:
        specs = json.loads(f.read())

    TL = EXOSIMS.Prototypes.TargetList.TargetList(**copy.deepcopy(specs))
    OS = TL.OpticalSystem

    scenario_name = []
    minangsep_lamD = []
    maxangsep_lamD = []
    minangsep_as = []
    maxangsep_as = []
    lam = []
    bandpass = []
    delta_lam = []
    lamD_as = []

    for jj, mode in enumerate(OS.observingModes):
        scenario_name.append(mode["Scenario"])
        minangsep_as.append(mode["IWA"].to_value(u.arcsec))
        maxangsep_as.append(mode["OWA"].to_value(u.arcsec))
        lamD_as.append(mode["syst"]["input_angle_unit_value"].to_value(u.arcsec))
        minangsep_lamD.append(
            (mode["IWA"] / mode["syst"]["input_angle_unit_value"]).value
        )
        maxangsep_lamD.append(
            (mode["OWA"] / mode["syst"]["input_angle_unit_value"]).value
        )
        lam.append(mode["lam"].to_value(u.nm))
        bandpass.append(mode["BW"])
        delta_lam.append(mode["deltaLam"].to_value(u.nm))

    out = pandas.DataFrame(
        {
            "scenario_name": scenario_name,
            "minangsep_lamD": minangsep_lamD,
            "maxangsep_lamD": maxangsep_lamD,
            "minangsep_as": minangsep_as,
            "maxangsep_as": maxangsep_as,
            "lam": lam,
            "bandpass": bandpass,
            "delta_lam": delta_lam,
            "lamD_as": lamD_as,
        }
    )

    return out
