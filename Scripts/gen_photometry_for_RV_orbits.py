from corgidb.photometry import (
    packagePhotometryData,
    loadPhotometryData,
    load_point_cloud,
)
import os
import glob
import re
import pandas as pd
import numpy as np


# load in model grid
# assumes that allphotdata.npz is in your working directory
# otherwise, update infile to specify full path to allphotdata.npz
infile = "allphotdata.npz"
if not os.path.exists(infile):
    # specify full path to db file
    dbfile = "AlbedoModels.db"  # from https://zenodo.org/records/2003949
    packagePhotometryData(dbfile=dbfile, outname=infile)

photdict = loadPhotometryData(infile=infile)


# load star catalog
stdatafile = "stdata_2025-02-25.p"
stars = pd.read_pickle(stdatafile)


# identify orbit fits available
# assumes that they're all in a subfolder of your working directory
# otherwise, updated to full path
orbfit_dir = "PointClouds_RVOnly"
orbfits = glob.glob(os.path.join(orbfit_dir, "*.pkl"))

# identify planet hosts
p = re.compile(
    r"^"
    r"(?P<prefix>[A-Za-z0-9]+(?:_[A-Za-z0-9]+)*)"  # tokens before the first date
    r"_"  # underscore before the first date
    r"(?P<start>\d{4}-\d{2}-\d{2})"  # start date: YYYY-MM-DD
    r"_to_"
    r"(?P<end>\d{4}-\d{2}-\d{2})"  # end date: YYYY-MM-DD
    r"_(?P<suffix>[A-Za-z0-9]+(?:_[A-Za-z0-9]+)*)"  # trailing tokens
    r"\.pkl"
    r"$"
)

targnames = np.array(
    [p.match(os.path.split(o)[-1]).group(1).replace("_", " ") for o in orbfits]
)
# one manual edit to our preferred naming
targnames[targnames == "pi Men"] = "HD 39091"

# verify that we know all of these stars
missing = list(set(targnames) - set(stars["st_name"]))
assert len(missing) == 0, f"Missing data for: {missing}"


for orbfit, targname in zip(orbfits, targnames):
    pointcloud = load_point_cloud(orbfit)


