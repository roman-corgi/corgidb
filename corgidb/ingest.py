import pandas
import keyring
import getpass
from sqlalchemy import create_engine, text
import os
import numpy as np


def gen_engine(username, db="plandb", server="127.0.0.1"):
    """Create an SQLalcehmy engine object. Saves password in local keyring and
    retrieves automatically for future connections.

    Args:
        username (str):
            username
        db (str):
            database name (defaults to plandb)
        server (str):
            Address of server to connec to. Defaults to 127.0.0.1

    Returns:
        sqlalchemy.engine.base.Engine:
            The engine object

    """

    # grab password from keyring (or save if not yet saved)
    passwd = keyring.get_password(f"plandb_{server}_login", username)
    if passwd is None:
        newpass = True
        passwd = getpass.getpass(f"Password for user {username} on server {server}:\n")
    else:
        newpass = False

    # create engine
    engine = create_engine(
        f"mysql+pymysql://{username}:{passwd}@{server}/{db}",
        echo=False,
        pool_pre_ping=True,
    )

    # open a connection to make sure everything works
    connection = engine.connect()
    _ = connection.execute(text("SHOW TABLES"))
    connection.close()

    if newpass:
        keyring.set_password(f"plandb_{server}_login", username, passwd)


def proc_col_req(fname, engine, comment="#"):
    """Process new table/key request spreadsheet

    Args:
        fname (str):
            Full path to input file on disk.
        engine (sqlalchemy.engine.base.Engine):
            Engine object returned by `gen_engine`
        comment (str):
            Comment symbol to assume.  Defaults to #

    Raises:
        NotImplementedError:
            For unsupported filetypes

    """

    # read in column description
    ext = os.path.splitext(fname)[-1].split(os.extsep)[-1].lower()
    if ext == "csv":
        data = pandas.read_csv(fname, comment=comment)
    elif ext in ["xlsx", "xls"]:
        data = pandas.read_excel(fname, comment=comment)
    else:
        raise NotImplementedError("Filetype must be csv, xls or xlsx.")

    # check for correct columns
    expected_cols = [
        "MY_COLNAME",
        "DB_COLNAME",
        "UNITS",
        "NEW_KEY",
        "DESCRIPTION",
        "TABLE",
    ]
    assert (
        len(set(data.columns).symmetric_difference(expected_cols)) == 0
    ), f'Input data must contain the columns: {", ".join(expected_cols)} ONLY'

    # verify that all rows are complete
    assert not (
        data[expected_cols[:3]].isna().values.any()
    ), "Input data is missing entries."

    # split UNITS col
    sqlu = []
    physu = []
    for val in data["UNITS"].values:
        if "," in val:
            tmp = val.split(",")
            sqlu.append(tmp[0].strip())
            physu.append(tmp[1].strip())
        else:
            sqlu.append(val)
            physu.append(None)
    data["UNITS"] = sqlu
    data["PHYSICALUNIT"] = physu

    # find all requested tables
    req_tables = data["TABLE"].unique()

    # grab all relevant tables and their current contents
    connection = engine.connect()
    tables = connection.execute(text("SHOW TABLES")).all()
    tables = [t[0] for t in tables]

    # these are the tables we are augmenting
    existing_tables = list(set(req_tables).intersection(tables))

    existing_keys = {}
    for t in existing_tables:
        tmp = connection.execute(text(f"SHOW COLUMNS IN {t}")).all()
        existing_keys[t] = np.array([k[0] for k in tmp])

    # identify all new requested keys and augemnt the tables
    for t in existing_tables:
        tmp = data.loc[data["TABLE"] == t]
        newkeys = tmp.loc[~tmp["DB_COLNAME"].isin(existing_keys[t])]
        assert newkeys["NEW_KEY"].all(), (
            f"Some keys requested for table {t} do not exist in DB, "
            "but are not marked as new."
        )

        for jj, row in newkeys.iterrows():
            if row.UNITS == "STRING":
                unit = "TEXT"
            else:
                unit = row.UNITS
            comm = (
                f"""ALTER TABLE {t} ADD COLUMN {row.DB_COLNAME} {unit} """
                f"""COMMENT "{row.DESCRIPTION}";"""
            )

            _ = connection.execute(text(comm))

    # these are the new tables
    new_tables = list(set(req_tables) - set(existing_tables))

    for t in new_tables:
        tmp = data.loc[data["TABLE"] == t]

        # generate create table text
        txt = []
        for jj, row in tmp.iterrows():
            if row.UNITS == "STRING":
                unit = "TEXT"
            else:
                unit = row.UNITS
            txt.append(f'''{row.DB_COLNAME} {unit} COMMENT "{row.DESCRIPTION}"''')

        comm = f"CREATE TABLE {t} ({', '.join(txt)});"

        _ = connection.execute(text(comm))
