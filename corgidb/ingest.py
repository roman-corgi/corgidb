import pandas
import keyring
import getpass
from sqlalchemy import create_engine, text, types
import os
import numpy as np
import re
import warnings


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

    return engine


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
        "TABLE",
        "UNITS",
        "NEW_KEY",
        "DESCRIPTION",
        "SQL_DATATYPE",
        "INDEX",
        "FOREIGNKEY",
    ]
    reqcols = [
        "MY_COLNAME",
        "DB_COLNAME",
        "TABLE",
        "NEW_KEY",
        "SQL_DATATYPE",
        "INDEX",
    ]

    assert (
        len(set(data.columns).symmetric_difference(expected_cols)) == 0
    ), f'Input data must contain the columns: {", ".join(expected_cols)} ONLY'

    # fill in missing DB cols for boolean keys
    for boolkey in ["NEW_KEY", "INDEX"]:
        data.loc[data[boolkey].isna(), boolkey] = False

    # new keys without db_colname entries should use my_colname
    inds = (data["NEW_KEY"]) & (data["DB_COLNAME"].isna())
    data.loc[inds, "DB_COLNAME"] = data.loc[inds, "MY_COLNAME"]

    # change any STRING datatypes to VARCHAR
    data.loc[data["SQL_DATATYPE"] == "STRING", "SQL_DATATYPE"] = "TEXT"

    # verify that all rows are complete
    assert not (data[reqcols].isna().values.any()), "Input data is missing entries."

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
            f"Some keys requested for table {t} do not exist in DB, "  # noqa
            "but are not marked as new."
        )

        for jj, row in newkeys.iterrows():
            comm = (
                f"""ALTER TABLE {t} ADD COLUMN {row.DB_COLNAME} {row.SQL_DATATYPE} """
                f"""COMMENT "{row.DESCRIPTION}";"""  # noqa
            )

            _ = connection.execute(text(comm))

        # add any requested indexes and foreignkeys
        indexes = tmp.loc[tmp["INDEX"], "DB_COLNAME"].values
        if len(indexes) > 0:
            add_indexes(connection, t, indexes)

        foreignkeys = tmp.loc[~tmp["FOREIGNKEY"].isna()]
        if len(foreignkeys) > 0:
            add_foreignkeys(
                connection,
                t,
                foreignkeys["DB_COLNAME"].values,
                foreignkeys["FOREIGNKEY"].values,
            )

    # these are the new tables
    new_tables = list(set(req_tables) - set(existing_tables))

    for t in new_tables:
        tmp = data.loc[data["TABLE"] == t]

        # generate create table text
        txt = []
        for jj, row in tmp.iterrows():
            txt.append(
                f'''{row.DB_COLNAME} {row.SQL_DATATYPE} COMMENT "{row.DESCRIPTION}"'''
            )

        comm = f"CREATE TABLE {t} ({', '.join(txt)});"  # noqa

        _ = connection.execute(text(comm))

        indexes = tmp.loc[tmp["INDEX"], "DB_COLNAME"].values
        if len(indexes) > 0:
            add_indexes(connection, t, indexes)

        foreignkeys = tmp.loc[~tmp["FOREIGNKEY"].isna()]
        if len(foreignkeys) > 0:
            add_foreignkeys(
                connection,
                t,
                foreignkeys["DB_COLNAME"].values,
                foreignkeys["FOREIGNKEY"].values,
            )


def gen_Scenarios_table(data, schema, engine):
    """Populate SaturationCurves table

    Args:
        data (pandas.DataFrame):
            Table data
        schema (pandas.DataFrame):
            Table of column names ('Column') and comments ('Comments')
        engine (sqlalchemy.engine.base.Engine):
            Engine

    Returns:
        None

    """
    connection = engine.connect()
    namemxchar = np.array([len(n) for n in data["scenario_name"].values]).max()

    _ = data.to_sql(
        "Scenarios",
        connection,
        chunksize=100,
        if_exists="replace",
        dtype={
            "scenario_name": types.String(namemxchar),
        },
        index=False,
    )

    updateSQLschema(connection, "Scenarios", schema)


def gen_SaturationCurves_table(data, schema, engine):
    """Populate SaturationCurves table

    Args:
        data (pandas.DataFrame):
            Table data
        schema (pandas.DataFrame):
            Table of column names ('Column') and comments ('Comments')
        engine (sqlalchemy.engine.base.Engine):
            Engine

    Returns:
        None

    """
    connection = engine.connect()
    namemxchar = np.array([len(n) for n in data["scenario_name"].values]).max()

    _ = data.to_sql(
        "SaturationCurves",
        connection,
        chunksize=100,
        if_exists="replace",
        dtype={
            "scenario_name": types.String(namemxchar),
        },
        index=False,
    )

    updateSQLschema(connection, "SaturationCurves", schema)


def add_indexes(connection, tablename, indexes):
    """Add indexes to an existing table

    Args:
        connection (sqlalchemy.engine.base.Connection):
            connection object
        tablename (str):
            Name of table
        indexes (list):
            List of columnnames to make into indexes

    Returns:
        None

    """

    _ = connection.execute(
        text(f"ALTER TABLE {tablename} ADD INDEX ({', '.join(indexes)})")
    )


def add_foreignkeys(connection, tablename, cols, foreignkeys):
    """Add foreign keys to table

    Args:
        connection (sqlalchemy.engine.base.Connection):
            connection object
        tablename (str):
            Name of table
        cols (list):
            List of columnnames to add as foreignkeys
        foreignkeys (list):
            List of foreign key specifications of the form "Table(column)"

    Returns:
        None


    """

    for col, fkey in zip(cols, foreignkeys):
        _ = connection.execute(
            text(
                f"ALTER TABLE {tablename} ADD FOREIGN KEY ({col}) "
                f"REFERENCES {fkey} ON DELETE NO ACTION "
                "ON UPDATE NO ACTION"
            )
        )


def updateSQLschema(connection, tablename, schema):
    """Add comments to an existing table and set indexes and foreign keys

    Args:
        connection (sqlalchemy.engine.base.Connection):
            connection object
        tablename (str):
            Name of table
        schema (pandas.DataFrame):
            Table of column names ('Column') and comments ('Comments')

    Returns:
        None

    """

    # handle indexes
    indexes = schema.loc[schema["Index"] == 1, "Column"].values
    if len(indexes) > 0:
        add_indexes(connection, tablename, indexes)

    # handle foregin keys
    foreignkeys = schema.loc[~schema["ForeignKey"].isna()]
    if len(foreignkeys) > 0:
        add_foreignkeys(
            connection,
            tablename,
            foreignkeys["Column"].values,
            foreignkeys["ForeignKey"].values,
        )

    # grab the original table definition
    result = connection.execute(text(f"show create table {tablename}"))
    res = result.fetchall()
    res = res[0][1]
    res = res.split("\n")

    # define regex for parsing col defs
    p = re.compile(r"`(\S+)`[\s\S]+")

    # loop through and find all col definitions without comments
    keys = []
    defs = []
    for r in res:
        r = r.strip().strip(",")
        if "COMMENT" in r:
            continue
        m = p.match(r)
        if m:
            keys.append(m.groups()[0])
            defs.append(r)

    # compare schema and table
    missing_from_schema = list(set(keys) - set(schema["Column"].values))
    if len(missing_from_schema) > 0:
        warnings.warn(
            (
                "Columns present in table but missing from schema: "
                f"{','.join(missing_from_schema)}"
            )
        )

    missing_from_table = list(set(schema["Column"].values) - set(keys))
    if len(missing_from_table) > 0:
        schema = schema.loc[~schema["Column"].isin(missing_from_table)].reset_index(
            drop=True
        )
        warnings.warn(
            (
                "Columns present in schema but missing from table (or already have "
                f"comments): {','.join(missing_from_table)}"
            )
        )

    for key, d in zip(keys, defs):
        comm = (
            f"ALTER TABLE `{tablename}` CHANGE `{key}` {d} COMMENT "
            f'"{schema.loc[schema["Column"] == key, "Comments"].values[0]}";'  # noqa
        )
        _ = connection.execute(text(comm))
