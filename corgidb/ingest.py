import pandas
import keyring
import getpass
from sqlalchemy import create_engine, text, types
import os
import numpy as np
import re
import warnings


def gen_engine(username, db="plandb", server="127.0.0.1"):
    """Creates a SQLAlchemy engine object.

    Saves password in local keyring and retrieves automatically for future
    connections.

    Args:
        username (str): The username for the database connection.
        db (str): The database name. Defaults to "plandb".
        server (str): The server address. Defaults to "127.0.0.1".

    Returns:
        sqlalchemy.engine.base.Engine:
            The configured SQLAlchemy engine object.
    """
    passwd = keyring.get_password(f"plandb_{server}_login", username)
    if passwd is None:
        newpass = True
        passwd = getpass.getpass(f"Password for user {username} on server {server}:\n")
    else:
        newpass = False

    engine = create_engine(
        f"mysql+pymysql://{username}:{passwd}@{server}/{db}",
        echo=False,
        pool_pre_ping=True,
    )

    connection = engine.connect()
    _ = connection.execute(text("SHOW TABLES"))
    connection.close()

    if newpass:
        keyring.set_password(f"plandb_{server}_login", username, passwd)

    return engine


def proc_col_req(fname, engine, comment="#"):
    """Processes a new table or key request from a spreadsheet.

    Args:
        fname (str): Full path to the input file (.csv, .xls, .xlsx).
        engine (sqlalchemy.engine.base.Engine): Engine object from `gen_engine`.
        comment (str): The comment character in the file. Defaults to "#".

    Raises:
        NotImplementedError: If the input file is not a supported type.
    """
    ext = os.path.splitext(fname)[-1].split(os.extsep)[-1].lower()
    if ext == "csv":
        data = pandas.read_csv(fname, comment=comment)
    elif ext in ["xlsx", "xls"]:
        data = pandas.read_excel(fname, comment=comment)
    else:
        raise NotImplementedError("Filetype must be csv, xls or xlsx.")

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

    for boolkey in ["NEW_KEY", "INDEX"]:
        data.loc[data[boolkey].isna(), boolkey] = False

    inds = (data["NEW_KEY"]) & (data["DB_COLNAME"].isna())
    data.loc[inds, "DB_COLNAME"] = data.loc[inds, "MY_COLNAME"]

    data.loc[data["SQL_DATATYPE"] == "STRING", "SQL_DATATYPE"] = "TEXT"

    assert not (data[reqcols].isna().values.any()), "Input data is missing entries."

    req_tables = data["TABLE"].unique()

    connection = engine.connect()
    tables = connection.execute(text("SHOW TABLES")).all()
    tables = [t[0] for t in tables]

    existing_tables = list(set(req_tables).intersection(tables))

    existing_keys = {}
    for t in existing_tables:
        tmp = connection.execute(text(f"SHOW COLUMNS IN {t}")).all()
        existing_keys[t] = np.array([k[0] for k in tmp])

    for t in existing_tables:
        tmp = data.loc[data["TABLE"] == t]
        newkeys = tmp.loc[~tmp["DB_COLNAME"].isin(existing_keys[t])]
        assert newkeys["NEW_KEY"].all(), (
            f"Some keys requested for table {t} do not exist in DB, "
            "but are not marked as new."
        )

        for jj, row in newkeys.iterrows():
            comm = (
                f"""ALTER TABLE {t} ADD COLUMN {row.DB_COLNAME} {row.SQL_DATATYPE} """
                f"""COMMENT "{row.DESCRIPTION}";"""
            )

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

    new_tables = list(set(req_tables) - set(existing_tables))

    for t in new_tables:
        tmp = data.loc[data["TABLE"] == t]

        txt = []
        for jj, row in tmp.iterrows():
            txt.append(
                f'''{row.DB_COLNAME} {row.SQL_DATATYPE} COMMENT "{row.DESCRIPTION}"'''
            )

        comm = f"CREATE TABLE {t} ({', '.join(txt)});"

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
    """Populates the Scenarios table.

    Args:
        data (pandas.DataFrame): DataFrame containing the table data.
        schema (pandas.DataFrame): DataFrame with column names ('Column') and
            comments ('Comments').
        engine (sqlalchemy.engine.base.Engine): The database engine.
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
    """Populates the SaturationCurves table.

    Args:
        data (pandas.DataFrame): DataFrame containing the table data.
        schema (pandas.DataFrame): DataFrame with column names ('Column') and
            comments ('Comments').
        engine (sqlalchemy.engine.base.Engine): The database engine.
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
    """Adds indexes to an existing table.

    Args:
        connection (sqlalchemy.engine.base.Connection): A database connection object.
        tablename (str): The name of the table to modify.
        indexes (list): A list of column names to be indexed.
    """
    _ = connection.execute(
        text(f"ALTER TABLE {tablename} ADD INDEX ({', '.join(indexes)})")
    )


def add_foreignkeys(connection, tablename, cols, foreignkeys):
    """Adds foreign key constraints to an existing table.

    Args:
        connection (sqlalchemy.engine.base.Connection): A database connection object.
        tablename (str): The name of the table to modify.
        cols (list): A list of column names in the current table to apply
            foreign keys to.
        foreignkeys (list): A list of foreign key specs, each in the
            format "ReferencedTable(ReferencedColumn)".
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
    """Adds comments to an existing table and sets indexes and foreign keys.

    Args:
        connection (sqlalchemy.engine.base.Connection): A database connection object.
        tablename (str): The name of the table to update.
        schema (pandas.DataFrame): DataFrame with column names ('Column') and
            the desired comments ('Comments').
    """
    indexes = schema.loc[schema["Index"] == 1, "Column"].values
    if len(indexes) > 0:
        add_indexes(connection, tablename, indexes)

    foreignkeys = schema.loc[~schema["ForeignKey"].isna()]
    if len(foreignkeys) > 0:
        add_foreignkeys(
            connection,
            tablename,
            foreignkeys["Column"].values,
            foreignkeys["ForeignKey"].values,
        )

    result = connection.execute(text(f"show create table {tablename}"))
    res = result.fetchall()
    res = res[0][1]
    res = res.split("\n")

    p = re.compile(r"`(\S+)`[\s\S]+")

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
            f'"{schema.loc[schema["Column"] == key, "Comments"].values[0]}";'
        )
        _ = connection.execute(text(comm))


def get_optimal_sql_datatypes(df: pandas.DataFrame) -> dict:
    """Generates optimal SQL datatypes for columns in a pandas DataFrame.

    Args:
        df (pandas.DataFrame): The input DataFrame.

    Returns:
        dict:
            A dictionary mapping column names to optimal SQL datatypes.
    """
    col_types = {}
    for col in df.columns:
        dtype = df[col].dtype
        if np.issubdtype(dtype, np.integer):
            min_val, max_val = df[col].min(), df[col].max()
            if -128 <= min_val and max_val <= 127:
                col_types[col] = "TINYINT"
            elif -32768 <= min_val and max_val <= 32767:
                col_types[col] = "SMALLINT"
            elif -2147483648 <= min_val and max_val <= 2147483647:
                col_types[col] = "INT"
            else:
                col_types[col] = "BIGINT"
        elif np.issubdtype(dtype, np.floating):
            col_types[col] = "DOUBLE"
        elif dtype == np.bool_:
            col_types[col] = "BOOLEAN"
        else:
            max_len = df[col].str.len().max()
            if pandas.isna(max_len):
                col_types[col] = "VARCHAR(255)"
            else:
                col_types[col] = f"VARCHAR({int(max_len)})"
    return col_types
