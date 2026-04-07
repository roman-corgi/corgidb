import pandas as pd
from corgidb import alias_check
from corgidb import ingest
from astroquery.simbad import Simbad
# script to generate the dataframe for all alias names from the stars table.

conn = ingest.gen_engine('plandb_user',db='plandb_scratch')
column = 'st_name'
table = 'Stars'
query = f"SELECT {column} FROM {table}"
names_df = pd.read_sql_query(query, conn)
names_list = names_df['st_name'].to_list()
alias_df, missing_df = alias_check.alias_check(names_list)
alias_df.to_csv('star_alias.csv')
missing_df.to_csv('missing_simbad_data.csv')
#TODO: rerun on plandb scratch and append code to include a list that cannot be resolved via simbad