import pandas as pd
from corgidb import alias_check

cars_df = pd.read_csv('car_targets.csv')
cars_df['Target Name'] = cars_df['Target Name'].str.split(',')
cars_df = cars_df.explode('Target Name').reset_index(drop=True)
cars_names_list = cars_df.values.tolist()
cars_names_list = [item for sublist in cars_names_list for item in sublist]
print(cars_names_list)
cars_alias_df, cars_missing_df = alias_check.alias_check(cars_names_list)
cars_alias_df.to_csv('cars_alias.csv')
cars_missing_df.to_csv('cars_missing.csv')