import pandas as pd
from astroquery.simbad import Simbad

def alias_check(name_list: list):
    """Checks the provided name agains the SIMBAD reference DB and returns all aliases 
    
    Args:
        name_list (list): list of input names to check for aliases
    
    Returns:
        aliasData (pandas.Dataframe): 2 x n dataframe of all aliases for the given name and the main identifier
    """
    alias_data = pd.DataFrame()
    missing_data = []
    for name in name_list:
        result = Simbad.query_object(name)
        if result is not None and len(result) > 0 :
            main_id = result['main_id'][0]
            alias_res = Simbad.query_objectids(main_id)
            if alias_res is not None:
                alias_list = alias_res['id'].tolist()
                main_list = [main_id] * len(alias_list)
                data = {'Alias': alias_list, 'Main ID': main_list}
                dataframe = pd.DataFrame(data)
                alias_data = pd.concat([alias_data, dataframe], ignore_index=True)
            else:
                print(f"Object {name} not found in Simbad MainIDs")
        else:
            print(f"Object {name} not found in Simbad")
            missing_data.append(name)
    missing_data = pd.DataFrame(missing_data)
    return alias_data, missing_data