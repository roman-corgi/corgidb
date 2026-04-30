from collections import Counter
import pandas as pd
from astroquery.simbad import Simbad


def alias_check(name_list: list, debug: bool=False):
    """Checks the provided name agains the SIMBAD reference DB and returns all aliases 
    
    Args:
        name_list (list): list of input names to check for aliases
        debug (bool): flag to print debug info about the data queried and returned
    
    Returns:
        alias_data (pandas.DataFrame): 2 x n dataframe of all aliases for the 
                                        given name and the main identifier
        missing_data (pandas.DataFrame): 2 x n dataframe of all st_name
                                         entries that are not in simbad
    """
    alias_data = pd.DataFrame()
    missing_data = []
    debug_data = {}
    checked_list = []
    for name in name_list:
        result = Simbad.query_object(name)
        if result is not None and len(result) > 0 :
            main_id = result['main_id'][0]
            if main_id not in checked_list:
                checked_list.append(main_id)
                alias_res = Simbad.query_objectids(main_id)
                if alias_res is not None:
                    alias_list = alias_res['id'].tolist()
                    main_list = [main_id] * len(alias_list)
                    if debug is True:
                        name_id = [name] * len(alias_list)
                        data = {'Alias': alias_list, 'Main ID': main_list, 'st_name': name_id}
                        dataframe = pd.DataFrame(data)
                        alias_data = pd.concat([alias_data, dataframe], ignore_index=True)
                    else:
                        data = {'Alias': alias_list, 'Main ID': main_list}
                        dataframe = pd.DataFrame(data)
                        alias_data = pd.concat([alias_data, dataframe], ignore_index=True)
                else:
                    print(f"Object {name} not found in Simbad ObjectIDs Query")
                    data = {'Alias': main_id, 'Main ID': main_id}
                    dataframe = pd.DataFrame(data)
                    alias_data = pd.concat([alias_data, dataframe], ignore_index=True)
            else:
                print(f"Skipping duplicate alias resolution for {main_id} with st_name {name}")
        else:
            print(f"Object {name} not found in Simbad")
            missing_data.append(name)
    missing_data = pd.DataFrame(missing_data)
    counts = Counter(checked_list)
    duplicates = [item for item, count in counts.items() if count > 1]
    st_missing = len(name_list) - len(checked_list)
    debug_data["duplicate_names"] = duplicates
    debug_data["st_name_missing"] = st_missing
    if debug is True:
        print(debug_data)
    return alias_data, missing_data
