import pandas as pd
from astroquery.simbad import Simbad

def alias_check(name: str):
    """Checks the provided name agains the SIMBAD reference DB and returns all aliases 
    
    Args:
        name (str): input name to check for aliases
    
    Returns:
        aliases (pandas.Dataframe): 2 x n dataframe of all aliases for the given name and the main identifier
    """
    result = Simbad.query_object(name)
    if result is not None:
        main_id = result['MAIN_ID'][0]
        alias_res = Simbad.query_objectids(main_id)
        if alias_res is not None:
            alias_list = alias_res['ID'] .tolist()
            main_list = [main_id] * len(alias_list)
            data = {'Alias': alias_list, 'Main ID': main_list}
            dataframe = pd.DataFrame(data)
        else:
            print("ID query failed, check entry")
            dataframe = None
    else:
        print("Object not found in Simbad")
        dataframe = None
    return dataframe