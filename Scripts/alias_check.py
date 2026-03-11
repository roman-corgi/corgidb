import pandas as pd
import time
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
        main_id = result['main_id'][0]
        alias_res = Simbad.query_objectids(main_id)
        if alias_res is not None:
            alias_list = alias_res['id'] .tolist()
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

def int_alias_check(nameList: list):
    """Checks the provided name agains the SIMBAD reference DB and returns all aliases 
    
    Args:
        nameList (list): list of input names to check for aliases
    
    Returns:
        aliases (pandas.Dataframe): 2 x n dataframe of all aliases for the given name and the main identifier
    """
    aliasData = pd.DataFrame()
    for name in nameList:
        result = Simbad.query_object(name)
        if result is not None:
            main_id = result['main_id'][0]
            alias_res = Simbad.query_objectids(main_id)
            if alias_res is not None:
                alias_list = alias_res['id'] .tolist()
                main_list = [main_id] * len(alias_list)
                data = {'Alias': alias_list, 'Main ID': main_list}
                dataframe = pd.DataFrame(data)
                aliasData = pd.concat([aliasData, dataframe], )
            else:
                print(f"ID query for {name} failed, check entry")
                dataframe = None
                break
        else:
            print(f"Object {name} not found in Simbad")
            dataframe = None
            break
    return aliasData

# test if it is more efficient to loop over single names or to query a list and then loop through aliases 

testList = ['Chi Uma', 'Omega Uma', '56 Uma', '58 Uma']
start_time1 = time.perf_counter()
all_data = pd.DataFrame()
for star in testList:
    data = alias_check(star)
    all_data = pd.concat([all_data, data], ignore_index=True)
stop_time1 = time.perf_counter()
print(all_data)
elapsed_time1 = stop_time1 - start_time1
print(elapsed_time1)
start_time2 = time.perf_counter()
all_data2 = int_alias_check(testList)
stop_time2 = time.perf_counter()
print(all_data2)
elapsed_time2 = stop_time2 - start_time2
print(elapsed_time2)