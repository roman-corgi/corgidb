import pandas as pd
import time
from corgidb import alias_check

testList = ['Chi Uma', 'Omega Uma', '56 Uma', '58 Uma']
start_time = time.perf_counter()
all_data = alias_check.alias_check(testList)
stop_time = time.perf_counter()
print(all_data)
elapsed_time = stop_time - start_time
print(elapsed_time)