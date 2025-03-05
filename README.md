# corgidb
Roman Coronagraph Instrument Target Database

Database Input Template
==============================
A database input is composed of two components: 	

1. The data (stored in a spreadsheet or any equivalent that can be read into pandas)	
2. A map between the data headers and db entries

To generate the map, use the Database Input Template.csv file. This file contains 6 columns, one header row, and 8 rows of sample data. 
	
For existing keys, the columns NEW_KEY, DESCRIPTION, TABLE should be blank	
For new keys: 	
	set NEW_KEY to TRUE
	set DESCRIPTION to a detailed description of the key, including units (as appropriate)
	set TABLE to the intended TABLE for the key
	
Current schema may be accessed here: https://plandb.sioslab.com/docs/plandbschema/index.html 	
	
This spreadsheet should not include anything other than 5 columns (i.e., these notes should be deleted)	
	
	
Set UNITS to string that is parseable by astropy.units	
If not units leave it blank	
