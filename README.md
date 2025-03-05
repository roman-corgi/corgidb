# corgidb
Roman Coronagraph Instrument Target Database

Database Input Template
==============================
A database input is composed of two components: 	

1. The data (stored in a spreadsheet or any equivalent that can be read into pandas)	
2. A map between the data headers and db entries

To generate the map, use the Database Input Template.csv file. This file contains 6 columns, one header row, and 8 rows of sample data. You should delete the sample data but **leave the header row unchanged**. For each column in your data set, enter the following in one row of the template:
1. MY_COLNAME: The exact name of your column (case sensitive)
2. DB_COLNAME: The exact name of the destination column in the database (case sensitive)
3. UNITS: If the column includes numerical data with associated units, specify the unit (otherwise leave blank).  The unit must be a string that can be parsed by the astropy units module.  To test your unit string, you can run something like the following (in Python):
	```
	import astropy.units as u
	mystring = 'm s^-2'
	print(u.Unit(mystring))
	```
	If this does not produce an error, then your string is good. Note that the unit refers to the units of your data. If the column is a new one, this will also be the unit used in the new column in the database.  If you are mapping to an existing column, the data will be transformed into the database's native units, as required. 
4. NEW_KEY: If this is an existing column in the database, leave this blank (or set to FALSE) and skip the remaining columns.  If this is intended to be a new column, set to TRUE and fill in the next two columns:
	* DESCRIPTION: Write a brief but complete description of the column, including the units (where applicable)
   	* TABLE: The name of the intended table in the database (case sensitive)

The spreadsheet should not include anything other than 6 columns - do not include any notes or any other extraneous information. 

The current database schema may be accessed here: https://plandb.sioslab.com/docs/plandbschema/index.html 	

To begin your database input, save a copy of the Database Input Tempalte.csv file and fill it out.  Then, go to: https://github.com/roman-corgi/corgidb/issues, click on 'New Issue', select 'Database input', and then fill out the issue template and submit. 
