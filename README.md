# corgidb
Roman Coronagraph Instrument Target Database

Database Input Template
==============================
A database input is composed of two components: 	

1. The data (stored in a spreadsheet or any equivalent that can be read into pandas)	
2. A map between the data headers and db entries

To generate the map, use the Database Input Template.csv file. This file contains 9 columns, one header row, and 8 rows of sample data. You should delete the sample data but **leave the header row unchanged and do not put anything in any other columns**. 

For each column in your data set, enter the following in one row of the template:
1. MY_COLNAME: The exact name of your column (case sensitive)
2. DB_COLNAME: The exact name of the destination column in the database (case sensitive).  If this entry is left blank, the value in MY_COLNAME will be copied over.
3. TABLE: The name of the table (either existing or new)
4. UNITS: If the column includes numerical data with associated units, specify the unit (otherwise leave blank).  The unit must be a string that can be parsed by the astropy units module.  To test your unit string, you can run something like the following (in Python):
	```
	import astropy.units as u
	mystring = 'm s^-2'
	print(u.Unit(mystring))
	```
	If this does not produce an error, then your string is good. Note that the unit refers to the units of your data. If the column is a new one, this will also be the unit used in the new column in the database.  If you are mapping to an existing column, the data will be transformed into the database's native units, as required. 
5. NEW_KEY: If this is an existing column in the database, leave this blank (or set to FALSE) and skip the remaining columns.  If this is intended to be a new column, set to TRUE and fill in the remaining columns:
6. DESCRIPTION: Write a brief but complete description of the column, including the units (where applicable)
7. SQL_DATATYPE: The SQL datatype to use for a new column.  See here for all available types: https://www.w3schools.com/sql/sql_datatypes.asp. Select the smallest type possible (e.g., do not use DOUBLE if INT will work). 
8. INDEX: Set to TRUE if you wish for this column to be an index (e.g., if this is the main thing you'll be querying on in this table). Set to FALSE or leave blank otherwise. 
9. FOREIGNKEY: If this column is an index from a different table, specify here as: `Table(column)`.  For example, if you wish for a column in your new table to reference st_id in table Stars, then the value in FOREIGNKEY for that column should be `Stars(st_id)`.

The spreadsheet should not include anything other than 9 columns. If you wish to add comments, add them in their own separate rows (do not leave any whitespace rows between data and comments).  All comment rows must have # as the leading character.  Do **not** use # in any data fields. 

The current database schema may be accessed here: https://plandb.sioslab.com/docs/plandbschema/index.html 	

To begin your database input, save a copy of the Database Input Tempalte.csv file and fill it out.  Then, go to: https://github.com/roman-corgi/corgidb/issues, click on 'New Issue', select 'Database input', and then fill out the issue template and submit. 

Database Snapshots
==============================
To simplify use, we will host snapshots of the databasde in Python pickle files that can be downloaded and directly manipulated. One snapshot is currently on the order of 2 GB. The snapshots will be available here: https://drive.google.com/drive/u/1/folders/1ibN_jW0ToOteTbTqz0niV4Adh-aWj9rU
