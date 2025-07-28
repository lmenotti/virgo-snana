# Startup instructions: 
To run, clone the repo, cd into the folder, and run process_data.py in your terminal. 
It will takes a minute, but should result in all snana files being generated in a new folder (snana_virgo_data).
Failure messages for individual parsers are ok, but keep an eye out for "WARNING: Could not parse file {file_name} with any available parser.", UNKNOWN bands, or unmapped bands (besides limits).


# File info: 
- Raw_virgo_data is where all the photometry is stored, per supernova. They are stored in .fit files, .csv files, and .txt files, as unedited as possible from the original source. Sources and information about each individual file can be found on the info spreadsheet.

- config.py: tells process_data.py which files to look at and what magnitude system to use for each supernova.

- parsers.py: collection of parsers that I use to get the data from each file.

- process_data.py: a script that parses through every supernova folder for data, parses the data, collects it, and generates a snana file for each supernova. a new folder will be created called snana_virgo_data where the snana files will be stored.

- supernova_info.csv: not uploaded yet, but will contain source info about every file stored in the repo, as well as info on degeneracies and edits.
