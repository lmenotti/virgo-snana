import io
import re
import numpy as np
import pandas as pd
from astropy.io import fits, ascii
from astropy.table import Table
from astropy.time import Time
from datetime import datetime


def _sanitize_and_standardize(tbl, column_map):
    """helper function to rename columns, ensure numeric types, and filter bad rows."""
    for generic_name, specific_name in column_map.items():
        if specific_name in tbl.colnames:
            tbl.rename_column(specific_name, generic_name)
            
    for col_name in ['time', 'mag', 'magerr']:
        if col_name in tbl.colnames:
            tbl[col_name] = pd.to_numeric(tbl[col_name], errors='coerce')
        elif col_name == 'magerr':
            tbl['magerr'] = np.nan

    if 'time' not in tbl.colnames or 'mag' not in tbl.colnames: return None
    tbl = tbl[~np.isnan(tbl['time']) & ~np.isnan(tbl['mag'])]
    
    if 'band' not in tbl.colnames: tbl['band'] = 'UNKNOWN'
    if 'reference' not in tbl.colnames: tbl['reference'] = 'N/A'
    
    tbl['band'] = tbl['band'].astype(str)
    
    return tbl['time', 'mag', 'magerr', 'band', 'reference'] if len(tbl) > 0 else None

def parse_txt(file_path):
    """
    Specifically for SN1960F/photometry_source2_data1.txt.
    Parses a tab-separated file with known headers.
    """
    try:
        with open(file_path, 'r') as f:
            header = f.readline()
        if not all(x in header for x in ['Julian Date', 'Gregorian Day', 'Magnitude', 'Indmag and Band']):
            return None

        # Use tab separator to match file format
        df = pd.read_csv(file_path, sep='\t')
        df = df.rename(columns={
            'Julian Date': 'time',
            'Magnitude': 'mag',
            'Indmag and Band': 'band'
        })

        tbl = Table.from_pandas(df)
        return _sanitize_and_standardize(tbl, {})
    except Exception as e:
        print(f"Failed to parse: {e}")
        return None

def parse_fit_vizier(file_path):
    """Parses VizieR FITS files, ignoring non-standard header issues."""
    try:
        with fits.open(file_path, ignore_missing_end=True) as hdul:
            if len(hdul) < 2: return None
            dat = Table(hdul[1].data)
            if 'band' not in dat.colnames: return None
            is_mag = np.array(['(' not in b for b in dat['band']])
            dat = dat[is_mag]
            return _sanitize_and_standardize(dat, {'time': 'JD', 'mag': 'm'})
    except Exception:
        return None

def parse_notes_and_limits(file_path):
    """
    Parses SN1980I photometry file, ignoring lines with upper/lower limits and handling 'null' values.
    """
    try:
        with open(file_path, 'r') as f:
            raw_lines = f.readlines()

        # Strip out comments and limits
        lines = [
            line for line in raw_lines
            if not line.strip().startswith('#') and '>' not in line and '<' not in line
        ]

        if not lines or len(lines) < 2:
            return None

        # First non-comment line should be the header
        header = lines[0].strip()
        data_lines = lines[1:]

        # Create a DataFrame from whitespace-separated values
        df = pd.read_csv(
            io.StringIO("".join(data_lines)),
            sep=r'\s{2,}',  # use two or more spaces to split columns
            header=None,
            engine='python'
        )

        # Manually assign column names based on file structure
        df.columns = ['time', 'gregorian', 'mag', 'magerr', 'band', 'reference', 'notes']

        # Clean up values
        df['mag'] = pd.to_numeric(df['mag'], errors='coerce')
        df['magerr'] = pd.to_numeric(df['magerr'].replace({'null': np.nan, 'nul': np.nan}), errors='coerce')
        df['time'] = pd.to_numeric(df['time'], errors='coerce')

        tbl = Table.from_pandas(df[['time', 'mag', 'magerr', 'band']])
        tbl['reference'] = df['reference']

        return _sanitize_and_standardize(tbl, {})
    except Exception as e:
        #print(f"Error: {e}")
        return None

def parse_csv(file_path):
    """
    Parses a hand-converted CSV version of SN1989M data from IAUC-style notes.
    Assumes columns: JD, Mag, Band, Reference
    """
    try:
        df = pd.read_csv(file_path)

        # Basic cleanup
        df = df.rename(columns={
            'Julian Date': 'time',
            'Gregorian Day': 'day',
            'Magnitude': 'mag',
            'Band': 'band',
            'Ref': 'reference',
            'Magerr': 'magerr'
        })

        df['mag'] = pd.to_numeric(df['mag'], errors='coerce')
        df['time'] = pd.to_numeric(df['time'], errors='coerce')
        df['magerr'] = pd.to_numeric(df['magerr'], errors='coerce')
        df['band'] = df['band'].astype(str).str.strip().str.strip("'").str.strip('"')

        # Leave band values like "B", "V", etc. — do not map here
        tbl = Table.from_pandas(df[['time', 'mag', 'magerr', 'band', 'reference']])
        return _sanitize_and_standardize(tbl, {})

    except Exception as e:
        #print(f"Failed to parse SN1989M CSV: {e}")
        return None



# The full list of parsers. New functions are added to the top to be tried first.
ALL_PARSERS = [
    parse_csv, 
    parse_txt,        
    parse_notes_and_limits,
    parse_fit_vizier, 
]