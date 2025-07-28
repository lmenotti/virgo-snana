import io
import re
import numpy as np
import pandas as pd
from astropy.io import fits, ascii
from astropy.table import Table
from astropy.time import Time

def _sanitize_and_standardize(tbl, column_map):
    """A robust helper function to rename columns, ensure numeric types, and filter bad rows."""
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

# --- NEW PARSERS ADDED FOR THIS ITERATION ---

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

        # Correct column renaming
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


# def parse_SN1965I_specific(file_path):
#     """
#     Specifically for SN1965I/photometry_source2_data2 (copy).txt.
#     This format is fixed-width with a complex header.
#     """
#     try:
#         with open(file_path, 'r') as f:
#             header = f.readline().lower()
#         if 'indmag and band' not in header or 'reference text' not in header:
#             return None
            
#         df = pd.read_fwf(file_path)
#         tbl = Table.from_pandas(df)
#         return _sanitize_and_standardize(tbl, {
#             'time': 'Julian Date', 'mag': 'Magnitude', 'band': 'Indmag and Band'
#         })
#     except Exception:
#         return None
# def parse_fit_vizier_lenient2(file_path):
#     """Parses VizieR FITS files with common quirks and non-mag data."""
#     try:
#         with fits.open(file_path, ignore_missing_end=True) as hdul:
#             if len(hdul) < 2:
#                 print("FITS file has no extension with data.")
#                 return None

#             dat = Table(hdul[1].data)

#             # Ensure 'band' is present
#             if 'band' not in dat.colnames or 'm' not in dat.colnames or 'JD' not in dat.colnames:
#                 print("Required columns missing:", dat.colnames)
#                 return None

#             # Convert FITS byte strings to native str
#             dat['band'] = dat['band'].astype(str)

#             # Filter out color indices like (B-V), (U-B), etc.
#             is_mag = ~dat['band'].str.contains(r'\(')
#             dat = dat[is_mag]

#             return _sanitize_and_standardize(dat, {'time': 'JD', 'mag': 'm'})

#     except Exception as e:
#         print(f"Failed to parse FITS file: {e}")
#         return None
    
# --- Baseline Parsers from the 10-Failure Run (Unchanged) ---

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

# def parse_simple_tab_separated(file_path):
#     """
#     Parses simple tab-separated files with 'Indmag and Band' column.
#     """
#     try:
#         with open(file_path, 'r') as f:
#             header = f.readline().lower()
#         if 'indmag and band' not in header or 'uncertainty' in header or 'reference' in header:
#             return None
            
#         tbl = Table.read(file_path, format='ascii.tab', header_start=0, data_start=1)
#         return _sanitize_and_standardize(tbl, {'time': 'Julian Date', 'mag': 'Magnitude', 'band': 'Indmag and Band'})
#     except Exception:
#         return None

# def parse_tab_with_uncertainty(file_path):
#     """
#     Parses tab-separated files with an 'Uncertainty' column.
#     """
#     try:
#         with open(file_path, 'r') as f:
#             header = f.readline().lower()
#         if 'uncertainty' not in header:
#             return None
        
#         tbl = Table.read(file_path, format='ascii.tab', header_start=0, data_start=1)
#         return _sanitize_and_standardize(tbl, {'time': 'Julian Date', 'mag': 'Magnitude', 'magerr': 'Uncertainty', 'band': 'Indmag and Band'})
#     except Exception:
#         return None

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



# def parse_tab_no_band(file_path):
#     """
#     Parses tab-separated files with a full reference column but NO band column.
#     """
#     try:
#         with open(file_path, 'r') as f:
#             header = f.readline().lower()
#         if 'reference text' not in header or 'band' in header:
#             return None

#         tbl = Table.read(file_path, format='ascii.tab', header_start=0, data_start=1)
#         return _sanitize_and_standardize(tbl, {
#             'time': 'Julian Date', 'mag': 'Magnitude', 'reference': 'Reference Text'
#         })
#     except Exception:
#         return None

# def parse_fwf_full_reference(file_path):
#     """
#     Parses fixed-width files with 'Indmag and Band' and a full reference column.
#     """
#     try:
#         with open(file_path, 'r') as f:
#             header = f.readline().lower()
#         if 'indmag and band' not in header or 'reference text' not in header:
#             return None

#         df = pd.read_fwf(file_path, comment='#')
#         tbl = Table.from_pandas(df)
#         return _sanitize_and_standardize(tbl, {
#             'time': 'Julian Date', 'mag': 'Magnitude', 'band': 'Indmag and Band'
#         })
#     except Exception:
#         return None

from astropy.time import Time
from datetime import datetime

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
        print(f"Failed to parse SN1989M CSV: {e}")
        return None

# def parse_SN1989M_specific_txt(file_path):
#     """
#     Specifically for SN1989M/photometry_source1_data1.txt.
#     Splits 'Indmag and Band' into band and magerr.
#     """
#     try:
#         with open(file_path, 'r') as f:
#             header = f.readline()
#         if not all(x in header for x in ['Julian Date', 'Gregorian Day', 'Magnitude', 'Indmag and Band']):
#             return None

#         df = pd.read_csv(file_path, sep='\t')

#         # Split "Indmag and Band Magerr" into two columns manually
#         if 'Indmag and Band' in df.columns:
#             split_cols = df['Indmag and Band'].astype(str).str.strip().str.split(r'\s+', expand=True)
#             if split_cols.shape[1] == 2:
#                 df['band'] = split_cols[0].str.strip().str.strip("'").str.strip('"')
#                 df['magerr'] = pd.to_numeric(split_cols[1], errors='coerce')
#             else:
#                 df['band'] = 'UNKNOWN'
#                 df['magerr'] = np.nan
#         else:
#             df['band'] = 'UNKNOWN'
#             df['magerr'] = np.nan

#         df = df.rename(columns={'Julian Date': 'time', 'Magnitude': 'mag'})

#         tbl = Table.from_pandas(df[['time', 'mag', 'magerr', 'band']])
#         tbl['reference'] = 'N/A'
#         return _sanitize_and_standardize(tbl, {})

#     except Exception as e:
#         print(f"Failed to parse SN1989M TXT: {e}")
#         return None


# The full list of parsers. New functions are added to the top to be tried first.
ALL_PARSERS = [
    parse_csv, 
    parse_txt,        
    #parse_SN1965I_specific,        
    parse_notes_and_limits,
    #parse_tab_with_uncertainty,
    #parse_tab_no_band,
    #parse_fwf_full_reference,
    #parse_simple_tab_separated,
    parse_fit_vizier, 
    #parse_fit_vizier_lenient2,
    #parse_SN1989M_specific_txt,

]