import os
import numpy as np
import pandas as pd
import sncosmo
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.simbad import Simbad
from config import SUPERNOVAE
from parsers import ALL_PARSERS

# Maps various band names found in the files to standard sncosmo registry names.
BAND_MAP = {
    'U': 'bessellux', 'B': 'bessellb', 'V': 'bessellv', 'R': 'bessellr', 'I': 'besselli',
    'pg': 'standard::b',
    'pv': 'standard::v',
    'm_v': 'bessellv',
    'm_pg': 'standard::b',
    'B_max': 'bessellb',
    'blue': 'bessellb',
    'red': 'bessellr',
    "'blue'": 'bessellb',
    "'red'": 'bessellr',
    'UNKNOWN': 'standard::u', 
    'C': 'standard::b',
}

# Hard-coded list of known bandpasses.
KNOWN_BANDS = [
    'skymapperr', 'lsstu', 'csphs', 'cspi', 'ps1::z', 'sdssu', 'ps1::i', 'f250m', 'sdss::g', 'keplercam::v', 'f1280w', 'uvf606w', 'f1500w', 'gotob', 'sdssz', 'f162m', 'f210m', 'bessellb', 'nicmos2::f110w', '2massks', 'f850lp', '2massj', 'uvf775w', 'ps1::w', 'f560w', 'gaia::grp', 'f140w', 'ztf::g', 'f360m', 'uvot::uvw2', 'f098m', 'lsstz', 'f127m', 'f182m', 'atlasc', 'f090w', 'ztfr', 'nicf160w', 'sdss::i', 'hsc::g', 'hsc::r', 'megacam6::z', 'f763m', 'f160w', 'f070w', 'hsc::r2', 'galex::fuv', 'f430m', 'hsc::z', 'cspjs', 'f438w', 'sdss::r', 'uvot::u', 'f2100w', 'f689m', 'nicf110w', 'f1800w', 'f350lp', 'cspu', 'acswf::f775w', 'cspv3009', 'standard::v', 'skymapperz', 'f335m', 'f2550w', 'keplercam::us', 'f625w', 'f184', 'f218w', 'acswf::f606w', 'lsstg', 'desr', 'gotog', 'swope2::h', 'f356w', 'swope2::v1', 'hsc::y', 'swope2::i', 'f480m', 'swope2::g', 'f606w', 'f105w', 'swope2::v2', 'ztfg', 'sdssi', 'besselli', 'megacam6::i', 'desg', 'f225w', 'f2300c', 'ps1::open', 'f555w', 'f140m', 'cspb', 'sdss::u', 'ztf::r', 'acswf::f850lp', 'uvf555w', 'ztf::i', '4shooter2::us', 'f845m', 'f770w', 'f153m', 'cspk', 'lssti', 'desi', 'f1550c', 'f129', '2massh', 'cspr', 'f390w', 'f125w', 'csphd', 'swope2::r', 'nicmos2::f160w', 'sdss::z', 'f300x', 'bessellr', 'skymapperg', 'f213', 'galex::nuv', 'sdssg', 'cspyd', 'swope2::u', 'cspg', 'desz', '4shooter2::r', 'f336w', 'gaia::grvs', 'f150w', 'megacam6::i2', 'f300m', 'f087', 'ps1::r', 'hsc::i', 'uvot::uvm2', 'lssty', 'swope2::j', '4shooter2::b', 'gaia::g', 'standard::r', 'megacam6::g', 'f275w', 'ztfi', 'swope2::b', 'cspv3014', 'standard::u', 'uvot::white', 'desu', 'uvot::b', 'atlaso', 'f110w', 'keplercam::b', 'f146', 'cspys', 'f435w', 'f115w', 'sdssr', 'gotol', 'tess', 'keplercam::i', 'f460m', 'cspv9844', 'f1065c', 'f1140c', 'desy', '4shooter2::i', 'f444w', 'hsc::i2', 'skymapperi', 'ps1::y', 'f200w', 'standard::i', 'kepler', 'f1130w', '4shooter2::v', 'uvf850lp', 'standard::b', 'f139m', 'f062', 'megacam6::r', 'f158', 'ps1::g', 'swope2::y', 'uvf814w', 'f1000w', 'f475w', 'f106', 'uvot::v', 'f775w', 'uvf625w', 'gotor', 'lsstr', 'keplercam::r', 'gaia::gbp', 'uvf475w', 'skymapperu', 'cspjd', 'f410m', 'bessellv', 'bessellux', 'uvot::uvw1', 'swope2::v', 'f277w'
]

def get_sn_metadata(sn_name):
    """Queries Simbad for supernova metadata."""
    try:
        from astroquery.ipac.irsa.irsa_dust import IrsaDust
        simbad = Simbad()
        simbad.add_votable_fields('ra', 'dec', 'rvz_redshift')
        query_result = simbad.query_object(sn_name)
        if query_result is None: return None, None, None, 0.0
        ra, dec, redshift = query_result['ra'][0], query_result['dec'][0], query_result['rvz_redshift'][0]
        if np.ma.is_masked(redshift): redshift = 0.0
        coords = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')
        try:
            mwebv = IrsaDust.get_query_table(coords, section='ebv')['E(B-V) IRSA'][0]
        except Exception: mwebv = 0.0
        return ra, dec, redshift, mwebv
    except Exception: return None, None, None, 0.0

def process_supernova(sn_name, sn_info):
    """Processes all photometry files for a single supernova."""
    print(f"--- Processing {sn_name} ---")
    photometry_tables = []
    base_dir = os.path.join("raw_virgo_data", sn_name, "Photometry")

    for file_name in sn_info.get("files", []):
        file_path = os.path.join(base_dir, file_name)
        if not os.path.exists(file_path): continue
        print(f"  Reading file: {file_name}")
        for parser_func in ALL_PARSERS:
            data = parser_func(file_path)
            if data is not None and len(data) > 0:
                print(f"    Successfully parsed with: {parser_func.__name__}")
                photometry_tables.append(data)
                break
        else:
            print(f"  WARNING: Could not parse file {file_name} with any available parser.")

    if not photometry_tables:
        print(f"--- No data processed for {sn_name}, skipping SNANA file creation. ---\n")
        return

    full_table = vstack(photometry_tables, join_type='outer')
    df = full_table.to_pandas().drop_duplicates(subset=['time', 'mag'], keep='first')
    
    df['band'] = df['band'].str.strip()
    df = df[~df['band'].str.contains(r'\(', na=False)]
    print("Unique bands before mapping:", df['band'].unique())
    unmapped = df[~df['band'].isin(BAND_MAP.keys())]
    print("Unmapped bands:", unmapped['band'].unique())
    df['band'] = df['band'].map(BAND_MAP).fillna(df['band'])
    
    df = df[df['band'].isin(KNOWN_BANDS)]
    
    if df.empty:
        print(f"  WARNING: No usable data remains for {sn_name} after parsing and band mapping.\n")
        return
        
    full_table = Table.from_pandas(df)
    print(f"  Found {len(full_table)} unique, usable photometric points for {sn_name}.")

    ra, dec, redshift, mwebv = get_sn_metadata(sn_name)
    if ra is None: return

    metadata = {
        "SURVEY": "VIRGO_PROJECT", "SNID": sn_name, "RA": f"{ra:.8f}", "DEC": f"{dec:.8f}",
        "MWEBV": f"{mwebv:.4f}", "REDSHIFT_HELIO": f"{redshift:.6f}",
        "FILTERS": " ".join(sorted(np.unique(full_table['band'])))
    }

    zp_system_str = sn_info.get("mag_system", "Vega").lower()
    magsys = sncosmo.get_magsystem(zp_system_str)

    flux_list, fluxerr_list = [], []
    for row in full_table:
        flux = magsys.band_mag_to_flux(row['mag'], row['band'])
        flux_list.append(flux)
        fluxerr = flux * 0.4 * np.log(10) * row['magerr'] if np.isfinite(row['magerr']) else -999.0
        fluxerr_list.append(fluxerr)

    final_table = Table({
        'time': np.atleast_1d(full_table['time']), 'band': np.atleast_1d(full_table['band']),
        'flux': np.atleast_1d(flux_list), 'fluxerr': np.atleast_1d(fluxerr_list),
        'mag': np.atleast_1d(full_table['mag']), 'magerr': np.nan_to_num(np.atleast_1d(full_table['magerr']), nan=-999.0),
        'zp': np.full(len(full_table), 25.0), 'zpsys': np.full(len(full_table), zp_system_str, dtype='<U10'),
    })
    final_table.meta = metadata

    output_dir = os.path.join("snana_virgo_data", sn_name, "Photometry")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{sn_name}.photometry.snana.dat")

    print(f"  Writing standard SNANA file to: {output_file}")
    sncosmo.write_lc(final_table, output_file, format='snana')
    print(f"--- Finished {sn_name} --- \n")

if __name__ == "__main__":
    for sn_name, sn_info in SUPERNOVAE.items():
        process_supernova(sn_name, sn_info)