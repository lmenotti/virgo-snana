# This file contains the configuration for processing supernova photometry data.
# For each supernova, specify the list of files to be processed and the magnitude system.

SUPERNOVAE = {
    "SN1939A": {
        "files": ["photometry_source1_data1.fit"],
        "mag_system": "Vega"
    },
    "SN1957B": {
        "files": ["photometry_source1_data1.fit"],
        "mag_system": "Vega"
    },
    "SN1960F": {
        "files": ["photometry_source1_data1.fit","photometry_source2_data1.csv"],
        "mag_system": "Vega"
    },
    "SN1960R": {
        "files": ["photometry_source1_data1.fit"],
        "mag_system": "Vega"
    },
    "SN1961H": {
        "files": ["photometry_source1_data1.fit"],
        "mag_system": "Vega"
    },
    "SN1963I": {
        "files": ["photometry_source1_data1.csv"],
        "mag_system": "Vega"
    },
    "SN1965I": {
        "files": ["photometry_source1_data1.fit","photometry_source2_data2.csv"],
        "mag_system": "Vega"
    },
    "SN1971G": {
        "files": ["photometry_source1_data1.fit","photometry_source1_data2.csv"],
        "mag_system": "Vega"
    },
    "SN1980I": {
        "files": ["photometry_source1_data1.txt"],
        "mag_system": "Vega"
    },
    "SN1981B": {
        "files": ["photometry_source1_data1.fit"],
        "mag_system": "Vega"
    },
    "SN1983G": {
        "files": ["photometry_source1_data1.fit","photometry_source1_data2.csv","photometry_source2_data1.csv"],
        "mag_system": "Vega"
    },
    "SN1984A": {
        "files": ["photometry_source1_data1.fit","photometry_source1_data2.csv","photometry_source2_data1.txt"],
        "mag_system": "Vega"
    },
    "SN1985B": {
        "files": ["photometry_source1_data1.csv","photometry_source2_data1.txt"],
        "mag_system": "Vega"
    },
    "SN1989M": {
        "files": ["photometry_source1_data1.csv","photometry_source2_data1.csv"],
        "mag_system": "Vega"
    },
    "SN1990N": {
        "files": ["photometry_source1_data1.txt"],
        "mag_system": "Vega"
    },
    "SN1991T": {
        "files": ["photometry_source1_data1.txt"],
        "mag_system": "Vega"
    },
    "SN1991bg": {
        "files": ["photometry_source1_data1.txt"],
        "mag_system": "Vega"
    },
    "SN1992P": {
        "files": ["photometry_source1_data1.txt"],
        "mag_system": "Vega"
    },
    "SN1994D": {
        "files": ["photometry_source1_data1.txt"],
        "mag_system": "Vega"
    },
    "SN1999cl": {
        "files": ["photometry_source1_data1.txt"],
        "mag_system": "Vega"
    }
    # Future supernovae can be added here following the same format.
}