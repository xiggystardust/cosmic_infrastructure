#!/usr/bin/env python3
""".

Initially prompted Claude Sonnet 4.5 on 20 July 2026 using HopGPT
platform.  Sarah then debugged and modified things by hand. Original
prompt at the end of the code.

Purpose: Cross-match astronomical sources with multiple catalogs.

This script reads a galaxy sample file and cross-matches each source with:
- BIGMAC catalog (binary AGN candidates)
- BOBCat (binary black hole catalog)
- RFC (Radio Fundamental Catalog)
- VLASS (VLA Sky Survey)
- FIRST (Faint Images of the Radio Sky at Twenty-cm)
- NVSS (NRAO VLA Sky Survey)

Output is a CSV file with all matched information.
"""

import pandas as pd
import numpy as np
from astroquery.simbad import Simbad
from astroquery.ipac.ned import Ned
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.vizier import Vizier
from astroquery.skyview import SkyView
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# Configuration
# ============================================================================

BASE = "/Users/sbs/Desktop/current_projects/horlaville/full_cross_match/source/"
INPUT_FILE = BASE + "Horlaville_apjae6662t2_mrt_no_dupes.txt"
RADEC_FILE = BASE + "nameradec_fullsample.txt"
BIGMAC_FILE = BASE + "BigMAC_maintable_DR0p9.csv"
#BOBCAT_FILE = BASE + "bobcat_07_20_2026.list"
BOBCAT_FILE = BASE + "bobcat_pos_z_uniq.csv"
RFC_FILE = BASE + "rfc_2026b_cat.txt"
OUTPUT_FILE = BASE + "new_cross_matched_galaxies.csv"
P2OUTPUT_FILE = BASE + "TARGET_FULL_LDA+2new_cross_matched_galaxies.csv"
N2OUTPUT_FILE = BASE + "LDA-1new_cross_matched_galaxies.csv"

# NEW: Nyland files
NYLAND_MASSIVE_FILE = BASE + "eden_H26_crossmatch_sorted.txt"
NYLAND_ATLAS_FILE = BASE + "nyland+2016-atlas3D.txt"

# NEW: extra LDA>2 outputs
KU_ONLY_LDA2_FILE = BASE + "TARGET_Ku_only_LDA+2new_cross_matched_galaxies.csv"
C_KU_LDA2_FILE = BASE + "TARGET_C+Ku_LDA+2new_cross_matched_galaxies.csv"


MATCH_RADIUS = 10.0 * u.arcsec  # Cross-match radius

# ============================================================================
# Function Definitions
# ============================================================================

def read_input_file(filename):
    """
    Read the input galaxy catalog file.
    
    Parameters:
    -----------
    filename : str
        Path to input file
    
    Returns:
    --------
    pandas.DataFrame
        DataFrame with galaxy information
    n_lda_positive how many have LDA >0
    """
    data = []
    
    n_lda_positive = 0
    n_lda_plusone = 0
    n_lda_negone = 0

    with open(filename, 'r') as f:
        for line in f:
            # Skip comments and empty lines
            if line.startswith('#') or len(line.strip()) == 0:
                continue
            
            # Parse fixed-width format
            try:
                name = line[0:13].strip()
                rank = int(line[14:17].strip())
                distance = float(line[18:23].strip())
                survey = line[24:31].strip()
                logMBH = float(line[32:37].strip())
                logh0 = float(line[38:44].strip())
                lda = float(line[45:50].strip())
                logh0norm = float(line[51:55].strip())
                ldanorm = float(line[56:60].strip())
                score = float(line[61:65].strip())
                profile = line[66:67].strip() if len(line) > 66 else ''

                if lda>0:
                    n_lda_positive += 1
                if lda>2:
                    n_lda_plusone += 1
                if lda<-2:
                    n_lda_negone += 1
                
                data.append({
                    'Name': name,
                    'Rank': rank,
                    'Distance_Mpc': distance,
                    'Survey': survey,
                    'logMBH': logMBH,
                    'logh0': logh0,
                    'LDA': lda,
                    'logh0norm': logh0norm,
                    'LDAnorm': ldanorm,
                    'Score': score,
                    'Profile': profile
                })
            except Exception as e:
                print(f"Warning: Could not parse line: {line.strip()}")
                print(f"Error: {e}")
                exit()
    
    return pd.DataFrame(data),n_lda_positive,n_lda_plusone,n_lda_negone


def read_name_set_from_horliville_format(filename):
    """
    Read a Horlaville-format fixed-width file and return a set of galaxy names.
    Only the Name column is used for cross-matching.
    """
    names = set()

    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or len(line.strip()) == 0:
                continue

            try:
                name = line[0:13].strip()
                if name:
                    names.add(name.upper())
            except Exception:
                continue

    return names

def read_nyland_atlas_names(filename):
    """
    Read NylandATLAS table and return a set of galaxy names.
    The file is whitespace-delimited; the first column is the galaxy name.
    """
    names = set()

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or len(line) == 0:
                continue

            parts = line.split()
            if len(parts) == 0:
                continue

            name = parts[0].strip()
            if name:
                names.add(name.upper())

    return names

def nyland_crossmatch_flag(galaxy_name, nyland_name_set):
    """
    Return 'Y' if galaxy_name is present in the provided Nyland name set, else 'N'.
    """
    return 'Y' if galaxy_name.upper() in nyland_name_set else 'N'


# I dropped this in from my other code but it doesn't differentiate NGC7436 and NGC7436B positions
#def coord_finder(object_name):
#    '''.
#    
#    Query CDS SESAME using object_name to return its J2000 ra and dec,
#    if object is found.
#
#    Inputs:
#        object_name = string of the name of the object the coordinates
#                      are needed for. Name must be resolvable by
#                      standard astro name query.
#
#    Outputs:
#        ra = J2000 right ascension, units = hms
#        dec = J2000 declination, units = dms
#
#    '''
#
#    # Basic user error check: name should be a string.
#    if not isinstance(object_name, str):
#        raise TypeError("coord_finder: Name must be a string.")
#
#    # Search for the coordinates of the object given the name. 
#    try:
#        coords = SkyCoord.from_name(object_name).to_string('hmsdms').split()
#    except Exception as err:
#        raise RuntimeError(f"Object {object_name} is not recognized by SESAME.")
#    
#    # Assign the ra and dec variables.
#    ra = coords[0]
#    dec = coords[1]
#    
#    # Return an array of the ra and dec. They are hms/dms strings.
#    return ra, dec


# This doesn't differentiate NGC7436 and NGC7436B positions
#def get_galaxy_coordinates(galaxy_name):
#    """
#    Query Simbad for galaxy coordinates.
#    
#    Parameters:
#    -----------
#    galaxy_name : str
#        Name of the galaxy
#    
#    Returns:
#    --------
#    SkyCoord or None
#        Coordinate object or None if not found
#    """
#    try:
#        result = Simbad.query_object(galaxy_name)
#        if result is not None:
#            ra = result['ra'][0]
#            dec = result['dec'][0]
#            coord = SkyCoord(ra, dec, unit=(u.deg), frame='icrs')
#            #print(f"found {coord} from {ra} and {dec}.")
#            return coord
#    except Exception as e:
#        print(f"Could not resolve {galaxy_name}: {e}")
#    return None
#def coord_finder(galaxy_name):
#    """
#    Query Ned for galaxy coordinates.
#    Sarah adapted to NED expected outputs.
#    
#    Parameters:
#    -----------
#    galaxy_name : str
#        Name of the galaxy
#    
#    Returns:
#    --------
#    SkyCoord or None
#        Coordinate object or None if not found
#    """
#    try:
#        result = Ned.query_object(galaxy_name)
#        if result is not None:
#            ra = result['RA'][0]
#            dec = result['DEC'][0]
#            coord = SkyCoord(ra, dec, unit=(u.deg), frame='icrs')
#            #print(f"found {coord} from {ra} and {dec}.")
#            return coord
#    except Exception as e:
#        print(f"Could not resolve {galaxy_name}: {e}")
#    return None

def coord_finder(galaxy_name):
    coord = "---"
    with open(RADEC_FILE, 'r') as f:
        for line in f:
            # Skip comments and empty lines
            if line.startswith('#') or len(line.strip()) == 0:
                continue

            # Remove frass
            line.strip()
            name,ra,dec = line.split()
            if name == galaxy_name:
                coord = SkyCoord(ra, dec, unit=(u.deg), frame='icrs')
                return coord
    if coord == "---":
        raise Exception(f"Could not resolve coordinates for {galaxy_name}")


def read_bigmac_catalog(filename):
    """
    Read the BIGMAC catalog.
    
    Parameters:
    -----------
    filename : str
        Path to BIGMAC CSV file
    
    Returns:
    --------
    pandas.DataFrame
        BIGMAC catalog data
    """
    try:
        #df = pd.read_csv(filename, comment='#')
        df = pd.read_csv(filename, header=1)
        return df
    except Exception as e:
        print(f"Error reading BIGMAC catalog: {e}")
        exit()
        #return pd.DataFrame()

        
def read_bobcat_catalog(filename):
    """
    Read the BOBCAT catalog containing RA, Dec, z.
    
    Parameters:
    -----------
    filename : str
        Path to BOBCAT csv file.
    
    Returns:
    --------
    pandas.DataFrame
        BOBCAT catalog data
    """
    try:
        df = pd.read_csv(filename, comment='#')
        return df
    except Exception as e:
        print(f"Error reading BOBCAT catalog: {e}")
        exit()
        #return pd.DataFrame()


def old_read_bobcat(filename):
    """
    Read the BOBCat list.
    
    Parameters:
    -----------
    filename : str
        Path to BOBCat file
    
    Returns:
    --------
    list
        List of source names in BOBCat
    """
    try:
        with open(filename, 'r') as f:
            sources = [line.strip() for line in f if not line.startswith('#') and len(line.strip()) > 0]
        return sources
    except Exception as e:
        print(f"Error reading BOBCat: {e}")
        exit()

def read_rfc_catalog(filename):
    """
    Read the Radio Fundamental Catalog.
    
    Parameters:
    -----------
    filename : str
        Path to RFC file
    
    Returns:
    --------
    pandas.DataFrame
        RFC catalog with coordinates and flux densities
    """
    data = []
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('#') or len(line.strip()) == 0:
                    continue
                
                try:
                    name = line[0:14].strip()
                    comnam = line[16:24].strip()
                    rah = int(line[26:28].strip())
                    ram = int(line[29:31].strip())
                    ras = float(line[32:41].strip())
                    dec_sign = line[42:43].strip()
                    decd = int(line[43:45].strip())
                    decm = int(line[46:48].strip())
                    decs = float(line[49:57].strip())
                    
                    # Parse flux densities (handle missing values)
                    def parse_flux(s):
                        try:
                            val = float(s.strip())
                            return val if val > -9.0 else np.nan
                        except:
                            return np.nan
                    
                    fs_s = parse_flux(line[104:110]) if len(line) > 110 else np.nan
                    fl_s = parse_flux(line[118:124]) if len(line) > 124 else np.nan
                    fs_c = parse_flux(line[126:132]) if len(line) > 132 else np.nan
                    fl_c = parse_flux(line[140:146]) if len(line) > 146 else np.nan
                    fs_x = parse_flux(line[148:154]) if len(line) > 154 else np.nan
                    fl_x = parse_flux(line[162:168]) if len(line) > 168 else np.nan
                    fs_u = parse_flux(line[170:176]) if len(line) > 176 else np.nan
                    fl_u = parse_flux(line[184:190]) if len(line) > 190 else np.nan
                    fs_k = parse_flux(line[192:198]) if len(line) > 198 else np.nan
                    fl_k = parse_flux(line[206:212]) if len(line) > 212 else np.nan
                    
                    # Create coordinate
                    ra_str = f"{rah}h{ram}m{ras}s"
                    dec_str = f"{dec_sign}{decd}d{decm}m{decs}s"
                    coord = SkyCoord(ra_str, dec_str, frame='icrs')
                    
                    data.append({
                        'RFC_Name': name,
                        'RFC_CommonName': comnam,
                        'RFC_RA': coord.ra.deg,
                        'RFC_Dec': coord.dec.deg,
                        'RFC_FsS': fs_s,
                        'RFC_FlS': fl_s,
                        'RFC_FsC': fs_c,
                        'RFC_FlC': fl_c,
                        'RFC_FsX': fs_x,
                        'RFC_FlX': fl_x,
                        'RFC_FsU': fs_u,
                        'RFC_FlU': fl_u,
                        'RFC_FsK': fs_k,
                        'RFC_FlK': fl_k
                    })
                except Exception as e:
                    continue
        
        return pd.DataFrame(data)
    except Exception as e:
        print(f"Error reading RFC catalog: {e}")
        exit()
        #return pd.DataFrame()


def match_bigmac(galaxy_name, coord, bigmac_df):
    """
    Match galaxy to BIGMAC catalog.
    """
    result = {
        'BIGMAC_detected': 'N',
        'BIGMAC_PrimaryType': '',
        'BIGMAC_SecondaryType': '',
        'BIGMAC_TertiaryType': ''
    }
    
    if bigmac_df.empty:
        return result

    for col in ['Name1', 'Name2']:
        if col in bigmac_df.columns:
            mask = bigmac_df[col].str.contains(galaxy_name, case=False, na=False)
            if mask.any():
                print(f"For BIGMAC I found name match for {galaxy_name}")
                match = bigmac_df[mask].iloc[0]
                result['BIGMAC_detected'] = 'Y'
                result['BIGMAC_PrimaryType'] = match.get('Primary System Type', '')
                result['BIGMAC_SecondaryType'] = match.get('Secondary System Type', '')
                result['BIGMAC_TertiaryType'] = match.get('Tertiary System Type', '')
                return result

    if coord is not None and 'RA1' in bigmac_df.columns and 'Dec1' in bigmac_df.columns:
        try:
            bigmac_coords = SkyCoord(bigmac_df['RA1'], bigmac_df['Dec1'], 
                                     unit=(u.hourangle,u.deg), frame='icrs')
            sep = coord.separation(bigmac_coords)
            min_idx = sep.argmin()
            if sep[min_idx] < MATCH_RADIUS:
                print(f"For BIGMAC I found coordinate match for {galaxy_name}")
                print("*"*80)
                match = bigmac_df.iloc[min_idx]
                result['BIGMAC_detected'] = 'Y'
                result['BIGMAC_PrimaryType'] = match.get('Primary System Type', '')
                result['BIGMAC_SecondaryType'] = match.get('Secondary System Type', '')
                result['BIGMAC_TertiaryType'] = match.get('Tertiary System Type', '')
        except:
            raise Exception
    
    return result


def match_bobcat(galaxy_name, coord, bobcat_df):
    """
    Match galaxy to BOBCAT catalog.
    """

    if bobcat_df.empty:
        return result
    
    if 'name' in bobcat_df.columns:
        mask = bobcat_df['name'].str.contains(galaxy_name, case=False, na=False)
        if mask.any():
            print(f"DEBUG: I found name match for {galaxy_name}")
            return 'Y'
    
    if coord is not None and 'ra' in bobcat_df.columns and 'dec' in bobcat_df.columns:
        bobcat_coords = SkyCoord(bobcat_df['ra'], bobcat_df['dec'],
                                    unit=u.deg, frame='icrs')
        sep = coord.separation(bobcat_coords)
        min_idx = sep.argmin()
        if sep[min_idx] < MATCH_RADIUS:
            match = bobcat_df.iloc[min_idx]
            print(f"DEBUG: I found a position match for {galaxy_name}")
            return 'Y'

    return 'N'

def old_check_bobcat(galaxy_name, bobcat_list):
    for source in bobcat_list:
        if galaxy_name.upper() in source.upper() or source.upper() in galaxy_name.upper():
            return 'Y'
    return 'N'


def match_rfc(coord, rfc_df):
    """
    Match coordinates to RFC catalog.
    """
    result = {
        'RFC_Match': 'N',
        'RFC_FsS': np.nan,
        'RFC_FlS': np.nan,
        'RFC_FsC': np.nan,
        'RFC_FlC': np.nan,
        'RFC_FsX': np.nan,
        'RFC_FlX': np.nan,
        'RFC_FsU': np.nan,
        'RFC_FlU': np.nan,
        'RFC_FsK': np.nan,
        'RFC_FlK': np.nan
    }
    
    if coord is None or rfc_df.empty:
        return result
    
    try:
        rfc_coords = SkyCoord(rfc_df['RFC_RA'], rfc_df['RFC_Dec'], 
                             unit=u.deg, frame='icrs')
        sep = coord.separation(rfc_coords)
        min_idx = sep.argmin()
        
        if sep[min_idx] < MATCH_RADIUS:
            match = rfc_df.iloc[min_idx]
            result['RFC_Match'] = 'Y'
            for key in ['RFC_FsS', 'RFC_FlS', 'RFC_FsC', 'RFC_FlC', 
                       'RFC_FsX', 'RFC_FlX', 'RFC_FsU', 'RFC_FlU', 
                       'RFC_FsK', 'RFC_FlK']:
                result[key] = match[key]
    except Exception as e:
        print(f"Error matching RFC: {e}")
        exit()
    
    return result


def query_vlass(coord):
    if coord is None:
        return np.nan
    
    try:
        v = Vizier(columns=["**"], row_limit=1)
        v.ROW_LIMIT = 1
        result = v.query_region(coord, radius=MATCH_RADIUS, 
                               catalog="J/ApJS/255/30/catalog")
        
        if len(result) > 0 and 'Ftot' in result[0].colnames:
            return result[0]['Ftot'][0]
    except Exception as e:
        print(f"VLASS query error: {e}")
        exit()
    
    return np.nan


def query_first(coord):
    if coord is None:
        return 'N', np.nan
    
    try:
        v = Vizier(columns=["**"], row_limit=1)
        result = v.query_region(coord, radius=MATCH_RADIUS, 
                               catalog="VIII/92/first")
        
        if len(result) > 0:
            flux = result[0]['Fpeak'][0] if 'Fpeak' in result[0].colnames else np.nan
            return 'Y', flux
        else:
            if coord.dec.deg > -10:
                return 'Y', -999.0
            else:
                return 'N', np.nan
    except Exception as e:
        print(f"FIRST query error: {e}")
        exit()
    
    return 'N', np.nan


def query_nvss(coord):
    if coord is None:
        return 'N', np.nan
    
    try:
        v = Vizier(columns=["**"], row_limit=1)
        result = v.query_region(coord, radius=MATCH_RADIUS, 
                               catalog="VIII/65/nvss")
        
        if len(result) > 0:
            flux = result[0]['S1.4'][0] if 'S1.4' in result[0].colnames else np.nan
            return 'Y', flux
        else:
            if coord.dec.deg > -40:
                return 'Y', -999.0
            else:
                return 'N', np.nan
    except Exception as e:
        print(f"NVSS query error: {e}")
        exit()
    
    return 'N', np.nan




def drop_source(df, df_id, column_name, value):
    """
    Drop a galaxy in a given df, if its column_name equals the value.
    Report it being dropped from dataframe df_id

    return true if it's there and false if it's not.
    """
    if column_name not in df.columns:
        print("*"*100)
        print(f"DEBUG: column name {column_name} not present in data frame.")
        return False

    if value not in df[column_name].values:
        return False

    print(f"Dropping {column_name}: {value} from data frame {df_id}.")
    df.drop(df[df[column_name] == value].index, inplace=True)

    return True


# ============================================================================
# Main Processing
# ============================================================================

def main():
    """Main processing function."""
    
    print("=" * 80)
    print("Cross-matching Galaxy Sample with Multiple Catalogs")
    print("=" * 80)
    
    # Read input file
    print("\n1. Reading input galaxy catalog...")
    galaxies_df,n_lda_positive,n_lda_plusone,n_lda_negone = read_input_file(INPUT_FILE)
    print(f"   Found {len(galaxies_df)} galaxies")
    
    # Read reference catalogs
    print("\n2. Loading reference catalogs...")
    print("   - BIGMAC catalog...")
    bigmac_df = read_bigmac_catalog(BIGMAC_FILE)
    print(f"     Loaded {len(bigmac_df)} sources")
    

    print("   - BOBCAT catalog...")
    bobcat_df = read_bobcat_catalog(BOBCAT_FILE)
    print(f"     Loaded {len(bobcat_df)} sources")
    
    print("   - RFC catalog...")
    rfc_df = read_rfc_catalog(RFC_FILE)
    print(f"     Loaded {len(rfc_df)} sources")

    # NEW: Nyland catalogs
    print("   - NylandMASSIVE catalog...")
    nyland_massive_names = read_name_set_from_horliville_format(NYLAND_MASSIVE_FILE)
    print(f"     Loaded {len(nyland_massive_names)} source names")

    print("   - NylandATLAS catalog...")
    nyland_atlas_names = read_nyland_atlas_names(NYLAND_ATLAS_FILE)
    print(f"     Loaded {len(nyland_atlas_names)} source names")
    
    # Process each galaxy
    print("\n3. Cross-matching galaxies...")
    results = []

    # NEW: counters
    nyland_atlas_count = 0
    nyland_massive_count = 0
    
    for idx, row in galaxies_df.iterrows():
        galaxy_name = row['Name']
        print(f"\n   Processing {idx+1}/{len(galaxies_df)}: {galaxy_name}")
        
        # Get coordinates
        print(f"      Resolving coordinates...")
        coord = coord_finder(galaxy_name)
        
        
        if coord is not None:
            print(f"      RA={coord.ra.deg:.4f}, Dec={coord.dec.deg:.4f}")
        else:
            print(f"      Warning: Could not resolve coordinates")
        
        # Start building result
        result = row.to_dict()
        result['RA_deg'] = coord.ra.deg if coord else np.nan
        result['Dec_deg'] = coord.dec.deg if coord else np.nan

        # NEW: Nyland flags

        print(f"      Matching BIGMAC...")        
        result['NylandATLAS'] = nyland_crossmatch_flag(galaxy_name, nyland_atlas_names)
        result['NylandMASSIVE'] = nyland_crossmatch_flag(galaxy_name, nyland_massive_names)

        if result['NylandATLAS'] == 'Y':
            nyland_atlas_count += 1
        if result['NylandMASSIVE'] == 'Y':
            nyland_massive_count += 1
        
        # BIGMAC match
        print(f"      Matching BIGMAC...")
        bigmac_match = match_bigmac(galaxy_name, coord, bigmac_df)
        result.update(bigmac_match)
        
        ## Old BOBCat check
        #print(f"      Checking BOBCat...")
        #result['BOBCat_Match'] = check_bobcat(galaxy_name, bobcat_list)

        # New BOBCAT search (ra/dec based)
        print(f"      Matching BOBCAT...")
        result['BOBCat_Match'] = match_bobcat(galaxy_name, coord, bobcat_df)
        
        # RFC match
        print(f"      Matching RFC...")
        rfc_match = match_rfc(coord, rfc_df)
        result.update(rfc_match)
        
        # VLASS query
        print(f"      Querying VLASS...")
        result['VLASS_Flux_mJy'] = query_vlass(coord)
        
        # FIRST query
        print(f"      Querying FIRST...")
        first_footprint, first_flux = query_first(coord)
        result['FIRST_InFootprint'] = first_footprint
        result['FIRST_Flux_mJy'] = first_flux
        
        # NVSS query
        print(f"      Querying NVSS...")
        nvss_footprint, nvss_flux = query_nvss(coord)
        result['NVSS_InFootprint'] = nvss_footprint
        result['NVSS_Flux_mJy'] = nvss_flux
        
        results.append(result)

        ##DEBUG
        #if idx>3:
        #    break

    # Create output DataFrame
    print("\n4. Creating output file...")
    output_df = pd.DataFrame(results)
    
    # Reorder columns
    column_order = [
        'Name', 'Rank', 'Distance_Mpc', 'Survey', 'logMBH', 'logh0', 
        'LDA', 'logh0norm', 'LDAnorm', 'Score', 'Profile',
        'RA_deg', 'Dec_deg','BIGMAC_detected',
        'BIGMAC_PrimaryType', 'BIGMAC_SecondaryType', 'BIGMAC_TertiaryType',
        'BOBCat_Match',
        'NylandATLAS', 'NylandMASSIVE',
        'RFC_Match', 'RFC_FsS', 'RFC_FlS', 'RFC_FsC', 'RFC_FlC', 
        'RFC_FsX', 'RFC_FlX', 'RFC_FsU', 'RFC_FlU', 'RFC_FsK', 'RFC_FlK',
        'VLASS_Flux_mJy',
        'FIRST_InFootprint', 'FIRST_Flux_mJy',
        'NVSS_InFootprint', 'NVSS_Flux_mJy'
    ]

    
    output_df = output_df[column_order]
    p2output_df = output_df[output_df['LDA'] > 2].copy().reset_index(drop=True)
    n2output_df = output_df[output_df['LDA'] < -1].copy().reset_index(drop=True)

    # Remove known dual AGN and recoils
    # NGC6338 - Comerford et al. dual.
    print("NGC6338 was identified by Comerford et al as a dual AGN.")
    drop_source(output_df,"full list",'Name',galaxy_name)
    drop_source(p2output_df,"full list",'Name',galaxy_name)
    drop_source(n2output_df,"full list",'Name',galaxy_name)

    """
    # NGC4162 - Lena et al. observed in VLBI
    print("NGC4162 was identified by Lena et al as a recoiling AGN.")
    drop_source(output_df,"full list",'Name',galaxy_name)
    drop_source(p2output_df,"full list",'Name',galaxy_name)
    drop_source(n2output_df,"full list",'Name',galaxy_name)

    # NGC4278 - Lena et al. observed in VLBI
    print("NGC4278 was  identified by Lena et al as a recoiling AGN.")
    drop_source(output_df,"full list",'Name',galaxy_name)
    drop_source(p2output_df,"full list",'Name',galaxy_name)
    drop_source(n2output_df,"full list",'Name',galaxy_name)

    # NGC4486 - Lena et al. observed in VLBI
    print("NGC4486 was  identified by Lena et al as a recoiling AGN.")
    drop_source(output_df,"full list",'Name',galaxy_name)
    drop_source(p2output_df,"full list",'Name',galaxy_name)
    drop_source(n2output_df,"full list",'Name',galaxy_name)

    # NGC4636 - Lena et al. observed in VLBI
    print("NGC4636 was  identified by Lena et al as a recoiling AGN.")
    drop_source(output_df,"full list",'Name',galaxy_name)
    drop_source(p2output_df,"full list",'Name',galaxy_name)
    drop_source(n2output_df,"full list",'Name',galaxy_name)

    # NGC5846 - Lena et al. observed in VLBI
    print("NGC5846 was  identified by Lena et al as a recoiling AGN.")
    drop_source(output_df,"full list",'Name',galaxy_name)
    drop_source(p2output_df,"full list",'Name',galaxy_name)
    drop_source(n2output_df,"full list",'Name',galaxy_name)
    """
    
    # NEW: additional LDA>2 split files
    ku_only_lda2_df = p2output_df[
        (p2output_df['NylandMASSIVE'] == 'Y') | (p2output_df['NylandATLAS'] == 'Y')
    ].copy().reset_index(drop=True)

    c_ku_lda2_df = p2output_df[
        (p2output_df['NylandMASSIVE'] == 'N') & (p2output_df['NylandATLAS'] == 'N')
    ].copy().reset_index(drop=True)
    
    # Create header comments
    header_comments = """# Cross-matched Galaxy Sample with Multiple Catalogs
# 
# Column Descriptions:
# -------------------
# Name: Galaxy name from input catalog
# Rank: Total score rank from input catalog
# Distance_Mpc: Luminosity distance in Mpc
# Survey: IFU survey name (MASSIVE or ATLAS3D)
# logMBH: Log of black hole mass in solar masses
# logh0: Log of hypothetical gravitational wave strain
# LDA: Linear Discriminant Analysis score (Equation 2)
# logh0norm: Normalized log of hypothetical GW strain
# LDAnorm: Normalized LDA score
# Score: Total score (Equation 4)
# Profile: Inner light profile (C=core, P=power-law, I=intermediate)
# RA_deg: Right Ascension in degrees (J2000) from Simbad
# Dec_deg: Declination in degrees (J2000) from Simbad
# BIGMAC_detected: Y/N indicator if source is in BIGMAC
# BIGMAC_PrimaryType: Primary system type from BIGMAC catalog
# BIGMAC_SecondaryType: Secondary system type from BIGMAC catalog
# BIGMAC_TertiaryType: Tertiary system type from BIGMAC catalog
# BOBCat_Match: Y/N indicator if source is in BOBCat
# NylandATLAS: Y/N indicator if source is present in nylandATLAS.txt
# NylandMASSIVE: Y/N indicator if source is present in eden_H26_crossmatch_sorted.txt
# RFC_Match: Y/N indicator if source is in Radio Fundamental Catalog
# RFC_FsS: RFC S-band short baseline flux density in Jy
# RFC_FlS: RFC S-band long baseline flux density in Jy
# RFC_FsC: RFC C-band short baseline flux density in Jy
# RFC_FlC: RFC C-band long baseline flux density in Jy
# RFC_FsX: RFC X-band short baseline flux density in Jy
# RFC_FlX: RFC X-band long baseline flux density in Jy
# RFC_FsU: RFC U-band short baseline flux density in Jy
# RFC_FlU: RFC U-band long baseline flux density in Jy
# RFC_FsK: RFC K-band short baseline flux density in Jy
# RFC_FlK: RFC K-band long baseline flux density in Jy
# VLASS_Flux_mJy: VLA Sky Survey flux density in mJy
# FIRST_InFootprint: Y/N indicator if source is in FIRST footprint
# FIRST_Flux_mJy: FIRST survey flux (mJy, -999 = non-detection in footprint)
# NVSS_InFootprint: Y/N indicator if source is in NVSS footprint
# NVSS_Flux_mJy: NVSS survey flux (mJy, -999 = non-detection in footprint)
#
"""
    
    # Write output with header
    with open(OUTPUT_FILE, 'w') as f:
        f.write(header_comments)
        
        # Write column numbers
        f.write("# Column Numbers:\n# ")
        f.write(",".join([str(i+1) for i in range(len(column_order))]))
        f.write("\n#\n")

    with open(P2OUTPUT_FILE, 'w') as f:
        f.write(header_comments)
        
        # Write column numbers
        f.write("# Column Numbers:\n# ")
        f.write(",".join([str(i+1) for i in range(len(column_order))]))
        f.write("\n#\n")

    with open(N2OUTPUT_FILE, 'w') as f:
        f.write(header_comments)
        
        # Write column numbers
        f.write("# Column Numbers:\n# ")
        f.write(",".join([str(i+1) for i in range(len(column_order))]))
        f.write("\n#\n")

    # NEW: write LDA>2 Nyland split files
    with open(KU_ONLY_LDA2_FILE, 'w') as f:
        f.write(header_comments)
        f.write("# Column Numbers:\n# ")
        f.write(",".join([str(i+1) for i in range(len(column_order))]))
        f.write("\n#\n")

    with open(C_KU_LDA2_FILE, 'w') as f:
        f.write(header_comments)
        f.write("# Column Numbers:\n# ")
        f.write(",".join([str(i+1) for i in range(len(column_order))]))
        f.write("\n#\n")
    
    # Append data
    output_df.to_csv(OUTPUT_FILE, mode='a', index=False)
    p2output_df.to_csv(P2OUTPUT_FILE, mode='a', index=False)
    n2output_df.to_csv(N2OUTPUT_FILE, mode='a', index=False)

    # NEW: write additional files
    ku_only_lda2_df.to_csv(KU_ONLY_LDA2_FILE, mode='a', index=False)
    c_ku_lda2_df.to_csv(C_KU_LDA2_FILE, mode='a', index=False)
    
    print(f"\n5. Output written to: {OUTPUT_FILE}")
    print(f"\n5. LDA>+2 Output written to: {P2OUTPUT_FILE}")
    print(f"\n5. LDA<-1 Output written to: {N2OUTPUT_FILE}")
    print(f"\n5. LDA>2 Nyland-present output written to: {KU_ONLY_LDA2_FILE}")
    print(f"\n5. LDA>2 Nyland-absent output written to: {C_KU_LDA2_FILE}")
    print(f"   Total sources processed: {len(output_df)}")
    print(f"   Total LDA>0: {n_lda_positive}")
    print(f"   Total LDA>2: {n_lda_plusone}")
    print(f"   Total LDA<-1: {n_lda_negone}")
    print(f"{nyland_massive_count} objects flagged for Ku-only due to presence in the NylandMASSIVE file.")
    print(f"{nyland_atlas_count} objects flagged for Ku-only due to presence in the NylandATLAS file.")
    print("\n" + "=" * 80)
    print("Processing complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
