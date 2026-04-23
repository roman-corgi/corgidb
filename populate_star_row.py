import pandas as pd
import numpy as np
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u
import warnings

# Suppress astroquery warnings for cleaner output
warnings.filterwarnings('ignore', category=UserWarning)


def populate_star_row(star_id):
    """
    Query astronomical catalogs and return a populated row for Stars table.

    """
    
    print(f"Querying {star_id}...")
    
    # Initialize data dictionary with all Stars table columns
    data = initialize_empty_row()
    
    # Query SIMBAD for basic astrometric and photometric data
    simbad_data = query_simbad_basic(star_id)
    if simbad_data is None:
        print(f"  ✗ {star_id} not found in SIMBAD")
        return None
    data.update(simbad_data)
    print(f"  ✓ SIMBAD basic data retrieved")
    
    # Query SIMBAD Collections for stellar parameters
    collections_data = query_simbad_collections(star_id)
    if collections_data is not None:
        data.update(collections_data)
        print(f"  ✓ SIMBAD Collections retrieved")
    
    # Query Vizier catalogs for additional photometry
    vizier_data = query_vizier_all(star_id)
    if vizier_data is not None:
        data.update(vizier_data)
        print(f"  ✓ Vizier catalogs queried")
    
    # Query NASA Exoplanet Archive for stellar parameters (planet hosts)
    exoplanet_data = query_exoplanet_archive(star_id)
    if exoplanet_data is not None:
        # Use Exoplanet Archive data (authoritative for planet hosts)
        data.update(exoplanet_data)
        print(f"  ✓ NASA Exoplanet Archive queried")
    
    # Calculate derived fields from retrieved data.
    derived_data = calculate_derived_fields(data)
    data.update(derived_data)
    print(f"  ✓ Derived fields calculated")
    
    # Convert to Series and return
    row = pd.Series(data)
    
    return row


def initialize_empty_row():
    """Initialize dictionary with all Stars table columns as NaN"""

  
    
    columns = [
        'st_id', 'simbad_name', 'st_name',
        'hd_name', 'hip_name', 'tic_id', 'gaia_id',
        'sy_snum', 'sy_pnum', 'sy_mnum', 'sy_dnum',
        'st_spectype', 'st_spectype_reflink',
        'st_teff', 'st_tefferr1', 'st_tefferr2', 'st_tefflim', 'st_teff_reflink',
        'st_rad', 'st_raderr1', 'st_raderr2', 'st_radlim', 'st_rad_reflink',
        'st_mass', 'st_masserr1', 'st_masserr2', 'st_masslim', 'st_mass_reflink',
        'st_met', 'st_meterr1', 'st_meterr2', 'st_metlim', 'st_met_reflink', 'st_metratio',
        'st_lum', 'st_lumerr1', 'st_lumerr2', 'st_lumlim', 'st_lum_reflink',
        'st_logg', 'st_loggerr1', 'st_loggerr2', 'st_logglim', 'st_logg_reflink',
        'st_age', 'st_ageerr1', 'st_ageerr2', 'st_agelim', 'st_age_reflink',
        'st_dens', 'st_denserr1', 'st_denserr2', 'st_denslim', 'st_dens_reflink',
        'st_vsin', 'st_vsinerr1', 'st_vsinerr2', 'st_vsinlim', 'st_vsin_reflink',
        'st_rotp', 'st_rotperr1', 'st_rotperr2', 'st_rotplim', 'st_rotp_reflink',
        'st_radv', 'st_radverr1', 'st_radverr2', 'st_radvlim', 'st_radv_reflink',
        'rastr', 'ra', 'decstr', 'dec',
        'raerr1', 'raerr2', 'decerr1', 'decerr2',
        'glon', 'glat', 'elon', 'elat',
        'glonerr1', 'glonerr2', 'glaterr1', 'glaterr2',
        'elonerr1', 'elonerr2', 'elaterr1', 'elaterr2',
        'ra_reflink',
        'sy_pm', 'sy_pmra', 'sy_pmdec',
        'sy_pmerr1', 'sy_pmerr2', 'sy_pmraerr1', 'sy_pmraerr2',
        'sy_pmdecerr1', 'sy_pmdecerr2', 'sy_pm_reflink',
        'sy_dist', 'sy_disterr1', 'sy_disterr2', 'sy_distlim', 'sy_dist_reflink',
        'sy_plx', 'sy_plxerr1', 'sy_plxerr2', 'sy_plxlim', 'sy_plx_reflink',
        'sy_bmag', 'sy_bmagerr1', 'sy_bmagerr2', 'sy_bmaglim', 'sy_bmag_reflink',
        'sy_vmag', 'sy_vmagerr1', 'sy_vmagerr2', 'sy_vmaglim', 'sy_vmag_reflink',
        'sy_jmag', 'sy_jmagerr1', 'sy_jmagerr2', 'sy_jmaglim', 'sy_jmag_reflink',
        'sy_hmag', 'sy_hmagerr1', 'sy_hmagerr2', 'sy_hmaglim', 'sy_hmag_reflink',
        'sy_kmag', 'sy_kmagerr1', 'sy_kmagerr2', 'sy_kmaglim', 'sy_kmag_reflink',
        'sy_umag', 'sy_umagerr1', 'sy_umagerr2', 'sy_umaglim', 'sy_umag_reflink',
        'sy_gmag', 'sy_gmagerr1', 'sy_gmagerr2', 'sy_gmaglim', 'sy_gmag_reflink',
        'sy_rmag', 'sy_rmagerr1', 'sy_rmagerr2', 'sy_rmaglim', 'sy_rmag_reflink',
        'sy_imag', 'sy_imagerr1', 'sy_imagerr2', 'sy_imaglim', 'sy_imag_reflink',
        'sy_zmag', 'sy_zmagerr1', 'sy_zmagerr2', 'sy_zmaglim', 'sy_zmag_reflink',
        'sy_w1mag', 'sy_w1magerr1', 'sy_w1magerr2', 'sy_w1maglim', 'sy_w1mag_reflink',
        'sy_w2mag', 'sy_w2magerr1', 'sy_w2magerr2', 'sy_w2maglim', 'sy_w2mag_reflink',
        'sy_w3mag', 'sy_w3magerr1', 'sy_w3magerr2', 'sy_w3maglim', 'sy_w3mag_reflink',
        'sy_w4mag', 'sy_w4magerr1', 'sy_w4magerr2', 'sy_w4maglim', 'sy_w4mag_reflink',
        'sy_gaiamag', 'sy_gaiamagerr1', 'sy_gaiamagerr2', 'sy_gaiamaglim', 'sy_gaiamag_reflink',
        'sy_icmag', 'sy_icmagerr1', 'sy_icmagerr2', 'sy_icmaglim', 'sy_icmag_reflink',
        'sy_tmag', 'sy_tmagerr1', 'sy_tmagerr2', 'sy_tmaglim', 'sy_tmag_reflink',
        'sy_kepmag', 'sy_kepmagerr1', 'sy_kepmagerr2', 'sy_kepmaglim', 'sy_kepmag_reflink',
        'reference_target', 'planet_host', 'disk_host',
        'calibration_target', 'engineering_target'
    ]
    
    return {col: np.nan for col in columns}


def query_simbad_basic(star_id):
    """
    Query SIMBAD for basic astrometric and photometric data

    """
    
    try:
        custom_simbad = Simbad()
        custom_simbad.add_votable_fields('ids', 'sp_type')
        custom_simbad.add_votable_fields('pmra', 'pmdec', 'plx_value', 'rvz_radvel')
        custom_simbad.add_votable_fields('B', 'V', 'J', 'H', 'K', 'I', 'G')
        
        result = custom_simbad.query_object(star_id)
        
        if result is None:
            return None
        
        data = {}
        
        # Main identifier and coordinates
        data['simbad_name'] = result['main_id'][0]
        data['ra'] = result['ra'][0]
        data['dec'] = result['dec'][0]
        
        # Format RA/Dec as sexagesimal strings
        coord = SkyCoord(ra=data['ra']*u.deg, dec=data['dec']*u.deg, frame='icrs')
        data['rastr'] = coord.ra.to_string(unit=u.hourangle, sep=' ', precision=2)
        data['decstr'] = coord.dec.to_string(unit=u.deg, sep=' ', precision=2)
        
        # Proper motion (check for masked values)
        if 'pmra' in result.colnames and not isinstance(result['pmra'][0], np.ma.core.MaskedConstant):
            data['sy_pmra'] = result['pmra'][0]
        if 'pmdec' in result.colnames and not isinstance(result['pmdec'][0], np.ma.core.MaskedConstant):
            data['sy_pmdec'] = result['pmdec'][0]
        
        # Parallax
        if 'plx_value' in result.colnames and not isinstance(result['plx_value'][0], np.ma.core.MaskedConstant):
            data['sy_plx'] = result['plx_value'][0]
        
        # Radial velocity
        if 'rvz_radvel' in result.colnames and not isinstance(result['rvz_radvel'][0], np.ma.core.MaskedConstant):
            data['st_radv'] = result['rvz_radvel'][0]
        
        # Spectral type
        if 'sp_type' in result.colnames and not isinstance(result['sp_type'][0], np.ma.core.MaskedConstant):
            data['st_spectype'] = result['sp_type'][0]
        
        # Photometry in multiple bands
        flux_map = {'B': 'sy_bmag', 'V': 'sy_vmag', 'J': 'sy_jmag', 
                    'H': 'sy_hmag', 'K': 'sy_kmag', 'I': 'sy_icmag', 'G': 'sy_gaiamag'}
        
        for band, field in flux_map.items():
            if band in result.colnames and not isinstance(result[band][0], np.ma.core.MaskedConstant):
                data[field] = result[band][0]
        
        # Parse catalog cross-identifications
        if 'ids' in result.colnames:
            catalog_ids = parse_catalog_ids(result['ids'][0])
            data.update(catalog_ids)
        
        return data
        
    except Exception as e:
        print(f"    Error in SIMBAD basic: {e}")
        return None


def parse_catalog_ids(ids_string):
    """Extract HD, HIP, TIC, Gaia IDs from SIMBAD IDS field"""
    
    ids = {}
    
    for id_str in ids_string.split('|'):
        id_str = id_str.strip()
        if id_str.startswith('HD '):
            ids['hd_name'] = id_str
        elif id_str.startswith('HIP '):
            ids['hip_name'] = id_str
        elif id_str.startswith('TIC '):
            ids['tic_id'] = id_str
        elif id_str.startswith('Gaia DR3 '):
            ids['gaia_id'] = id_str
    
    return ids


def query_simbad_collections(star_id):
    """
    Query SIMBAD Collections for stellar parameters

    """
    
    try:
        collections_simbad = Simbad()
        collections_simbad.add_votable_fields('mesfe_h', 'mesrot')
        
        result = collections_simbad.query_object(star_id)
        
        if result is None:
            return None
        
        data = {}
        
        # Convert result to pandas DataFrame.
        measurements_df = result.to_pandas()
        
        # Stellar parameters from Fe_H measurements table. Newest measurement with complete parameter set.
        # Filter for rows with complete parameter sets (Teff + log g + [Fe/H]).
        complete_params = measurements_df[
            (measurements_df['mesfe_h.teff'].notna()) &
            (~measurements_df['mesfe_h.teff'].apply(lambda x: isinstance(x, np.ma.core.MaskedConstant))) &
            (measurements_df['mesfe_h.log_g'].notna()) &
            (~measurements_df['mesfe_h.log_g'].apply(lambda x: isinstance(x, np.ma.core.MaskedConstant))) &
            (measurements_df['mesfe_h.fe_h'].notna()) &
            (~measurements_df['mesfe_h.fe_h'].apply(lambda x: isinstance(x, np.ma.core.MaskedConstant)))
        ]
        
        if len(complete_params) > 0:
            # Sort by bibcode (newest first - bibcodes are YYYYJJJJJ format). Take the most recent measurement with complete parameter set.
            best_measurement = complete_params.sort_values(
                'mesfe_h.bibcode', 
                ascending=False
            ).iloc[0]
            
            # Extract parameters from single measurement.
            data['st_teff'] = best_measurement['mesfe_h.teff']
            data['st_logg'] = best_measurement['mesfe_h.log_g']
            data['st_met'] = best_measurement['mesfe_h.fe_h']
            
            if 'mesfe_h.compstar' in measurements_df.columns:
                if not isinstance(best_measurement['mesfe_h.compstar'], np.ma.core.MaskedConstant):
                    data['st_metratio'] = best_measurement['mesfe_h.compstar']
        
        else:
            # No complete sets available; fall back to taking newest available values.
            if 'mesfe_h.teff' in measurements_df.columns:
                teff_data = measurements_df.dropna(subset=['mesfe_h.teff'])
                if len(teff_data) > 0:
                    newest_teff = teff_data.sort_values('mesfe_h.bibcode', ascending=False).iloc[0]
                    if not isinstance(newest_teff['mesfe_h.teff'], np.ma.core.MaskedConstant):
                        data['st_teff'] = newest_teff['mesfe_h.teff']
            
            if 'mesfe_h.log_g' in measurements_df.columns:
                logg_data = measurements_df.dropna(subset=['mesfe_h.log_g'])
                if len(logg_data) > 0:
                    newest_logg = logg_data.sort_values('mesfe_h.bibcode', ascending=False).iloc[0]
                    if not isinstance(newest_logg['mesfe_h.log_g'], np.ma.core.MaskedConstant):
                        data['st_logg'] = newest_logg['mesfe_h.log_g']
            
            if 'mesfe_h.fe_h' in measurements_df.columns:
                feh_data = measurements_df.dropna(subset=['mesfe_h.fe_h'])
                if len(feh_data) > 0:
                    newest_feh = feh_data.sort_values('mesfe_h.bibcode', ascending=False).iloc[0]
                    if not isinstance(newest_feh['mesfe_h.fe_h'], np.ma.core.MaskedConstant):
                        data['st_met'] = newest_feh['mesfe_h.fe_h']
        
        # Rotational velocity from rotation measurements table. Rotation period not available in mesrot table.
        if 'mesrot.vsini' in measurements_df.columns:
            vsini_data = measurements_df.dropna(subset=['mesrot.vsini'])
            if len(vsini_data) > 0:
                newest_vsini = vsini_data.sort_values('mesrot.bibcode', ascending=False).iloc[0]
                if not isinstance(newest_vsini['mesrot.vsini'], np.ma.core.MaskedConstant):
                    data['st_vsin'] = newest_vsini['mesrot.vsini']
                    
                    # Include error if available
                    if 'mesrot.vsini_err' in measurements_df.columns:
                        err = newest_vsini['mesrot.vsini_err']
                        if not isinstance(err, np.ma.core.MaskedConstant):
                            data['st_vsinerr1'] = err
                            data['st_vsinerr2'] = -err  # Symmetric error
        
        return data
        
    except Exception as e:
        print(f"    Error in SIMBAD Collections: {e}")
        return None


def query_vizier_all(star_id):
    """
    Query Vizier catalogs for additional photometry and stellar parameters

    """
    
    data = {}
    
    # AllWISE (II/328) - WISE infrared photometry
    try:
        v = Vizier(columns=["**"], row_limit=1)
        result = v.query_object(star_id, catalog="II/328")
        
        if result and len(result) > 0:
            for i, band in enumerate(['W1mag', 'W2mag', 'W3mag', 'W4mag'], 1):
                if band in result[0].colnames:
                    data[f'sy_w{i}mag'] = result[0][band][0]
                    err_col = f'e_{band}'
                    if err_col in result[0].colnames:
                        err = result[0][err_col][0]
                        data[f'sy_w{i}magerr1'] = err
                        data[f'sy_w{i}magerr2'] = -err
    except:
        pass
    
    # SDSS DR16 (V/154) - Sloan optical photometry
    try:
        v = Vizier(columns=["**"], row_limit=1)
        result = v.query_object(star_id, catalog="V/154")
        
        if result and len(result) > 0:
            sloan_map = {'umag': 'sy_umag', 'gmag': 'sy_gmag', 'rmag': 'sy_rmag',
                        'imag': 'sy_imag', 'zmag': 'sy_zmag'}
            for sdss_col, our_col in sloan_map.items():
                if sdss_col in result[0].colnames:
                    data[our_col] = result[0][sdss_col][0]
                    err_col = f'e_{sdss_col}'
                    if err_col in result[0].colnames:
                        err = result[0][err_col][0]
                        data[f'{our_col}err1'] = err
                        data[f'{our_col}err2'] = -err
    except:
        pass
    
    # TIC v8 (IV/38) - TESS magnitudes and stellar parameters
    try:
        v = Vizier(columns=["**"], row_limit=1)
        result = v.query_object(star_id, catalog="IV/38/tic")
        
        if result and len(result) > 0:
            # TESS magnitude
            if 'Tmag' in result[0].colnames:
                data['sy_tmag'] = result[0]['Tmag'][0]
            
            # Stellar mass
            if 'Mass' in result[0].colnames:
                data['st_mass'] = result[0]['Mass'][0]
                if 'E_Mass' in result[0].colnames:
                    err = result[0]['E_Mass'][0]
                    data['st_masserr1'] = err
                    data['st_masserr2'] = -err
            
            # Stellar radius
            if 'Rad' in result[0].colnames:
                data['st_rad'] = result[0]['Rad'][0]
                if 'E_Rad' in result[0].colnames:
                    err = result[0]['E_Rad'][0]
                    data['st_raderr1'] = err
                    data['st_raderr2'] = -err
    except:
        pass
    
    return data if data else None


def query_exoplanet_archive(star_id):
    """
    Query NASA Exoplanet Archive for stellar mass and radius.    
    Attempts to retrieve stellar parameters for planet-hosting stars.
    """
    
    try:
        from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
        
        # Query the Planetary Systems Composite Parameters table.
        result = NasaExoplanetArchive.query_object(star_id, table='pscomppars')
        
        if result is not None and len(result) > 0:
            data = {}
            
            # Stellar mass
            if 'st_mass' in result.colnames:
                mass_val = result['st_mass'][0]
                # Extract value from Quantity object if needed
                if hasattr(mass_val, 'value'):
                    mass_val = mass_val.value
                
                if not np.isnan(mass_val):
                    data['st_mass'] = float(mass_val)
                    
                    # Include errors if available
                    if 'st_masserr1' in result.colnames:
                        err1 = result['st_masserr1'][0]
                        if hasattr(err1, 'value'):
                            err1 = err1.value
                        if not np.isnan(err1):
                            data['st_masserr1'] = float(err1)
                    
                    if 'st_masserr2' in result.colnames:
                        err2 = result['st_masserr2'][0]
                        if hasattr(err2, 'value'):
                            err2 = err2.value
                        if not np.isnan(err2):
                            data['st_masserr2'] = float(err2)
            
            # Stellar radius
            if 'st_rad' in result.colnames:
                rad_val = result['st_rad'][0]
                # Extract value from Quantity object if needed
                if hasattr(rad_val, 'value'):
                    rad_val = rad_val.value
                
                if not np.isnan(rad_val):
                    data['st_rad'] = float(rad_val)
                    
                    # Include errors if available
                    if 'st_raderr1' in result.colnames:
                        err1 = result['st_raderr1'][0]
                        if hasattr(err1, 'value'):
                            err1 = err1.value
                        if not np.isnan(err1):
                            data['st_raderr1'] = float(err1)
                    
                    if 'st_raderr2' in result.colnames:
                        err2 = result['st_raderr2'][0]
                        if hasattr(err2, 'value'):
                            err2 = err2.value
                        if not np.isnan(err2):
                            data['st_raderr2'] = float(err2)
            
            return data if data else None
        
        return None
        
    except Exception as e:
        # Star may not be in archive or module not available.
        return None
    
    
def calculate_derived_fields(data):
    """
    Calculate derived astronomical quantities

    """
    
    derived = {}
    
    # Transform coordinates to Galactic and ecliptic systems
    if not np.isnan(data.get('ra', np.nan)) and not np.isnan(data.get('dec', np.nan)):
        coord = SkyCoord(ra=data['ra']*u.deg, dec=data['dec']*u.deg, frame='icrs')
        
        # Galactic coordinates
        galactic = coord.galactic
        derived['glon'] = galactic.l.deg
        derived['glat'] = galactic.b.deg
        
        # Ecliptic coordinates
        ecliptic = coord.transform_to('geocentrictrueecliptic')
        derived['elon'] = ecliptic.lon.deg
        derived['elat'] = ecliptic.lat.deg
    
    # Total proper motion from components
    pmra = data.get('sy_pmra', np.nan)
    pmdec = data.get('sy_pmdec', np.nan)
    if not (np.isnan(pmra) or np.isnan(pmdec)):
        derived['sy_pm'] = np.sqrt(pmra**2 + pmdec**2)
    
    # Distance from parallax (convert mas to pc)
    plx = data.get('sy_plx', np.nan)
    if not np.isnan(plx) and plx > 0:
        derived['sy_dist'] = 1000.0 / plx
    
    # Stellar density from mass and radius
    mass = data.get('st_mass', np.nan)
    radius = data.get('st_rad', np.nan)
    if not (np.isnan(mass) or np.isnan(radius)) and radius > 0:
        volume = (4.0/3.0) * np.pi * (radius**3)
        derived['st_dens'] = mass / volume
    
    return derived


# Test/demonstration code
if __name__ == "__main__":
    print("Testing populate_star_row function...")
    print("=" * 60)
    
    row = populate_star_row("HD 209458")
    
    if row is not None:
        print(f"\n{'='*60}")
        print(f"SUCCESS! Retrieved data for {row['simbad_name']}")
        print(f"{'='*60}")
        print(f"\nCoordinates:")
        print(f"  RA:  {row['ra']:.6f}° ({row['rastr']})")
        print(f"  Dec: {row['dec']:.6f}° ({row['decstr']})")
        print(f"  Galactic: l={row['glon']:.2f}°, b={row['glat']:.2f}°")
        print(f"  Ecliptic: lon={row['elon']:.2f}°, lat={row['elat']:.2f}°")
        
        print(f"\nStellar Parameters:")
        if not np.isnan(row['st_teff']):
            print(f"  Teff: {row['st_teff']:.0f} K")
        if not np.isnan(row['st_logg']):
            print(f"  log g: {row['st_logg']:.2f}")
        if not np.isnan(row['st_met']):
            print(f"  [Fe/H]: {row['st_met']:.2f}")
        if not np.isnan(row['st_mass']) and not isinstance(row['st_mass'], np.ma.core.MaskedConstant):
            print(f"  Mass: {row['st_mass']:.3f} M☉")
        if not np.isnan(row['st_rad']) and not isinstance(row['st_rad'], np.ma.core.MaskedConstant):
            print(f"  Radius: {row['st_rad']:.3f} R☉")
        if not np.isnan(row['st_vsin']):
            print(f"  v*sin(i): {row['st_vsin']:.2f} km/s")
        
        print(f"\nPhotometry:")
        bands = [('sy_vmag', 'V'), ('sy_bmag', 'B'), ('sy_jmag', 'J'), 
                 ('sy_hmag', 'H'), ('sy_kmag', 'K'), ('sy_gaiamag', 'G'),
                 ('sy_w1mag', 'W1'), ('sy_tmag', 'T')]
        for field, name in bands:
            if not np.isnan(row[field]):
                print(f"  {name}: {row[field]:.3f}")
        
        populated = row.notna().sum()
        total = len(row)
        print(f"\n{'='*60}")
        print(f"Fields populated: {populated}/{total} ({100*populated/total:.1f}%)")
        print(f"{'='*60}")
    else:
        print("\n✗ Error: Failed to retrieve data")
