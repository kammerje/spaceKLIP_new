from __future__ import division

import matplotlib
matplotlib.rcParams.update({'font.size': 14})


# =============================================================================
# IMPORTS
# =============================================================================

import os
import pdb
import sys

import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np

import importlib
import webbpsf_ext

import astropy.units as u

from astroquery.svo_fps import SvoFps
from synphot import Observation, SourceSpectrum, SpectralElement
from synphot.models import Empirical1D
from synphot.units import convert_flux

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)


# =============================================================================
# MAIN
# =============================================================================

def read_spec_file(starfile):
    """
    Read a spectrum from a TXT file.
    
    Parameters
    ----------
    starfile : path
        Path of two column TXT file with wavelength (micron) and flux (Jy).
    
    Returns
    -------
    sed : synphot.SourceSpectrum
        Spectrum of the source.
    
    """
    
    try:
        data = np.genfromtxt(starfile).transpose()
        model_wave = data[0]
        model_flux = data[1]    
        sed = SourceSpectrum(Empirical1D, points=model_wave << u.Unit('micron'), lookup_table=model_flux << u.Unit('Jy'))
        sed.meta['name'] = starfile.split('/')[-1]
    except:
        raise ValueError('Unable to read provided starfile; ensure format is in two columns with wavelength (microns), flux (Jy)')
    
    return sed

def get_stellar_magnitudes(starfile,
                           spectral_type,
                           instrume,
                           return_si=False,
                           output_dir=None):
    """
    Get the source brightness and zero point fluxes in each filter of the JWST
    instrument in use.
    
    Parameters
    ----------
    starfile : path
        Path of VizieR VOTable containing host star photometry or two
        column TXT file with wavelength (micron) and flux (Jy).
    spectral_type : str, optional
        Host star spectral type for the stellar model SED. The default is
        'G2V'.
    instrume : 'NIRCAM', 'NIRISS', or 'MIRI'
        JWST instrument in use.
    return_si : bool, optional
        Return the filter zero point in SI units in addition to Jy? The default
        is False.
    output_dir : path, optional
        Path of the directory where the SED plot shall be saved. The default is
        None.
    
    Returns
    -------
    mstar : dict
        Dictionary of the source brightness (vegamag) in each filter of the
        JWST instrument in use.
    fzero : dict
        Dictionary of the zero point flux (Jy) of each filter of the JWST
        instrument in use.
    fzero_si : dict, optional
        Dictionary of the zero point flux (erg/cm^2/s/A) of each filter of the
        JWST instrument in use.
    
    """
    
    # VOTable.
    if starfile[-4:] == '.vot':
        
        # Initialize SED in random bandpass and at random magnitude.
        bp_k = webbpsf_ext.bp_2mass('k')
        bp_mag = 5.
        try:
            spec = webbpsf_ext.spectra.source_spectrum(name='Input Data & SED', sptype=spectral_type, mag_val=bp_mag, bp=bp_k, votable_file=starfile)
        except:
            spec = webbpsf_ext.spectra.source_spectrum(name='Input Data & SED', sptype=spectral_type, mag_val=bp_mag, bp=bp_k, votable_input=starfile)
        
        # Split between NIR and MIR exposures.
        if instrume == 'MIRI':
            wlim = [5., 30.]
        else:
            wlim = [1., 5.]
        
        # Fit SED to input photometry and plot SED.
        spec.fit_SED(x0=[1.], wlim=wlim, use_err=False, verbose=False)
        if output_dir is not None:
            spec.plot_SED()
            plt.savefig(os.path.join(output_dir, 'sed.pdf'))
            plt.close()
        
        # Convert units to photlam.
        input_flux = u.Quantity(spec.sp_model.flux, str(spec.sp_model.fluxunits))
        photlam_flux = convert_flux(spec.sp_model.wave, input_flux, out_flux_unit='photlam')
        sed = SourceSpectrum(Empirical1D, points=spec.sp_model.wave << u.Unit(str(spec.sp_model.waveunits)), lookup_table=photlam_flux << u.Unit('photlam'))
    
    # TXT file.
    else:
        sed = read_spec_file(starfile)
    
    # Load respective filters from the SVO Filter Profile Service.
    # http://svo2.cab.inta-csic.es/theory/fps/
    filts = []
    zeros = []
    filter_list = SvoFps.get_filter_list(facility='JWST', instrument=instrume)
    for i in range(len(filter_list)):
        filts += [filter_list['filterID'][i].split('.')[-1]]
        zeros += [filter_list['ZeroPoint'][i]]
    zero_points_si = {'F182M': 7.44007e-11,
                      'F210M': 4.69758e-11,
                      'F250M': 2.41440e-11,
                      'F300M': 1.24029e-11,
                      'F335M': 7.92772e-12,
                      'F356W': 6.38971e-12,
                      'F444W': 2.84527e-12}
    
    # Compute magnitude in each filter.
    mstar = {}  # vegamag
    fzero = {}  # Jy
    fzero_si = {}  # erg/cm^2/s/A
    for i, filt in enumerate(filts):
        
        # Read bandpass.
        try:
            with importlib.resources.open_text(f'spaceKLIP.resources.PCEs.{instrume}', f'{filt}.txt') as bandpass_file:
                bandpass_data = np.genfromtxt(bandpass_file).transpose()
                bandpass_wave = bandpass_data[0] * 1e4  # Angstrom
                bandpass_throughput = bandpass_data[1]
        except FileNotFoundError:
            continue
        
        # Create bandpass object.
        bandpass = SpectralElement(Empirical1D, points=bandpass_wave, lookup_table=bandpass_throughput)
        
        # Compute magnitude.
        obs = Observation(sed, bandpass, binset=bandpass.waveset)
        vegased = SourceSpectrum.from_vega()
        mag = obs.effstim(flux_unit='vegamag', vegaspec=vegased).value
        mstar[filt.upper()] = mag
        fzero[filt.upper()] = zeros[i]
        try:
            fzero_si[filt.upper()] = zero_points_si[filt.upper()]
        except KeyError:
            fzero_si[filt.upper()] = np.nan
    
    if return_si:
        return mstar, fzero, fzero_si
    else:
        return mstar, fzero
