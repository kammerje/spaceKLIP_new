from __future__ import division

import matplotlib
matplotlib.rcParams.update({'font.size': 14})


# =============================================================================
# IMPORTS
# =============================================================================

import os
import pdb
import sys
import functools

import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)


# =============================================================================
# MAIN
# =============================================================================

def read_obs(fitsfile):
    
    # Read FITS file.
    hdul = pyfits.open(fitsfile)
    data = hdul['SCI'].data
    erro = hdul['ERR'].data
    pxdq = hdul['DQ'].data
    head_pri = hdul[0].header
    head_sci = hdul['SCI'].header
    hdul.close()
    is2d = False
    if data.ndim == 2:
        data = data[np.newaxis, :]
        erro = erro[np.newaxis, :]
        pxdq = pxdq[np.newaxis, :]
        is2d = True
    if data.ndim != 3:
        raise UserWarning('Requires 2D/3D data cube')
    
    return data, erro, pxdq, head_pri, head_sci, is2d

def write_obs(fitsfile,
              output_dir,
              data,
              erro,
              pxdq,
              head_pri,
              head_sci,
              is2d):
    
    # Write FITS file.
    hdul = pyfits.open(fitsfile)
    if is2d:
        hdul['SCI'].data = data[0]
        hdul['ERR'].data = erro[0]
        hdul['DQ'].data = pxdq[0]
    else:
        hdul['SCI'].data = data
        hdul['ERR'].data = erro
        hdul['DQ'].data = pxdq
    hdul[0].header = head_pri
    hdul['SCI'].header = head_sci
    fitsfile = os.path.join(output_dir, os.path.split(fitsfile)[1])
    hdul.writeto(fitsfile, output_verify='fix', overwrite=True)
    hdul.close()
    
    return fitsfile

def read_red(fitsfile):
    
    # Read FITS file.
    hdul = pyfits.open(fitsfile)
    data = hdul[0].data
    if data is None:
        try:
            data = hdul['SCI'].data
        except:
            raise UserWarning('Could not find any data')
    head_pri = hdul[0].header
    try:
        head_sci = hdul['SCI'].header
    except:
        head_sci = None
    hdul.close()
    is2d = False
    if data.ndim == 2:
        data = data[np.newaxis, :]
        is2d = True
    if data.ndim != 3:
        raise UserWarning('Requires 2D/3D data cube')
    
    return data, head_pri, head_sci, is2d


@functools.lru_cache()
def lookup_pixscale(instrument, detector):
    """Retrieve the pixel scale of a given high contrast mode from the SIAF.
    Used in importing files into the database.

    Approximations: We ignore the slight departures from ideal square pixels and just
    return the X scale. We also ignore the slight wavelength dependence of the pixel scale.

    Alternatively we could read pixel scales in from FITS WCS headers, with much the same result.

    """

    import pysiaf
    siaf = pysiaf.Siaf(instrument)

    if detector.upper() == 'NRCALONG':
        return siaf.apertures['NRCA5_MASK430R'].XSciScale
    elif detector.upper() == 'NRCA2':
        return siaf.apertures['NRCA2_MASK210R'].XSciScale
    elif detector.upper() == 'NRCA4':
        return siaf.apertures['NRCA4_MASKSWB'].XSciScale
    elif detector.upper() == 'NIRISS':
        return siaf.apertures['NIS_AMI1'].XSciScale
    else:
        return siaf.apertures['MIRIM_MASK1065'].XSciScale
