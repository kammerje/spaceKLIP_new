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

from copy import deepcopy
from jwst.pipeline import Detector1Pipeline, Image2Pipeline, Coron3Pipeline
from scipy.ndimage import fourier_shift, median_filter
from scipy.ndimage import shift as spline_shift
from scipy.optimize import leastsq, minimize
from spaceKLIP import utils as ut

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

from webbpsf_ext import robust


# =============================================================================
# MAIN
# =============================================================================

class ImageTools():
    """
    The spaceKLIP image manipulation tools class.
    """
    
    def __init__(self,
                 database):
        
        # Make an internal alias of the spaceKLIP database class.
        self.database = database
        
        pass
    
    def remove_frames(self,
                      index=[0],
                      frame=None,
                      types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                      subdir='removed'):
        """ Remove the first frame due to reset switch charge delay. Only required
               for MIRI.

        Parameters
        ----------
        index
        frame
        types
        subdir

        Returns None, but writes out new files to subdir and updates database.
        """

        # Check input.
        if isinstance(index, int):
            index = [index]
        
        # Set output directory.
        output_dir = os.path.join(self.database.output_dir, subdir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Loop through concatenations.
        for i, key in enumerate(self.database.obs.keys()):
            log.info('--> Concatenation ' + key)
            
            # Loop through FITS files.
            nfitsfiles = len(self.database.obs[key])
            for j in range(nfitsfiles):
                
                # Read FITS file.
                fitsfile = self.database.obs[key]['FITSFILE'][j]
                data, erro, pxdq, head_pri, head_sci, is2d = ut.read_obs(fitsfile)
                nints = self.database.obs[key]['NINTS'][j]
                
                # Skip file types that are not in the list of types.
                if self.database.obs[key]['TYPE'][j] in types:
                    
                    # Remove frames.
                    try:
                        index_temp = frame[key][j]
                    except:
                        index_temp = index.copy()
                    head, tail = os.path.split(fitsfile)
                    log.info('  --> Frame removal: ' + tail)
                    log.info('  --> Frame removal: removing frame(s) ' + str(index_temp))
                    data = np.delete(data, index_temp, axis=0)
                    erro = np.delete(erro, index_temp, axis=0)
                    pxdq = np.delete(pxdq, index_temp, axis=0)
                    nints = data.shape[0]
                
                # Write FITS file.
                head_pri['NINTS'] = nints
                fitsfile = ut.write_obs(fitsfile, output_dir, data, erro, pxdq, head_pri, head_sci, is2d)
                
                # Update spaceKLIP database.
                self.database.update_obs(key, j, fitsfile, nints=nints)
        
        pass
    
    def crop_frames(self,
                    npix=1,
                    types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                    subdir='cropped'):
        """ Crop frames

        Parameters
        ----------
        npix
        types
        subdir

        Returns None, but writes out new files to subdir and updates database.

        """
        
        # Check input.
        if isinstance(npix, int):
            npix = [npix, npix, npix, npix]  # left, right, bottom, top
        if len(npix) != 4:
            raise UserWarning('Parameter npix must either be an int or a list of four int (left, right, bottom, top)')
        
        # Set output directory.
        output_dir = os.path.join(self.database.output_dir, subdir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Loop through concatenations.
        for i, key in enumerate(self.database.obs.keys()):
            log.info('--> Concatenation ' + key)
            
            # Loop through FITS files.
            nfitsfiles = len(self.database.obs[key])
            for j in range(nfitsfiles):
                
                # Read FITS file.
                fitsfile = self.database.obs[key]['FITSFILE'][j]
                data, erro, pxdq, head_pri, head_sci, is2d = ut.read_obs(fitsfile)
                crpix1 = self.database.obs[key]['CRPIX1'][j]
                crpix2 = self.database.obs[key]['CRPIX2'][j]
                
                # Skip file types that are not in the list of types.
                if self.database.obs[key]['TYPE'][j] in types:
                    
                    # Crop frames.
                    head, tail = os.path.split(fitsfile)
                    log.info('  --> Frame cropping: ' + tail)
                    sh = data.shape
                    data = data[:, npix[2]:-npix[3], npix[0]:-npix[1]]
                    erro = erro[:, npix[2]:-npix[3], npix[0]:-npix[1]]
                    pxdq = pxdq[:, npix[2]:-npix[3], npix[0]:-npix[1]]
                    crpix1 -= npix[0]
                    crpix2 -= npix[2]
                    log.info('  --> Frame cropping: old shape = ' + str(sh[1:]) + ', new shape = ' + str(data.shape[1:]))
                
                # Write FITS file.
                head_sci['CRPIX1'] = crpix1
                head_sci['CRPIX2'] = crpix2
                fitsfile = ut.write_obs(fitsfile, output_dir, data, erro, pxdq, head_pri, head_sci, is2d)
                
                # Update spaceKLIP database.
                self.database.update_obs(key, j, fitsfile, crpix1=crpix1, crpix2=crpix2)
        
        pass
    
    def pad_frames(self,
                   npix=1,
                   cval=np.nan,
                   types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                   subdir='padded'):
        """ Pad frames with additional pixels around the border

        Parameters
        ----------
        npix : int, or iterable of 4 ints
            number of additional pixels to add per side. If scalar, same number of rows/columns will be added to
            all 4 sides. Alternatively provide an iterable to set [left, right, top, bottom] independently
        cval : float
            Constant value to write to the padded pixels
        types
        subdir

        Returns None, but writes out new files to subdir and updates database.

        """
        
        # Check input.
        if isinstance(npix, int):
            npix = [npix, npix, npix, npix]  # left, right, bottom, top
        if len(npix) != 4:
            raise UserWarning('Parameter npix must either be an int or a list of four int (left, right, bottom, top)')
        
        # Set output directory.
        output_dir = os.path.join(self.database.output_dir, subdir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Loop through concatenations.
        for i, key in enumerate(self.database.obs.keys()):
            log.info('--> Concatenation ' + key)
            
            # Loop through FITS files.
            nfitsfiles = len(self.database.obs[key])
            for j in range(nfitsfiles):
                
                # Read FITS file.
                fitsfile = self.database.obs[key]['FITSFILE'][j]
                data, erro, pxdq, head_pri, head_sci, is2d = ut.read_obs(fitsfile)
                crpix1 = self.database.obs[key]['CRPIX1'][j]
                crpix2 = self.database.obs[key]['CRPIX2'][j]
                
                # Skip file types that are not in the list of types.
                if self.database.obs[key]['TYPE'][j] in types:
                    
                    # Crop frames.
                    head, tail = os.path.split(fitsfile)
                    log.info('  --> Frame padding: ' + tail)
                    sh = data.shape
                    data = np.pad(data, ((0, 0), (npix[2], npix[3]), (npix[0], npix[1])), mode='constant', constant_values=cval)
                    erro = np.pad(erro, ((0, 0), (npix[2], npix[3]), (npix[0], npix[1])), mode='constant', constant_values=cval)
                    pxdq = np.pad(pxdq, ((0, 0), (npix[2], npix[3]), (npix[0], npix[1])), mode='constant', constant_values=0)
                    crpix1 += npix[0]
                    crpix2 += npix[2]
                    log.info('  --> Frame padding: old shape = ' + str(sh[1:]) + ', new shape = ' + str(data.shape[1:]) + ', fill value = %.2f' % cval)
                
                # Write FITS file.
                head_sci['CRPIX1'] = crpix1
                head_sci['CRPIX2'] = crpix2
                fitsfile = ut.write_obs(fitsfile, output_dir, data, erro, pxdq, head_pri, head_sci, is2d)
                
                # Update spaceKLIP database.
                self.database.update_obs(key, j, fitsfile, crpix1=crpix1, crpix2=crpix2)
        
        pass
    
    def coadd_frames(self,
                     nframes=None,
                     types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                     subdir='coadded'):
        """ Coadd frames

        Parameters
        ----------
        nframes
        types
        subdir

        Returns None, but writes out new files to subdir and updates database.

        """
        
        # Set output directory.
        output_dir = os.path.join(self.database.output_dir, subdir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        #The starting value
        nframes0 = nframes
        
        # Loop through concatenations.
        for i, key in enumerate(self.database.obs.keys()):
            log.info('--> Concatenation ' + key)
            
            # Loop through FITS files.
            nfitsfiles = len(self.database.obs[key])
            for j in range(nfitsfiles):
                
                # Read FITS file.
                fitsfile = self.database.obs[key]['FITSFILE'][j]
                data, erro, pxdq, head_pri, head_sci, is2d = ut.read_obs(fitsfile)
                nints = self.database.obs[key]['NINTS'][j]
                effinttm = self.database.obs[key]['EFFINTTM'][j]
                
                #If nframes is not provided, collapse everything
                if nframes0 is None:
                    nframes = nints
                
                # Skip file types that are not in the list of types.
                if self.database.obs[key]['TYPE'][j] in types:
                    
                    # Coadd frames.
                    head, tail = os.path.split(fitsfile)
                    log.info('  --> Frame coadding: ' + tail)
                    ncoadds = data.shape[0] // nframes
                    data = np.nanmedian(data[:nframes * ncoadds].reshape((nframes, ncoadds, data.shape[-2], data.shape[-1])), axis=0)
                    erro_reshape = erro[:nframes * ncoadds].reshape((nframes, ncoadds, erro.shape[-2], erro.shape[-1]))
                    nsample = np.sum(np.logical_not(np.isnan(erro_reshape)), axis=0)
                    erro = np.true_divide(np.sqrt(np.nansum(erro_reshape**2, axis=0)), nsample)
                    pxdq_temp = pxdq[:nframes * ncoadds].reshape((nframes, ncoadds, pxdq.shape[-2], pxdq.shape[-1]))
                    
                    pxdq = pxdq_temp[0]
                    for k in range(1, nframes):
                        pxdq = np.bitwise_or(pxdq, pxdq_temp[k])
                    nints = data.shape[0]
                    effinttm *= nframes
                    log.info('  --> Frame coadding: %.0f coadd(s) of %.0f frame(s)' % (ncoadds, nframes))
                
                # Write FITS file.
                head_pri['NINTS'] = nints
                head_pri['EFFINTTM'] = effinttm
                fitsfile = ut.write_obs(fitsfile, output_dir, data, erro, pxdq, head_pri, head_sci, is2d)
                
                # Update spaceKLIP database.
                self.database.update_obs(key, j, fitsfile, nints=nints, effinttm=effinttm)
        
        pass
    
    def subtract_median(self,
                        types=['SCI', 'SCI_TA', 'SCI_BG', 'REF', 'REF_TA', 'REF_BG'],
                        subdir='medsub'):
        """ Median-subtract each frame to mitigate uncalibrated bias drifts.
        Generally only required for NIRCam.

        Returns None, but writes out new files to subdir and updates database.
        """
        
        # Set output directory.
        output_dir = os.path.join(self.database.output_dir, subdir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Loop through concatenations.
        for i, key in enumerate(self.database.obs.keys()):
            log.info('--> Concatenation ' + key)
            
            # Loop through FITS files.
            nfitsfiles = len(self.database.obs[key])
            for j in range(nfitsfiles):
                
                # Read FITS file.
                fitsfile = self.database.obs[key]['FITSFILE'][j]
                data, erro, pxdq, head_pri, head_sci, is2d = ut.read_obs(fitsfile)
                
                # Skip file types that are not in the list of types.
                if self.database.obs[key]['TYPE'][j] in types:
                    
                    # Subtract median.
                    head, tail = os.path.split(fitsfile)
                    log.info('  --> Median subtraction: ' + tail)
                    data_temp = data.copy()
                    # if self.database.obs[key]['TELESCOP'][j] == 'JWST' and self.database.obs[key]['INSTRUME'][j] == 'NIRCAM':
                    data_temp[pxdq != 0] = np.nan
                    # else:
                    #     data_temp[pxdq & 1 == 1] = np.nan
                    median = np.nanmedian(data_temp, axis=(1, 2), keepdims=True)
                    data -= median
                    log.info('  --> Median subtraction: mean of frame median = %.2f' % np.mean(median))
                
                # Write FITS file.
                fitsfile = ut.write_obs(fitsfile, output_dir, data, erro, pxdq, head_pri, head_sci, is2d)
                
                # Update spaceKLIP database.
                self.database.update_obs(key, j, fitsfile)
    
    def subtract_background(self,
                            subdir='bgsub'):
        """ Perform a background subtraction to remove the MIRI glowstick.

        Only required for MIRI.
        Parameters
        ----------
        subdir : str
            Output subdirectory

        Returns None, but writes out new files and updates database

        """
        
        # Set output directory.
        output_dir = os.path.join(self.database.output_dir, subdir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Loop through concatenations.
        for i, key in enumerate(self.database.obs.keys()):
            log.info('--> Concatenation ' + key)
            
            # Find science and reference files.
            ww_sci = np.where(self.database.obs[key]['TYPE'] == 'SCI')[0]
            ww_sci_bg = np.where(self.database.obs[key]['TYPE'] == 'SCI_BG')[0]
            ww_ref = np.where(self.database.obs[key]['TYPE'] == 'REF')[0]
            ww_ref_bg = np.where(self.database.obs[key]['TYPE'] == 'REF_BG')[0]
            
            # Loop through science background files.
            if len(ww_sci_bg) != 0:
                sci_bg_data = []
                sci_bg_erro = []
                sci_bg_pxdq = []
                for j in ww_sci_bg:
                    
                    # Read science background file.
                    fitsfile = self.database.obs[key]['FITSFILE'][j]
                    data, erro, pxdq, head_pri, head_sci, is2d = ut.read_obs(fitsfile)
                    
                    # Compute median science background.
                    sci_bg_data += [data]
                    sci_bg_erro += [erro]
                    sci_bg_pxdq += [pxdq]
                sci_bg_data = np.concatenate(sci_bg_data)
                sci_bg_erro = np.concatenate(sci_bg_erro)
                sci_bg_pxdq = np.concatenate(sci_bg_pxdq)
                sci_bg_data = np.nanmedian(sci_bg_data, axis=0)
                nsample = np.sum(np.logical_not(np.isnan(sci_bg_erro)), axis=0)
                sci_bg_erro = np.true_divide(np.sqrt(np.nansum(sci_bg_erro**2, axis=0)), nsample)
                # if self.database.obs[key]['TELESCOP'][j] == 'JWST' and self.database.obs[key]['INSTRUME'][j] == 'NIRCAM':
                #     sci_bg_pxdq = np.sum(sci_bg_pxdq != 0, axis=0) != 0
                # else:
                sci_bg_pxdq = np.sum(sci_bg_pxdq & 1 == 1, axis=0) != 0
            else:
                sci_bg_data = None
            
            # Loop through reference background files.
            if len(ww_ref_bg) != 0:
                ref_bg_data = []
                ref_bg_erro = []
                ref_bg_pxdq = []
                for j in ww_ref_bg:
                    
                    # Read reference background file.
                    fitsfile = self.database.obs[key]['FITSFILE'][j]
                    data, erro, pxdq, head_pri, head_sci, is2d = ut.read_obs(fitsfile)
                    
                    # Compute median reference background.
                    ref_bg_data += [data]
                    ref_bg_erro += [erro]
                    ref_bg_pxdq += [pxdq]
                ref_bg_data = np.concatenate(ref_bg_data)
                ref_bg_erro = np.concatenate(ref_bg_erro)
                ref_bg_pxdq = np.concatenate(ref_bg_pxdq)
                ref_bg_data = np.nanmedian(ref_bg_data, axis=0)
                nsample = np.sum(np.logical_not(np.isnan(ref_bg_erro)), axis=0)
                ref_bg_erro = np.true_divide(np.sqrt(np.nansum(ref_bg_erro**2, axis=0)), nsample)
                # if self.database.obs[key]['TELESCOP'][j] == 'JWST' and self.database.obs[key]['INSTRUME'][j] == 'NIRCAM':
                #     ref_bg_pxdq = np.sum(ref_bg_pxdq != 0, axis=0) != 0
                # else:
                ref_bg_pxdq = np.sum(ref_bg_pxdq & 1 == 1, axis=0) != 0
            else:
                ref_bg_data = None
            
            # Check input.
            if sci_bg_data is None and ref_bg_data is None:
                raise UserWarning('Could not find any background files')
            
            # Loop through science and reference files.
            for j in np.append(ww_sci, ww_ref):
                if j in ww_sci:
                    sci = True
                else:
                    sci = False
                
                # Read FITS file.
                fitsfile = self.database.obs[key]['FITSFILE'][j]
                data, erro, pxdq, head_pri, head_sci, is2d = ut.read_obs(fitsfile)
                
                # Subtract background.
                head, tail = os.path.split(fitsfile)
                log.info('  --> Background subtraction: ' + tail)
                if sci and sci_bg_data is not None:
                    
                    # test = []
                    # for k in np.logspace(-0.5, 0.5, 100):
                    #     temp = data[0] - k * sci_bg_data
                    #     test += [temp]
                    # test = np.array(test)
                    # hdu0 = pyfits.PrimaryHDU(test)
                    # hdul = pyfits.HDUList([hdu0])
                    # hdul.writeto(os.path.join(output_dir, tail[:-5] + '_test.fits'), output_verify='fix', overwrite=True)
                    # hdul.close()
                    
                    data = data - sci_bg_data
                    erro = np.sqrt(erro**2 + sci_bg_erro**2)
                    # if self.database.obs[key]['TELESCOP'][j] == 'JWST' and self.database.obs[key]['INSTRUME'][j] == 'NIRCAM':
                    #     pxdq[(pxdq == 0) & (sci_bg_pxdq != 0)] += 1
                    # else:
                    pxdq[np.logical_not(pxdq & 1 == 1) & (sci_bg_pxdq != 0)] += 1
                elif sci and sci_bg_data is None:
                    log.warning('  --> Could not find science background, attempting to use reference background')
                    data = data - ref_bg_data
                    erro = np.sqrt(erro**2 + ref_bg_erro**2)
                    # if self.database.obs[key]['TELESCOP'][j] == 'JWST' and self.database.obs[key]['INSTRUME'][j] == 'NIRCAM':
                    #     pxdq[(pxdq == 0) & (ref_bg_pxdq != 0)] += 1
                    # else:
                    pxdq[np.logical_not(pxdq & 1 == 1) & (ref_bg_pxdq != 0)] += 1
                elif not sci and ref_bg_data is not None:
                    data = data - ref_bg_data
                    erro = np.sqrt(erro**2 + ref_bg_erro**2)
                    # if self.database.obs[key]['TELESCOP'][j] == 'JWST' and self.database.obs[key]['INSTRUME'][j] == 'NIRCAM':
                    #     pxdq[(pxdq == 0) & (ref_bg_pxdq != 0)] += 1
                    # else:
                    pxdq[np.logical_not(pxdq & 1 == 1) & (ref_bg_pxdq != 0)] += 1
                elif not sci and ref_bg_data is None:
                    log.warning('  --> Could not find reference background, attempting to use science background')
                    data = data - sci_bg_data
                    erro = np.sqrt(erro**2 + sci_bg_erro**2)
                    # if self.database.obs[key]['TELESCOP'][j] == 'JWST' and self.database.obs[key]['INSTRUME'][j] == 'NIRCAM':
                    #     pxdq[(pxdq == 0) & (sci_bg_pxdq != 0)] += 1
                    # else:
                    pxdq[np.logical_not(pxdq & 1 == 1) & (sci_bg_pxdq != 0)] += 1
                
                # Write FITS file.
                fitsfile = ut.write_obs(fitsfile, output_dir, data, erro, pxdq, head_pri, head_sci, is2d)
                
                # Update spaceKLIP database.
                self.database.update_obs(key, j, fitsfile)
        
        pass
    
    def fix_bad_pixels(self,
                       method='timemed+dqmed+medfilt',
                       bpclean_kwargs={},
                       custom_kwargs={},
                       timemed_kwargs={},
                       dqmed_kwargs={},
                       medfilt_kwargs={},
                       types=['SCI', 'SCI_TA', 'SCI_BG', 'REF', 'REF_TA', 'REF_BG'],
                       subdir='bpcleaned'):
        """Fix bad pixels using custom spaceKLIP routines.

        Multiple routines can be  combined in a custom order by joining them with a + sign.
        - bpclean: use sigma clipping to find additional bad pixels.
        - custom: use custom map to find additional bad pixels.
        - timemed: replace pixels which are only bad in some frames with their
                   median value from the good frames.
        - dqmed:   replace bad pixels with the median of surrounding good
                   pixels.
        - medfilt: replace bad pixels with an image plane median filter.


        Parameters
        ----------
        method : str
        bpclean_kwargs : dict
        custom_kwargs : dict
        timemed_kwargs : dict
        dqmed_kwargs : dict
        medfilt_kwargs : dict
        types : list of str
        subdir : str

        Returns None, but writes out new files to the specified subdirectory, and updates the database object.

        """

        # Set output directory.
        output_dir = os.path.join(self.database.output_dir, subdir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Loop through concatenations.
        for i, key in enumerate(self.database.obs.keys()):
            log.info('--> Concatenation ' + key)
            
            # Loop through FITS files.
            nfitsfiles = len(self.database.obs[key])
            for j in range(nfitsfiles):
                
                # Read FITS file.
                fitsfile = self.database.obs[key]['FITSFILE'][j]
                data, erro, pxdq, head_pri, head_sci, is2d = ut.read_obs(fitsfile)
                
                # Skip file types that are not in the list of types.
                if self.database.obs[key]['TYPE'][j] in types:
                    
                    # Call bad pixel cleaning routines.
                    pxdq_temp = pxdq.copy()
                    # if self.database.obs[key]['TELESCOP'][j] == 'JWST' and self.database.obs[key]['INSTRUME'][j] == 'NIRCAM':
                    #     pxdq_temp = (pxdq_temp != 0) & np.logical_not(pxdq_temp & 512 == 512)
                    # else:
                    pxdq_temp = (np.isnan(data) | (pxdq_temp & 1 == 1)) & np.logical_not(pxdq_temp & 512 == 512)
                    method_split = method.split('+')
                    for k in range(len(method_split)):
                        head, tail = os.path.split(fitsfile)
                        if method_split[k] == 'bpclean':
                            log.info('  --> Method ' + method_split[k] + ': ' + tail)
                            self.find_bad_pixels_bpclean(data, erro, pxdq_temp, pxdq & 512 == 512, bpclean_kwargs)
                        elif method_split[k] == 'custom':
                            log.info('  --> Method ' + method_split[k] + ': ' + tail)
                            self.find_bad_pixels_custom(data, erro, pxdq_temp, key, custom_kwargs)
                        elif method_split[k] == 'timemed':
                            log.info('  --> Method ' + method_split[k] + ': ' + tail)
                            self.fix_bad_pixels_timemed(data, erro, pxdq_temp, timemed_kwargs)
                        elif method_split[k] == 'dqmed':
                            log.info('  --> Method ' + method_split[k] + ': ' + tail)
                            self.fix_bad_pixels_dqmed(data, erro, pxdq_temp, dqmed_kwargs)
                        elif method_split[k] == 'medfilt':
                            log.info('  --> Method ' + method_split[k] + ': ' + tail)
                            self.fix_bad_pixels_medfilt(data, erro, pxdq_temp, medfilt_kwargs)
                        else:
                            log.info('  --> Unknown method ' + method_split[k] + ': skipped')
                    # if self.database.obs[key]['TELESCOP'][j] == 'JWST' and self.database.obs[key]['INSTRUME'][j] == 'NIRCAM':
                    #     pxdq[(pxdq != 0) & np.logical_not(pxdq & 512 == 512) & (pxdq_temp == 0)] = 0
                    # else:
                    # pxdq[(pxdq & 1 == 1) & np.logical_not(pxdq & 512 == 512) & (pxdq_temp == 0)] = 0
                
                # Write FITS file.
                fitsfile = ut.write_obs(fitsfile, output_dir, data, erro, pxdq, head_pri, head_sci, is2d)
                
                # Update spaceKLIP database.
                self.database.update_obs(key, j, fitsfile)
        
        pass
    
    def find_bad_pixels_bpclean(self,
                                data,
                                erro,
                                pxdq,
                                NON_SCIENCE,
                                bpclean_kwargs={}):
        """Implements the 'bpclean' method for find_bad_pixels.


        Parameters
        ----------
        data
        erro
        pxdq
        NON_SCIENCE
        bpclean_kwargs

        Returns
        -------

        """
        
        # Check input.
        if 'sigclip' not in bpclean_kwargs.keys():
            bpclean_kwargs['sigclip'] = 5
        if 'shift_x' not in bpclean_kwargs.keys():
            bpclean_kwargs['shift_x'] = [-1, 0, 1]
        if 'shift_y' not in bpclean_kwargs.keys():
            bpclean_kwargs['shift_y'] = [-1, 0, 1]
        if 0 not in bpclean_kwargs['shift_x']:
            bpclean_kwargs['shift_x'] += [0]
        if 0 not in bpclean_kwargs['shift_y']:
            bpclean_kwargs['shift_y'] += [0]
        
        # Pad data.
        pad_left = np.abs(np.min(bpclean_kwargs['shift_x']))
        pad_right = np.abs(np.max(bpclean_kwargs['shift_x']))
        if pad_right == 0:
            right = None
        else:
            right = -pad_right
        pad_bottom = np.abs(np.min(bpclean_kwargs['shift_y']))
        pad_top = np.abs(np.max(bpclean_kwargs['shift_y']))
        if pad_top == 0:
            top = None
        else:
            top = -pad_top
        pad_vals = ((pad_bottom, pad_top), (pad_left, pad_right))
        
        # Find bad pixels using median of neighbors.
        pxdq_orig = pxdq.copy()
        ww = pxdq != 0
        data_temp = data.copy()
        data_temp[ww] = np.nan
        erro_temp = erro.copy()
        erro_temp[ww] = np.nan
        for i in range(ww.shape[0]):
            
            # Get median background and standard deviation.
            bg_med = np.nanmedian(data_temp[i])
            bg_std = robust.medabsdev(data_temp[i])
            bg_ind = data[i] < (bg_med + 10. * bg_std)  # clip bright PSFs for final calculation
            bg_med = np.nanmedian(data_temp[i][bg_ind])
            bg_std = robust.medabsdev(data_temp[i][bg_ind])
            
            # Create initial mask of large negative values.
            ww[i] = ww[i] | (data[i] < bg_med - bpclean_kwargs['sigclip'] * bg_std)
            
            # Loop through max 10 iterations.
            for it in range(10):
                data_temp[i][ww[i]] = np.nan
                erro_temp[i][ww[i]] = np.nan
                
                # Shift data. 
                pad_data = np.pad(data_temp[i], pad_vals, mode='edge')
                pad_erro = np.pad(erro_temp[i], pad_vals, mode='edge')
                data_arr = []
                erro_arr = []
                for ix in bpclean_kwargs['shift_x']:
                    for iy in bpclean_kwargs['shift_y']:
                        if ix != 0 or iy != 0:
                            data_arr += [np.roll(pad_data, (iy, ix), axis=(0, 1))]
                            erro_arr += [np.roll(pad_erro, (iy, ix), axis=(0, 1))]
                data_arr = np.array(data_arr)
                data_arr = data_arr[:, pad_bottom:top, pad_left:right]
                data_med = np.nanmedian(data_arr, axis=0)
                diff = data[i] - data_med
                data_std = np.nanstd(data_arr, axis=0)
                # data_std = robust.medabsdev(data_arr, axis=0)
                mask_new = diff > bpclean_kwargs['sigclip'] * data_std
                nmask_new = np.sum(mask_new & np.logical_not(ww[i]))
                # print('Iteration %.0f: %.0f bad pixels identified, %.0f are new' % (it + 1, np.sum(mask_new), nmask_new))
                sys.stdout.write('\rFrame %.0f/%.0f, iteration %.0f' % (i + 1, ww.shape[0], it + 1))
                sys.stdout.flush()
                if it > 0 and nmask_new == 0:
                    break
                ww[i] = ww[i] | mask_new
            ww[i][NON_SCIENCE[i]] = 0
            pxdq[i][ww[i]] = 1
        print('')
        log.info('  --> Method bpclean: identified %.0f additional bad pixel(s) -- %.2f%%' % (np.sum(pxdq) - np.sum(pxdq_orig), 100. * (np.sum(pxdq) - np.sum(pxdq_orig)) / np.prod(pxdq.shape)))
        
        pass
    
    def find_bad_pixels_custom(self,
                               data,
                               erro,
                               pxdq,
                               key,
                               custom_kwargs={}):
        """Implements a custom method for find_bad_pixels

        Parameters
        ----------
        data
        erro
        pxdq
        key
        custom_kwargs

        Returns
        -------

        """
        
        # Find bad pixels using median of neighbors.
        pxdq_orig = pxdq.copy()
        pxdq_custom = custom_kwargs[key] != 0
        pxdq_custom = np.array([pxdq_custom] * pxdq.shape[0])
        pxdq[pxdq_custom] = 1
        log.info('  --> Method custom: flagged %.0f additional bad pixel(s) -- %.2f%%' % (np.sum(pxdq) - np.sum(pxdq_orig), 100. * (np.sum(pxdq) - np.sum(pxdq_orig)) / np.prod(pxdq.shape)))
        
        pass
    
    def fix_bad_pixels_timemed(self,
                               data,
                               erro,
                               pxdq,
                               timemed_kwargs={}):
        """Implements the timemed method for find_bad_pixels

        Parameters
        ----------
        data
        erro
        pxdq
        timemed_kwargs

        Returns
        -------

        """
        
        # Fix bad pixels using time median.
        ww = pxdq != 0
        ww_all_bad = np.array([np.sum(ww, axis=0) == ww.shape[0]] * ww.shape[0])
        ww_not_all_bad = ww & np.logical_not(ww_all_bad)
        log.info('  --> Method timemed: fixing %.0f bad pixel(s) -- %.2f%%' % (np.sum(ww_not_all_bad), 100. * np.sum(ww_not_all_bad) / np.prod(ww_not_all_bad.shape)))
        data[ww_not_all_bad] = np.nan
        data[ww_not_all_bad] = np.array([np.nanmedian(data, axis=0)] * data.shape[0])[ww_not_all_bad]
        erro[ww_not_all_bad] = np.nan
        erro[ww_not_all_bad] = np.array([np.nanmedian(erro, axis=0)] * erro.shape[0])[ww_not_all_bad]
        pxdq[ww_not_all_bad] = 0
        
        pass
    
    def fix_bad_pixels_dqmed(self,
                             data,
                             erro,
                             pxdq,
                             dqmed_kwargs={}):
        """ Implements the dqmed method for find_bad_pixels

        Parameters
        ----------
        data
        erro
        pxdq
        dqmed_kwargs

        Returns
        -------

        """

        
        # Check input.
        if 'shift_x' not in dqmed_kwargs.keys():
            dqmed_kwargs['shift_x'] = [-1, 0, 1]
        if 'shift_y' not in dqmed_kwargs.keys():
            dqmed_kwargs['shift_y'] = [-1, 0, 1]
        if 0 not in dqmed_kwargs['shift_x']:
            dqmed_kwargs['shift_x'] += [0]
        if 0 not in dqmed_kwargs['shift_y']:
            dqmed_kwargs['shift_y'] += [0]
        
        # Pad data.
        pad_left = np.abs(np.min(dqmed_kwargs['shift_x']))
        pad_right = np.abs(np.max(dqmed_kwargs['shift_x']))
        if pad_right == 0:
            right = None
        else:
            right = -pad_right
        pad_bottom = np.abs(np.min(dqmed_kwargs['shift_y']))
        pad_top = np.abs(np.max(dqmed_kwargs['shift_y']))
        if pad_top == 0:
            top = None
        else:
            top = -pad_top
        pad_vals = ((0, 0), (pad_bottom, pad_top), (pad_left, pad_right))
        
        # Fix bad pixels using median of neighbors.
        ww = pxdq != 0
        data_temp = data.copy()
        data_temp[ww] = np.nan
        pad_data = np.pad(data_temp, pad_vals, mode='edge')
        erro_temp = erro.copy()
        erro_temp[ww] = np.nan
        pad_erro = np.pad(erro_temp, pad_vals, mode='edge')
        for i in range(ww.shape[0]):
            data_arr = []
            erro_arr = []
            for ix in dqmed_kwargs['shift_x']:
                for iy in dqmed_kwargs['shift_y']:
                    if ix != 0 or iy != 0:
                        data_arr += [np.roll(pad_data[i], (iy, ix), axis=(0, 1))]
                        erro_arr += [np.roll(pad_erro[i], (iy, ix), axis=(0, 1))]
            data_arr = np.array(data_arr)
            data_arr = data_arr[:, pad_bottom:top, pad_left:right]
            data_med = np.nanmedian(data_arr, axis=0)
            ww[i][np.isnan(data_med)] = 0
            data[i][ww[i]] = data_med[ww[i]]
            erro_arr = np.array(erro_arr)
            erro_arr = erro_arr[:, pad_bottom:top, pad_left:right]
            erro_med = np.nanmedian(erro_arr, axis=0)
            erro[i][ww[i]] = erro_med[ww[i]]
            pxdq[i][ww[i]] = 0
        log.info('  --> Method dqmed: fixing %.0f bad pixel(s) -- %.2f%%' % (np.sum(ww), 100. * np.sum(ww) / np.prod(ww.shape)))
        
        pass
    
    def fix_bad_pixels_medfilt(self,
                               data,
                               erro,
                               pxdq,
                               medfilt_kwargs={}):
        """Implements the medfilt method for find_bad_pixels

        Parameters
        ----------
        data
        erro
        pxdq
        medfilt_kwargs

        Returns
        -------

        """
        
        # Check input.
        if 'size' not in medfilt_kwargs.keys():
            medfilt_kwargs['size'] = 4
        
        # Fix bad pixels using median filter.
        ww = pxdq != 0
        log.info('  --> Method medfilt: fixing %.0f bad pixel(s) -- %.2f%%' % (np.sum(ww), 100. * np.sum(ww) / np.prod(ww.shape)))
        data_temp = data.copy()
        data_temp[np.isnan(data_temp)] = 0.
        erro_temp = erro.copy()
        erro_temp[np.isnan(erro_temp)] = 0.
        for i in range(ww.shape[0]):
            data[i][ww[i]] = median_filter(data_temp[i], **medfilt_kwargs)[ww[i]]
            erro[i][ww[i]] = median_filter(erro_temp[i], **medfilt_kwargs)[ww[i]]
            pxdq[i][ww[i]] = 0
        
        pass
    
    def replace_nans(self,
                     cval=0.,
                     types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                     subdir='nanreplaced'):
        """Replace NaN pixels with a constant value.


        Parameters
        ----------
        cval : float
            Constant value to write in place of each NaN
        types : list of str
        subdir : str
            Output subdirectory

        Returns None, but writes out new files to subdir and updates database.


        """
        
        # Set output directory.
        output_dir = os.path.join(self.database.output_dir, subdir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Loop through concatenations.
        for i, key in enumerate(self.database.obs.keys()):
            log.info('--> Concatenation ' + key)
            
            # Loop through FITS files.
            nfitsfiles = len(self.database.obs[key])
            for j in range(nfitsfiles):
                
                # Read FITS file.
                fitsfile = self.database.obs[key]['FITSFILE'][j]
                data, erro, pxdq, head_pri, head_sci, is2d = ut.read_obs(fitsfile)
                
                # Skip file types that are not in the list of types.
                if self.database.obs[key]['TYPE'][j] in types:
                    
                    # Replace nans.
                    head, tail = os.path.split(fitsfile)
                    log.info('  --> Nan replacement: ' + tail)
                    ww = np.isnan(data)
                    data[ww] = cval
                    log.info('  --> Nan replacement: replaced %.0f nan pixel(s) with value ' % (np.sum(ww)) + str(cval) + ' -- %.2f%%' % (100. * np.sum(ww)/np.prod(ww.shape)))
                
                # Write FITS file.
                fitsfile = ut.write_obs(fitsfile, output_dir, data, erro, pxdq, head_pri, head_sci, is2d)
                
                # Update spaceKLIP database.
                self.database.update_obs(key, j, fitsfile)
        
        pass
    
    def align_frames(self,
                     method='fourier',
                     kwargs={},
                     subdir='aligned'):
        """ Use image registration to align all frames in a concatenation to the
        first science frame in that concatenation.


        Parameters
        ----------
        method
        kwargs
        subdir

        Returns None, but writes out new files to subdir and updates database.

        """
        
        # Set output directory.
        output_dir = os.path.join(self.database.output_dir, subdir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Loop through concatenations.
        database_temp = deepcopy(self.database.obs)
        for i, key in enumerate(self.database.obs.keys()):
            log.info('--> Concatenation ' + key)
            
            # Find science and reference files.
            ww_sci = np.where(self.database.obs[key]['TYPE'] == 'SCI')[0]
            if len(ww_sci) == 0:
                raise UserWarning('Could not find any science files')
            ww_ref = np.where(self.database.obs[key]['TYPE'] == 'REF')[0]
            ww_all = np.append(ww_sci, ww_ref)
            
            # Loop through FITS files.
            shifts_all = []
            for j in ww_all:
                
                # Read FITS file.
                fitsfile = self.database.obs[key]['FITSFILE'][j]
                data, erro, pxdq, head_pri, head_sci, is2d = ut.read_obs(fitsfile)
                
                # Align frames.
                head, tail = os.path.split(fitsfile)
                log.info('  --> Align frames: ' + tail)
                if np.sum(np.isnan(data)) != 0:
                    raise UserWarning('Please replace nan pixels before attempting to align frames')
                shifts = []
                mask = None
                for k in range(data.shape[0]):
                    if j == ww_sci[0] and k == 0:
                        ref_image = data[k].copy()
                        pp = np.array([0., 0., 1.])
                        xoffset = self.database.obs[key]['XOFFSET'][j]
                        yoffset = self.database.obs[key]['YOFFSET'][j]
                        crpix1 = self.database.obs[key]['CRPIX1'][j]
                        crpix2 = self.database.obs[key]['CRPIX2'][j]
                    else:
                        p0 = np.array([((crpix1 + xoffset) - (self.database.obs[key]['CRPIX1'][j] + self.database.obs[key]['XOFFSET'][j])) / self.database.obs[key]['PIXSCALE'][j], ((crpix2 + yoffset) - (self.database.obs[key]['CRPIX2'][j] + self.database.obs[key]['YOFFSET'][j])) / self.database.obs[key]['PIXSCALE'][j], 1.])
                        pp = leastsq(self.alignlsq,
                                     p0,
                                     args=(data[k], ref_image, mask, method, kwargs))[0]
                    shifts += [np.array([pp[0], pp[1], pp[2]])]
                    if j != ww_sci[0] or k != 0:
                        data[k] = self.imshift(data[k], [shifts[k][0], shifts[k][1]], method, kwargs)
                        erro[k] = self.imshift(erro[k], [shifts[k][0], shifts[k][1]], method, kwargs)
                shifts = np.array(shifts)
                shifts_all += [shifts]
                
                # Compute shift distances.
                dist = np.sqrt(np.sum(shifts[:, :2]**2, axis=1))  # pix
                dist *= self.database.obs[key]['PIXSCALE'][j]  # mas
                if j == ww_sci[0]:
                    dist = dist[1:]
                log.info('  --> Align frames: median required shift = %.2f mas' % np.median(dist))
                if self.database.obs[key]['TELESCOP'][j] == 'JWST':
                    ww = (dist < 1e-5) | (dist > 100.)
                else:
                    ww = (dist < 1e-5)
                if np.sum(ww) != 0:
                    if j == ww_sci[0]:
                        ww = np.append(np.array([False]), ww)
                    ww = np.where(ww == True)[0]
                    log.warning('  --> The following frames might not be properly aligned: '+str(ww))
                
                # Write FITS file.
                head_pri['XOFFSET'] = xoffset
                head_pri['YOFFSET'] = yoffset
                head_sci['CRPIX1'] = crpix1
                head_sci['CRPIX2'] = crpix2
                fitsfile = ut.write_obs(fitsfile, output_dir, data, erro, pxdq, head_pri, head_sci, is2d)
                
                # Update spaceKLIP database.
                self.database.update_obs(key, j, fitsfile, xoffset=xoffset, yoffset=yoffset, crpix1=crpix1, crpix2=crpix2)
            
            # Plot science frame alignment.
            colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
            f = plt.figure(figsize=(6.4, 4.8))
            ax = plt.gca()
            for index, j in enumerate(ww_sci):
                ax.scatter(shifts_all[index][:, 0] * self.database.obs[key]['PIXSCALE'][j], shifts_all[index][:, 1] * self.database.obs[key]['PIXSCALE'][j], s=5, color=colors[index], marker='o', label='PA = %.0f deg' % self.database.obs[key]['ROLL_REF'][j])
            ax.axhline(0., color='gray', lw=1)
            ax.axvline(0., color='gray', lw=1)
            ax.set_aspect('equal')
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            xrng = xlim[1]-xlim[0]
            yrng = ylim[1]-ylim[0]
            if xrng > yrng:
                ax.set_xlim(np.mean(xlim) - xrng, np.mean(xlim) + xrng)
                ax.set_ylim(np.mean(ylim) - xrng, np.mean(ylim) + xrng)
            else:
                ax.set_xlim(np.mean(xlim) - yrng, np.mean(xlim) + yrng)
                ax.set_ylim(np.mean(ylim) - yrng, np.mean(ylim) + yrng)
            ax.set_xlabel('x-shift [mas]')
            ax.set_ylabel('y-shift [mas]')
            ax.legend(loc='upper right')
            ax.set_title('Science frame alignment')
            output_file = os.path.join(output_dir, key + '_align_sci.pdf')
            plt.savefig(output_file)
            plt.close()
            
            # Plot reference frame alignment.
            colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
            f = plt.figure(figsize=(6.4, 4.8))
            ax = plt.gca()
            seen = []
            reps = []
            syms = ['o', 'v', '^', '<', '>']
            add = len(ww_sci)
            for index, j in enumerate(ww_ref):
                this = '%.0f_%.0f' % (database_temp[key]['XOFFSET'][j], database_temp[key]['YOFFSET'][j])
                if this not in seen:
                    ax.scatter(shifts_all[index + add][:, 0] * self.database.obs[key]['PIXSCALE'][j], shifts_all[index + add][:, 1] * self.database.obs[key]['PIXSCALE'][j], s=5, color=colors[len(seen)], marker=syms[0], label='dpos %.0f' % (len(seen) + 1))
                    ax.hlines(-database_temp[key]['YOFFSET'][j] + yoffset, -database_temp[key]['XOFFSET'][j] + xoffset - 4., -database_temp[key]['XOFFSET'][j] + xoffset + 4., color=colors[len(seen)], lw=1)
                    ax.vlines(-database_temp[key]['XOFFSET'][j] + xoffset, -database_temp[key]['YOFFSET'][j] + yoffset - 4., -database_temp[key]['YOFFSET'][j] + yoffset + 4., color=colors[len(seen)], lw=1)
                    seen += [this]
                    reps += [1]
                else:
                    ww = np.where(np.array(seen) == this)[0][0]
                    ax.scatter(shifts_all[index + add][:, 0] * self.database.obs[key]['PIXSCALE'][j], shifts_all[index + add][:, 1] * self.database.obs[key]['PIXSCALE'][j], s=5, color=colors[ww], marker=syms[reps[ww]])
                    reps[ww] += 1
            ax.set_aspect('equal')
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            xrng = xlim[1]-xlim[0]
            yrng = ylim[1]-ylim[0]
            if xrng > yrng:
                ax.set_xlim(np.mean(xlim) - xrng, np.mean(xlim) + xrng)
                ax.set_ylim(np.mean(ylim) - xrng, np.mean(ylim) + xrng)
            else:
                ax.set_xlim(np.mean(xlim) - yrng, np.mean(xlim) + yrng)
                ax.set_ylim(np.mean(ylim) - yrng, np.mean(ylim) + yrng)
            ax.set_xlabel('x-shift [mas]')
            ax.set_ylabel('y-shift [mas]')
            ax.legend(loc='upper right')
            ax.set_title('Reference frame alignment')
            output_file = os.path.join(output_dir, key + '_align_ref.pdf')
            plt.savefig(output_file)
            plt.close()
        
        pass
    
    def recenter_frames(self,
                        method='fourier',
                        kwargs={},
                        subdir='recentered'):
        """

        Parameters
        ----------
        method
        kwargs
        subdir

        Returns None, but writes out new files to subdir and updates database.

        """
        
        # Set output directory.
        output_dir = os.path.join(self.database.output_dir, subdir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Loop through concatenations.
        for i, key in enumerate(self.database.obs.keys()):
            log.info('--> Concatenation ' + key)
            
            # Find TA files.
            ww_ta = np.append(np.where(self.database.obs[key]['TYPE'] == 'SCI_TA')[0], np.where(self.database.obs[key]['TYPE'] == 'REF_TA')[0])
            
            # Loop through FITS files.
            shifts_all = []
            for j in ww_ta:
                
                # Read FITS file.
                fitsfile = self.database.obs[key]['FITSFILE'][j]
                data, erro, pxdq, head_pri, head_sci, is2d = ut.read_obs(fitsfile)
                
                # Recenter frames.
                head, tail = os.path.split(fitsfile)
                log.info('  --> Recenter frames: ' + tail)
                if np.sum(np.isnan(data)) != 0:
                    raise UserWarning('Please replace nan pixels before attempting to recenter frames')
                shifts = []
                xoffset = 0.  # mas
                yoffset = 0.  # mas
                crpix1 = data.shape[-1]//2 + 1  # 1-indexed
                crpix2 = data.shape[-2]//2 + 1  # 1-indexed
                for k in range(data.shape[0]):
                    p0 = np.array([0., 0.])
                    pp = minimize(self.recenterlsq,
                                  p0,
                                  args=(data[k], method, kwargs))['x']
                    shifts += [np.array([pp[0], pp[1]])]
                    data[k] = self.imshift(data[k], [shifts[k][0], shifts[k][1]], method, kwargs)
                    erro[k] = self.imshift(erro[k], [shifts[k][0], shifts[k][1]], method, kwargs)
                    ww_max = np.unravel_index(np.argmax(data[k]), data[k].shape)
                    if ww_max != (data.shape[-2]//2, data.shape[-1]//2):
                        dx, dy = data.shape[-1]//2 - ww_max[1], data.shape[-2]//2 - ww_max[0]
                        shifts[-1][0] += dx
                        shifts[-1][1] += dy
                        data[k] = np.roll(np.roll(data[k], dx, axis=1), dy, axis=0)
                        erro[k] = np.roll(np.roll(erro[k], dx, axis=1), dy, axis=0)
                shifts = np.array(shifts)
                shifts_all += [shifts]
                
                # Compute shift distances.
                dist = np.sqrt(np.sum(shifts[:, :2]**2, axis=1))  # pix
                dist *= self.database.obs[key]['PIXSCALE'][j]  # mas
                head, tail = os.path.split(self.database.obs[key]['FITSFILE'][j])
                log.info('  --> Align frames: ' + tail)
                log.info('  --> Align frames: median required shift = %.2f mas' % np.median(dist))
                ww = dist > 300.
                if np.sum(ww) != 0:
                    ww = np.where(ww == True)[0]
                    log.warning('The following frames might not be properly aligned (0-indexed): '+str(ww))
                
                # Write FITS file.
                head_pri['XOFFSET'] = xoffset
                head_pri['YOFFSET'] = yoffset
                head_sci['CRPIX1'] = crpix1
                head_sci['CRPIX2'] = crpix2
                fitsfile = ut.write_obs(fitsfile, output_dir, data, erro, pxdq, head_pri, head_sci, is2d)
                
                # Update spaceKLIP database.
                self.database.update_obs(key, j, fitsfile, xoffset=xoffset, yoffset=yoffset, crpix1=crpix1, crpix2=crpix2)
        
        pass
    
    def imshift(self,
                image,
                shift,
                method='fourier',
                kwargs={}):
        """Shift one image

        Parameters
        ----------
        image
        shift : ndarray of length 2
            Shift to apply, in [dX,dY] order
        method : str
            'fourier' or 'spline'
        kwargs : dict
            additional kwargs to pass to the spline_shift function

        Returns shifted image

        """
        
        if method == 'fourier':
            return np.fft.ifftn(fourier_shift(np.fft.fftn(image), shift[::-1])).real
        elif method == 'spline':
            return spline_shift(image, shift[::-1], **kwargs)
        else:
            raise UserWarning('Image shift method "' + method + '" is not known')
    
    def alignlsq(self,
                 shift,
                 image,
                 ref_image,
                 mask=None,
                 method='fourier',
                 kwargs={}):
        """

        Parameters
        ----------
        shift
        image
        ref_image
        mask
        method
        kwargs

        Returns
        -------

        """
        
        if mask is None:
            return (ref_image - shift[2] * self.imshift(image, shift[:2], method, kwargs)).ravel()
        else:
            return ((ref_image - shift[2] * self.imshift(image, shift[:2], method, kwargs)) * mask).ravel()
    
    def recenterlsq(self,
                    shift,
                    image,
                    method='fourier',
                    kwargs={}):
        """

        Parameters
        ----------
        shift
        image
        method
        kwargs

        Returns
        -------

        """
        
        return 1. / np.nanmax(self.imshift(image, shift, method, kwargs))
