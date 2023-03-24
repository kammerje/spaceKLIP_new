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

import copy
import webbpsf

from astropy.table import Table
from astroquery.svo_fps import SvoFps
from jwst.pipeline import Detector1Pipeline, Image2Pipeline, Coron3Pipeline

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

from . import utils

# =============================================================================
# MAIN
# =============================================================================

# Load NIRCam, NIRISS, and MIRI filters from the SVO Filter Profile Service.
# http://svo2.cab.inta-csic.es/theory/fps/
wave_nircam = {}
weff_nircam = {}
filter_list = SvoFps.get_filter_list(facility='JWST', instrument='NIRCAM')
for i in range(len(filter_list)):
    name = filter_list['filterID'][i]
    name = name[name.rfind('.') + 1:]
    wave_nircam[name] = filter_list['WavelengthMean'][i] / 1e4 # micron
    weff_nircam[name] = filter_list['WidthEff'][i] / 1e4 # micron
wave_niriss = {}
weff_niriss = {}
filter_list = SvoFps.get_filter_list(facility='JWST', instrument='NIRISS')
for i in range(len(filter_list)):
    name = filter_list['filterID'][i]
    name = name[name.rfind('.') + 1:]
    wave_niriss[name] = filter_list['WavelengthMean'][i] / 1e4 # micron
    weff_niriss[name] = filter_list['WidthEff'][i] / 1e4 # micron
wave_miri = {}
weff_miri = {}
filter_list = SvoFps.get_filter_list(facility='JWST', instrument='MIRI')
for i in range(len(filter_list)):
    name = filter_list['filterID'][i]
    name = name[name.rfind('.') + 1:]
    wave_miri[name] = filter_list['WavelengthMean'][i] / 1e4 # micron
    weff_miri[name] = filter_list['WidthEff'][i] / 1e4 # micron
wave_miri['FND'] = 13. # micron
weff_miri['FND'] = 10. # micron
del filter_list

class Database():
    """
    The central spaceKLIP database class.
    """
    
    def __init__(self,
                 output_dir):
        
        # Output directory for saving manipulated files.
        self.output_dir = output_dir
        
        # Initialize observations dictionary which contains the individual
        # concatenations.
        self.obs = {}
        
        # Initialize reductions dictionary which contains the individual
        # concatenations.
        self.red = {}
        
        # Print information.
        self.verbose = True
        
        pass
    
    def read_jwst_s012_data(self,
                            datapaths,
                            psflibpaths=None,
                            bgpaths=None):
        
        # Check input.
        if isinstance(datapaths, str):
            datapaths = [datapaths]
        if len(datapaths) == 0:
            raise UserWarning('Could not find any data paths')
        if isinstance(psflibpaths, str):
            psflibpaths = [psflibpaths]
        if isinstance(bgpaths, str):
            bgpaths = [bgpaths]
        if bgpaths is not None:
            for i in range(len(bgpaths)):
                if psflibpaths is not None:
                    if bgpaths[i] not in datapaths and bgpaths[i] not in psflibpaths:
                        raise UserWarning('Background path ' + bgpaths[i] + ' does not occur in data or PSF library paths')
                else:
                    if bgpaths[i] not in datapaths:
                        raise UserWarning('Background path ' + bgpaths[i] + ' does not occur in data paths')
        
        # Read FITS headers.
        DATAMODL = []
        TELESCOP = []
        TARGPROP = []
        TARG_RA = [] # deg
        TARG_DEC = [] # deg
        INSTRUME = []
        DETECTOR = []
        FILTER = []
        CWAVEL = [] # micron
        DWAVEL = [] # micron
        PUPIL = []
        CORONMSK = []
        EXP_TYPE = []
        EXPSTART = [] # MJD
        NINTS = []
        EFFINTTM = [] # s
        IS_PSF = []
        SELFREF = []
        SUBARRAY = []
        NUMDTHPT = []
        XOFFSET = [] # mas
        YOFFSET = [] # mas
        APERNAME = []
        PIXSCALE = [] # mas
        BUNIT = []
        CRPIX1 = [] # pix
        CRPIX2 = [] # pix
        VPARITY = []
        V3I_YANG = [] # deg
        RA_REF = [] # deg
        DEC_REF = [] # deg
        ROLL_REF = [] # deg
        HASH = []
        if psflibpaths is not None:
            allpaths = np.array(datapaths + psflibpaths)
        else:
            allpaths = np.array(datapaths)
        Nallpaths = len(allpaths)
        for i in range(Nallpaths):
            hdul = pyfits.open(allpaths[i])
            head = hdul[0].header
            data = hdul['SCI'].data
            if 'uncal' in allpaths[i]:
                DATAMODL += ['STAGE0']
            elif 'rate' in allpaths[i] or 'rateints' in allpaths[i]:
                DATAMODL += ['STAGE1']
            elif 'cal' in allpaths[i] or 'calints' in allpaths[i]:
                DATAMODL += ['STAGE2']
            else:
                raise UserWarning('File name must contain one of the following: uncal, rate, rateints, cal, calints')
            TELESCOP += [head.get('TELESCOP', 'JWST')]
            TARGPROP += [head.get('TARGPROP', 'UNKNOWN')]
            TARG_RA += [head.get('TARG_RA', np.nan)]
            TARG_DEC += [head.get('TARG_DEC', np.nan)]
            INSTRUME += [head.get('INSTRUME', 'UNKNOWN')]
            DETECTOR += [head.get('DETECTOR', 'UNKNOWN')]
            FILTER += [head['FILTER']]
            if TELESCOP[-1] == 'JWST':
                if INSTRUME[-1] == 'NIRCAM':
                    CWAVEL += [wave_nircam[FILTER[-1]]]
                    DWAVEL += [weff_nircam[FILTER[-1]]]
                elif INSTRUME[-1] == 'NIRISS':
                    CWAVEL += [wave_niriss[FILTER[-1]]]
                    DWAVEL += [weff_niriss[FILTER[-1]]]
                elif INSTRUME[-1] == 'MIRI':
                    CWAVEL += [wave_miri[FILTER[-1]]]
                    DWAVEL += [weff_miri[FILTER[-1]]]
                else:
                    raise UserWarning('Data originates from unknown JWST instrument')
            else:
                raise UserWarning('Data originates from unknown telescope')
            PUPIL += [head.get('PUPIL', 'NONE')]
            CORONMSK += [head.get('CORONMSK', 'NONE')]
            EXP_TYPE += [head.get('EXP_TYPE', 'UNKNOWN')]
            EXPSTART += [head.get('EXPSTART', np.nan)]
            NINTS += [head.get('NINTS', data.shape[0] if data.ndim == 3 else 1)]
            EFFINTTM += [head.get('EFFINTTM', np.nan)]
            IS_PSF += [str(head.get('IS_PSF', 'NONE'))]
            SELFREF += [str(head.get('SELFREF', 'NONE'))]
            SUBARRAY += [head.get('SUBARRAY', 'UNKNOWN')]
            NUMDTHPT += [head.get('NUMDTHPT', 1)]
            XOFFSET += [1e3 * head.get('XOFFSET', 0.)]
            YOFFSET += [1e3 * head.get('YOFFSET', 0.)]
            APERNAME += [head.get('APERNAME', 'UNKNOWN')]
            if TELESCOP[-1] == 'JWST':
                PIXSCALE += [utils.lookup_pixscale(INSTRUME[-1], DETECTOR[-1]) * 1e3]
            else:
                raise UserWarning('Data originates from unknown telescope')
            head = hdul['SCI'].header
            BUNIT += [head.get('BUNIT', 'NONE')]
            CRPIX1 += [head.get('CRPIX1', np.nan)]
            CRPIX2 += [head.get('CRPIX2', np.nan)]
            VPARITY += [head.get('VPARITY', -1)]
            V3I_YANG += [head.get('V3I_YANG', 0.)]
            RA_REF += [head.get('RA_REF', np.nan)]
            DEC_REF += [head.get('DEC_REF', np.nan)]
            ROLL_REF += [head.get('ROLL_REF', 0.)]
            HASH += [TELESCOP[-1] + '_' + INSTRUME[-1] + '_' + DETECTOR[-1] + '_' + FILTER[-1] + '_' + PUPIL[-1] + '_' + CORONMSK[-1] + '_' + SUBARRAY[-1]]
            hdul.close()
        DATAMODL = np.array(DATAMODL)
        TELESCOP = np.array(TELESCOP)
        TARGPROP = np.array(TARGPROP)
        TARG_RA = np.array(TARG_RA)
        TARG_DEC = np.array(TARG_DEC)
        INSTRUME = np.array(INSTRUME)
        DETECTOR = np.array(DETECTOR)
        FILTER = np.array(FILTER)
        CWAVEL = np.array(CWAVEL)
        DWAVEL = np.array(DWAVEL)
        PUPIL = np.array(PUPIL)
        CORONMSK = np.array(CORONMSK)
        EXP_TYPE = np.array(EXP_TYPE)
        EXPSTART = np.array(EXPSTART)
        NINTS = np.array(NINTS)
        EFFINTTM = np.array(EFFINTTM)
        IS_PSF = np.array(IS_PSF)
        SELFREF = np.array(SELFREF)
        SUBARRAY = np.array(SUBARRAY)
        NUMDTHPT = np.array(NUMDTHPT)
        XOFFSET = np.array(XOFFSET)
        YOFFSET = np.array(YOFFSET)
        APERNAME = np.array(APERNAME)
        PIXSCALE = np.array(PIXSCALE)
        BUNIT = np.array(BUNIT)
        CRPIX1 = np.array(CRPIX1)
        CRPIX2 = np.array(CRPIX2)
        VPARITY = np.array(VPARITY)
        V3I_YANG = np.array(V3I_YANG)
        RA_REF = np.array(RA_REF)
        DEC_REF = np.array(DEC_REF)
        ROLL_REF = np.array(ROLL_REF)
        HASH = np.array(HASH)
        
        # Find unique concatenations.
        HASH_unique = np.unique(HASH)
        NHASH_unique = len(HASH_unique)
        
        # Associate TA files with science or reference files.
        ww = []
        for i in range(NHASH_unique):
            if 'MASKRND_NONE' in HASH_unique[i] or 'MASKBAR_NONE' in HASH_unique[i] or 'NONE_MASK1065' in HASH_unique[i] or 'NONE_MASK1140' in HASH_unique[i] or 'NONE_MASK1550' in HASH_unique[i] or 'NONE_MASKLYOT' in HASH_unique[i]:
                ww += [-1]
                HASH_i_split = HASH_unique[i].split('_')
                for j in range(NHASH_unique):
                    HASH_j_split = HASH_unique[j].split('_')
                    if HASH_j_split[0] == HASH_i_split[0] and HASH_j_split[1] == HASH_i_split[1] and HASH_j_split[2] == HASH_i_split[2] and HASH_j_split[4] == HASH_i_split[4] and HASH_j_split[5] != 'NONE' and HASH_j_split[6][-4:] in HASH_i_split[6]:
                        ww[-1] = i
                        break
                if ww[-1] != -1:
                    HASH[HASH == HASH_unique[i]] = HASH_unique[j]
                else:
                    raise UserWarning('Could not associate TA files with science or reference files')
        HASH_unique = np.delete(HASH_unique, ww)
        NHASH_unique = len(HASH_unique)
        
        # Loop through concatenations.
        for i in range(NHASH_unique):
            ww = HASH == HASH_unique[i]
            
            # Find science and reference files.
            if psflibpaths is not None:
                ww_sci = []
                ww_ref = []
                for j in range(len(allpaths[ww])):
                    if allpaths[ww][j] in psflibpaths:
                        ww_ref += [j]
                    else:
                        ww_sci += [j]
                ww_sci = np.array(ww_sci)
                ww_ref = np.array(ww_ref)
            else:
                is_psf = IS_PSF[ww]
                exp_type = EXP_TYPE[ww]
                for j in range(len(exp_type)):
                    if 'TA' in exp_type[j]:
                        is_psf[j] = 'False'
                if 'NONE' not in is_psf:
                    ww_sci = np.where(is_psf == 'False')[0]
                    ww_ref = np.where(is_psf == 'True')[0]
                else:
                    log.warning('  --> Could not find IS_PSF header keyword')
                    numdthpt = NUMDTHPT[ww]
                    numdthpt_unique = np.unique(numdthpt)
                    if len(numdthpt_unique) == 2 and numdthpt_unique[0] == 1:
                        ww_sci = np.where(numdthpt == numdthpt_unique[0])[0]
                        ww_ref = np.where(numdthpt == numdthpt_unique[1])[0]
                    else:
                        log.warning('  --> Could not identify science and reference files based on dither pattern')
                        raise UserWarning('Please use psflibpaths to specify reference files')
            
            # Make Astropy tables for concatenations.
            tab = Table(names=('TYPE',
                               'EXP_TYPE',
                               'DATAMODL',
                               'TELESCOP',
                               'TARGPROP',
                               'TARG_RA',
                               'TARG_DEC',
                               'INSTRUME',
                               'DETECTOR',
                               'FILTER',
                               'CWAVEL',
                               'DWAVEL',
                               'PUPIL',
                               'CORONMSK',
                               'EXPSTART',
                               'NINTS',
                               'EFFINTTM',
                               'SUBARRAY',
                               'NUMDTHPT',
                               'XOFFSET',
                               'YOFFSET',
                               'APERNAME',
                               'PIXSCALE',
                               'BUNIT',
                               'CRPIX1',
                               'CRPIX2',
                               'RA_REF',
                               'DEC_REF',
                               'ROLL_REF',
                               'FITSFILE'),
                        dtype=('object',
                               'object',
                               'object',
                               'object',
                               'object',
                               'float',
                               'float',
                               'object',
                               'object',
                               'object',
                               'float',
                               'float',
                               'object',
                               'object',
                               'float',
                               'int',
                               'float',
                               'object',
                               'int',
                               'float',
                               'float',
                               'object',
                               'float',
                               'object',
                               'float',
                               'float',
                               'float',
                               'float',
                               'float',
                               'object'))
            for j in np.append(ww_sci, ww_ref):
                if j in ww_sci:
                    sci = True
                else:
                    sci = False
                if 'TA' in EXP_TYPE[ww][j]:
                    if sci:
                        tt = 'SCI_TA'
                    else:
                        tt = 'REF_TA'
                elif 'BG' in TARGPROP[ww][j].upper() or 'BACKGROUND' in TARGPROP[ww][j].upper():
                    if sci:
                        tt = 'SCI_BG'
                    else:
                        tt = 'REF_BG'
                else:
                    if sci:
                        tt = 'SCI'
                    else:
                        tt = 'REF'
                if bgpaths is not None:
                    if allpaths[ww][j] in bgpaths:
                        if sci:
                            tt = 'SCI_BG'
                        else:
                            tt = 'REF_BG'
                tab.add_row((tt,
                             EXP_TYPE[ww][j],
                             DATAMODL[ww][j],
                             TELESCOP[ww][j],
                             TARGPROP[ww][j],
                             TARG_RA[ww][j],
                             TARG_DEC[ww][j],
                             INSTRUME[ww][j],
                             DETECTOR[ww][j],
                             FILTER[ww][j],
                             CWAVEL[ww][j],
                             DWAVEL[ww][j],
                             PUPIL[ww][j],
                             CORONMSK[ww][j],
                             EXPSTART[ww][j],
                             NINTS[ww][j],
                             EFFINTTM[ww][j],
                             SUBARRAY[ww][j],
                             NUMDTHPT[ww][j],
                             XOFFSET[ww][j],
                             YOFFSET[ww][j],
                             APERNAME[ww][j],
                             PIXSCALE[ww][j],
                             BUNIT[ww][j],
                             CRPIX1[ww][j],
                             CRPIX2[ww][j],
                             RA_REF[ww][j],
                             DEC_REF[ww][j],
                             ROLL_REF[ww][j] - V3I_YANG[ww][j] * VPARITY[ww][j],
                             allpaths[ww][j]))
            self.obs[HASH_unique[i]] = tab.copy()
            del tab
            
            # Associate background files with science or reference files.
            for j in range(len(self.obs[HASH_unique[i]])):
                if self.obs[HASH_unique[i]]['TYPE'][j] == 'SCI_BG':
                    if self.obs[HASH_unique[i]]['EFFINTTM'][j] not in self.obs[HASH_unique[i]]['EFFINTTM'][self.obs[HASH_unique[i]]['TYPE'] == 'SCI']:
                        if self.obs[HASH_unique[i]]['EFFINTTM'][j] in self.obs[HASH_unique[i]]['EFFINTTM'][self.obs[HASH_unique[i]]['TYPE'] == 'REF']:
                            self.obs[HASH_unique[i]]['TYPE'][j] = 'REF_BG'
                        else:
                            raise UserWarning('Background exposure ' + self.obs[HASH_unique[i]]['FITSFILE'][j] + ' could not be matched with PSF')
                elif self.obs[HASH_unique[i]]['TYPE'][j] == 'REF_BG':
                    if self.obs[HASH_unique[i]]['EFFINTTM'][j] not in self.obs[HASH_unique[i]]['EFFINTTM'][self.obs[HASH_unique[i]]['TYPE'] == 'REF']:
                        if self.obs[HASH_unique[i]]['EFFINTTM'][j] in self.obs[HASH_unique[i]]['EFFINTTM'][self.obs[HASH_unique[i]]['TYPE'] == 'SCI']:
                            self.obs[HASH_unique[i]]['TYPE'][j] = 'SCI_BG'
                        else:
                            raise UserWarning('Background exposure ' + self.obs[HASH_unique[i]]['FITSFILE'][j] + ' could not be matched with PSF')
            
            # Reassociate TA and background files with science or reference files based on target name.
            for j in range(len(self.obs[HASH_unique[i]])):
                if self.obs[HASH_unique[i]]['TYPE'][j] in ['SCI_TA', 'SCI_BG']:
                    targprop = self.obs[HASH_unique[i]]['TARGPROP'][self.obs[HASH_unique[i]]['TYPE'] == 'SCI']
                    ww = np.array([s in self.obs[HASH_unique[i]]['TARGPROP'][j] for s in targprop])
                    if np.sum(ww) == 0:
                        targprop = self.obs[HASH_unique[i]]['TARGPROP'][self.obs[HASH_unique[i]]['TYPE'] == 'REF']
                        ww = np.array([s in self.obs[HASH_unique[i]]['TARGPROP'][j] for s in targprop])
                        if np.sum(ww) != 0:
                            self.obs[HASH_unique[i]]['TYPE'][j] = self.obs[HASH_unique[i]]['TYPE'][j].replace('SCI', 'REF')
                if self.obs[HASH_unique[i]]['TYPE'][j] in ['REF_TA', 'REF_BG']:
                    targprop = self.obs[HASH_unique[i]]['TARGPROP'][self.obs[HASH_unique[i]]['TYPE'] == 'REF']
                    ww = np.array([s in self.obs[HASH_unique[i]]['TARGPROP'][j] for s in targprop])
                    if np.sum(ww) == 0:
                        targprop = self.obs[HASH_unique[i]]['TARGPROP'][self.obs[HASH_unique[i]]['TYPE'] == 'SCI']
                        ww = np.array([s in self.obs[HASH_unique[i]]['TARGPROP'][j] for s in targprop])
                        if np.sum(ww) != 0:
                            self.obs[HASH_unique[i]]['TYPE'][j] = self.obs[HASH_unique[i]]['TYPE'][j].replace('REF', 'SCI')
        
        # Print Astropy tables for concatenations.
        if self.verbose:
            self.print_obs()
        
        pass
    
    def read_jwst_s3_data(self,
                          datapaths):
        
        # Check input.
        if isinstance(datapaths, str):
            datapaths = [datapaths]
        if len(datapaths) == 0:
            raise UserWarning('Could not find any data paths')
        
        # Read FITS headers.
        TYPE = []
        DATAMODL = []
        TELESCOP = []
        TARGPROP = []
        TARG_RA = [] # deg
        TARG_DEC = [] # deg
        INSTRUME = []
        DETECTOR = []
        FILTER = []
        CWAVEL = [] # micron
        DWAVEL = [] # micron
        PUPIL = []
        CORONMSK = []
        EXP_TYPE = []
        EXPSTART = [] # MJD
        NINTS = []
        EFFINTTM = [] # s
        SUBARRAY = []
        APERNAME = []
        PIXSCALE = [] # mas
        MODE = []
        ANNULI = []
        SUBSECTS = []
        KLMODES = []
        BUNIT = []
        CRPIX1 = [] # pix
        CRPIX2 = [] # pix
        HASH = []
        Ndatapaths = len(datapaths)
        for i in range(Ndatapaths):
            hdul = pyfits.open(datapaths[i])
            head = hdul[0].header
            if datapaths[i].endswith('i2d.fits'):
                TYPE += ['CORON3']
            elif datapaths[i].endswith('KLmodes-all.fits'):
                TYPE += ['PYKLIP']
            else:
                raise UserWarning('File must have one of the following endings: i2d.fits, KLmodes-all.fits')
            DATAMODL += ['STAGE3']
            TELESCOP += [head.get('TELESCOP', 'JWST')]
            TARGPROP += [head.get('TARGPROP', 'UNKNOWN')]
            TARG_RA += [head.get('TARG_RA', np.nan)]
            TARG_DEC += [head.get('TARG_DEC', np.nan)]
            INSTRUME += [head.get('INSTRUME', 'UNKNOWN')]
            DETECTOR += [head.get('DETECTOR', 'UNKNOWN')]
            FILTER += [head['FILTER']]
            if TELESCOP[-1] == 'JWST':
                if INSTRUME[-1] == 'NIRCAM':
                    CWAVEL += [wave_nircam[FILTER[-1]]]
                    DWAVEL += [weff_nircam[FILTER[-1]]]
                elif INSTRUME[-1] == 'NIRISS':
                    CWAVEL += [wave_niriss[FILTER[-1]]]
                    DWAVEL += [weff_niriss[FILTER[-1]]]
                elif INSTRUME[-1] == 'MIRI':
                    CWAVEL += [wave_miri[FILTER[-1]]]
                    DWAVEL += [weff_miri[FILTER[-1]]]
                else:
                    raise UserWarning('Data originates from unknown JWST instrument')
            else:
                raise UserWarning('Data originates from unknown telescope')
            PUPIL += [head.get('PUPIL', 'NONE')]
            CORONMSK += [head.get('CORONMSK', 'NONE')]
            EXP_TYPE += [head.get('EXP_TYPE', 'UNKNOWN')]
            EXPSTART += [head.get('EXPSTART', np.nan)]
            NINTS += [head.get('NINTS', 1)]
            EFFINTTM += [head.get('EFFINTTM', np.nan)]
            SUBARRAY += [head.get('SUBARRAY', 'UNKNOWN')]
            APERNAME += [head.get('APERNAME', 'UNKNOWN')]
            if TELESCOP[-1] == 'JWST':
                PIXSCALE += [utils.lookup_pixscale(INSTRUME[-1], DETECTOR[-1]) * 1e3]
            else:
                raise UserWarning('Data originates from unknown telescope')
            if TYPE[-1] == 'CORON3':
                MODE += ['RDI']
                ANNULI += [1]
                SUBSECTS += [1]
                try:
                    KLMODES += [str(head['KLMODE0'])]
                except:
                    log.warning('  --> Could not find KL mode in header, assuming default value of 50')
                    KLMODES += ['50']
            elif TYPE[-1] == 'PYKLIP':
                MODE += [head['MODE']]
                ANNULI += [head['ANNULI']]
                SUBSECTS += [head['SUBSECTS']]
                klmodes = str(head['KLMODE0'])
                j = 1
                while True:
                    try:
                        klmodes += ',' + str(head['KLMODE{0}'.format(j)])
                    except:
                        break
                    j += 1
                KLMODES += [klmodes]
            else:
                raise UserWarning('File must have one of the following types: CORON3, PYKLIP')
            if TYPE[-1] == 'CORON3':
                head = hdul['SCI'].header
            BUNIT += [head.get('BUNIT', 'NONE')]
            CRPIX1 += [head.get('CRPIX1', np.nan)]
            CRPIX2 += [head.get('CRPIX2', np.nan)]
            HASH += [TELESCOP[-1] + '_' + INSTRUME[-1] + '_' + DETECTOR[-1] + '_' + FILTER[-1] + '_' + PUPIL[-1] + '_' + CORONMSK[-1] + '_' + SUBARRAY[-1]]
            hdul.close()
        TYPE = np.array(TYPE)
        DATAMODL = np.array(DATAMODL)
        TELESCOP = np.array(TELESCOP)
        TARGPROP = np.array(TARGPROP)
        TARG_RA = np.array(TARG_RA)
        TARG_DEC = np.array(TARG_DEC)
        INSTRUME = np.array(INSTRUME)
        DETECTOR = np.array(DETECTOR)
        FILTER = np.array(FILTER)
        CWAVEL = np.array(CWAVEL)
        DWAVEL = np.array(DWAVEL)
        PUPIL = np.array(PUPIL)
        CORONMSK = np.array(CORONMSK)
        EXP_TYPE = np.array(EXP_TYPE)
        EXPSTART = np.array(EXPSTART)
        NINTS = np.array(NINTS)
        EFFINTTM = np.array(EFFINTTM)
        SUBARRAY = np.array(SUBARRAY)
        APERNAME = np.array(APERNAME)
        PIXSCALE = np.array(PIXSCALE)
        MODE = np.array(MODE)
        ANNULI = np.array(ANNULI)
        SUBSECTS = np.array(SUBSECTS)
        KLMODES = np.array(KLMODES)
        BUNIT = np.array(BUNIT)
        CRPIX1 = np.array(CRPIX1)
        CRPIX2 = np.array(CRPIX2)
        HASH = np.array(HASH)
        
        # Find unique concatenations.
        HASH_unique = np.unique(HASH)
        NHASH_unique = len(HASH_unique)
        
        # Loop through concatenations.
        for i in range(NHASH_unique):
            ww = np.where(HASH == HASH_unique[i])[0]
            
            # Make Astropy tables for concatenations.
            if HASH_unique[i] not in self.red.keys():
                tab = Table(names=('TYPE',
                                   'EXP_TYPE',
                                   'DATAMODL',
                                   'TELESCOP',
                                   'TARGPROP',
                                   'TARG_RA',
                                   'TARG_DEC',
                                   'INSTRUME',
                                   'DETECTOR',
                                   'FILTER',
                                   'CWAVEL',
                                   'DWAVEL',
                                   'PUPIL',
                                   'CORONMSK',
                                   'EXPSTART',
                                   'NINTS',
                                   'EFFINTTM',
                                   'SUBARRAY',
                                   'APERNAME',
                                   'PIXSCALE',
                                   'MODE',
                                   'ANNULI',
                                   'SUBSECTS',
                                   'KLMODES',
                                   'BUNIT',
                                   'FITSFILE'),
                            dtype=('object',
                                   'object',
                                   'object',
                                   'object',
                                   'object',
                                   'float',
                                   'float',
                                   'object',
                                   'object',
                                   'object',
                                   'float',
                                   'float',
                                   'object',
                                   'object',
                                   'float',
                                   'int',
                                   'float',
                                   'object',
                                   'object',
                                   'float',
                                   'object',
                                   'int',
                                   'int',
                                   'object',
                                   'object',
                                   'object'))
            else:
                tab = self.red[HASH_unique[i]].copy()
            for j in range(len(ww)):
                tab.add_row((TYPE[ww[j]],
                             EXP_TYPE[ww[j]],
                             DATAMODL[ww[j]],
                             TELESCOP[ww[j]],
                             TARGPROP[ww[j]],
                             TARG_RA[ww[j]],
                             TARG_DEC[ww[j]],
                             INSTRUME[ww[j]],
                             DETECTOR[ww[j]],
                             FILTER[ww[j]],
                             CWAVEL[ww[j]],
                             DWAVEL[ww[j]],
                             PUPIL[ww[j]],
                             CORONMSK[ww[j]],
                             EXPSTART[ww[j]],
                             NINTS[ww[j]],
                             EFFINTTM[ww[j]],
                             SUBARRAY[ww[j]],
                             APERNAME[ww[j]],
                             PIXSCALE[ww[j]],
                             MODE[ww[j]],
                             ANNULI[ww[j]],
                             SUBSECTS[ww[j]],
                             KLMODES[ww[j]],
                             BUNIT[ww[j]],
                             datapaths[ww[j]]))
            self.red[HASH_unique[i]] = tab.copy()
            del tab
        
        # Print Astropy tables for concatenations.
        if self.verbose:
            self.print_red()
        
        pass
    
    def print_obs(self,
                  include_fitsfiles=False):
        
        # Print Astropy tables for concatenations.
        log.info('--> Identified %.0f concatenation(s)' % len(self.obs))
        for i, key in enumerate(self.obs.keys()):
            log.info('  --> Concatenation %.0f: ' % (i + 1) + key)
            print_tab = copy.deepcopy(self.obs[key])
            if include_fitsfiles:
                print_tab.remove_columns(['TARG_RA', 'TARG_DEC', 'EXPSTART', 'APERNAME', 'CRPIX1', 'CRPIX2', 'RA_REF', 'DEC_REF'])
            else:
                print_tab.remove_columns(['TARG_RA', 'TARG_DEC', 'EXPSTART', 'APERNAME', 'CRPIX1', 'CRPIX2', 'RA_REF', 'DEC_REF', 'FITSFILE'])
            print_tab['XOFFSET'] = np.round(print_tab['XOFFSET'])
            print_tab['XOFFSET'][print_tab['XOFFSET'] == 0.] = 0.
            print_tab['YOFFSET'] = np.round(print_tab['YOFFSET'])
            print_tab['YOFFSET'][print_tab['YOFFSET'] == 0.] = 0.
            print_tab.pprint()
        
        pass
    
    def print_red(self,
                  include_fitsfiles=False):
        
        # Print Astropy tables for concatenations.
        log.info('--> Identified %.0f concatenation(s)' % len(self.red))
        for i, key in enumerate(self.red.keys()):
            log.info('  --> Concatenation %.0f: ' % (i + 1) + key)
            print_tab = copy.deepcopy(self.red[key])
            if include_fitsfiles:
                print_tab.remove_columns(['TARG_RA', 'TARG_DEC', 'EXPSTART', 'APERNAME'])
            else:
                print_tab.remove_columns(['TARG_RA', 'TARG_DEC', 'EXPSTART', 'APERNAME', 'FITSFILE'])
            print_tab.pprint()
        
        pass
    
    def update_obs(self,
                   key,
                   index,
                   fitsfile,
                   nints=None,
                   effinttm=None,
                   xoffset=None,
                   yoffset=None,
                   crpix1=None,
                   crpix2=None):
        
        # Update spaceKLIP database.
        if 'uncal' in fitsfile:
            DATAMODL = 'STAGE0'
        elif 'rate' in fitsfile or 'rateints' in fitsfile:
            DATAMODL = 'STAGE1'
        elif 'cal' in fitsfile or 'calints' in fitsfile:
            DATAMODL = 'STAGE2'
        else:
            raise UserWarning('File name must contain one of the following: uncal, rate, rateints, cal, calints')
        hdul = pyfits.open(fitsfile)
        self.obs[key]['DATAMODL'][index] = DATAMODL
        if nints is not None:
            self.obs[key]['NINTS'][index] = nints
        if effinttm is not None:
            self.obs[key]['EFFINTTM'][index] = effinttm
        self.obs[key]['BUNIT'][index] = hdul['SCI'].header['BUNIT']
        if xoffset is not None:
            self.obs[key]['XOFFSET'][index] = xoffset
        if yoffset is not None:
            self.obs[key]['YOFFSET'][index] = yoffset
        if crpix1 is not None:
            self.obs[key]['CRPIX1'][index] = crpix1
        if crpix2 is not None:
            self.obs[key]['CRPIX2'][index] = crpix2
        self.obs[key]['FITSFILE'][index] = fitsfile
        hdul.close()
        
        pass

    def summarize(self):
        """Succinctly summarize the contents of a database

        How many files are present at each level of reduction, and what kind (sci, ref, TA), etc
        """

        def short_concat_name(concat_name):
            """Return a shorter, less redundant name for a coronagraphic mode, useful for display"""
            parts = concat_name.split('_')
            return "_".join([parts[1], parts[3], parts[5], ])

        for mode in self.obs:
            print(short_concat_name(mode))

            tab = self.obs[mode]

            for stage in [0, 1, 2, 3]:
                stagetab = tab[tab['DATAMODL'] == f'STAGE{stage}']
                if len(stagetab):
                    nsci = np.sum(stagetab['TYPE'] == 'SCI')
                    npsf = np.sum(stagetab['TYPE'] == 'REF')
                    nta = np.sum((stagetab['TYPE'] == 'SCI_TA') | (stagetab['TYPE'] == 'REF_TA'))

                    summarystr = f"\tSTAGE{stage}: {len(stagetab)} files.\t{nsci} SCI, {npsf} PSF"
                    if nta: summarystr += f", {nta} TA"
                    print(summarystr)
