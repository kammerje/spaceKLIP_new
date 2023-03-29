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

from jwst.pipeline import Detector1Pipeline, Image2Pipeline, Coron3Pipeline

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)


# =============================================================================
# MAIN
# =============================================================================

class Coron3Pipeline_spaceKLIP(Coron3Pipeline):
    """
    The spaceKLIP JWST stage 3 pipeline class.

    For now, this is simply identical to the regular Coron3Pipeline !
    """



def run_obs(database,
            steps={},
            subdir='stage3'):
    
    # Set output directory.
    output_dir = os.path.join(database.output_dir, subdir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Loop through concatenations.
    datapaths = []
    for i, key in enumerate(database.obs.keys()):
        log.info('--> Concatenation ' + key)
        
        # Make ASN file.
        log.info('  --> Coron3Pipeline: processing ' + key)
        asnpath = make_asn_file(database, output_dir, key)
        
        # Initialize Coron1Pipeline.
        pipeline = Coron3Pipeline_spaceKLIP(output_dir=output_dir)
        pipeline.save_results = True
        
        # Set step parameters.
        for key1 in steps.keys():
            for key2 in steps[key1].keys():
                setattr(getattr(pipeline, key1), key2, steps[key1][key2])
        
        # Run Coron3Pipeline.
        pipeline.run(asnpath)
        
        # Get reduction path.
        datapath = asnpath.replace('asn.json', 'i2d.fits')
        datapaths += [datapath]
        
        # Update reduction header.
        hdul = pyfits.open(datapath)
        hdul[0].header['KLMODE0'] = pipeline.klip.truncate
        hdul.writeto(datapath, output_verify='fix', overwrite=True)
        hdul.close()
    
    # Read reductions into database.
    database.read_jwst_s3_data(datapaths)
    
    pass


def make_asn_file(database,
                  output_dir,
                  key):
    """ Make an association file, for use with stage 3 pipeline processing.

    Parameters
    ----------
    database
    output_dir
    key

    Returns
    -------

    """
    
    # Find science and reference files.
    ww_sci = np.where(database.obs[key]['TYPE'] == 'SCI')[0]
    if len(ww_sci) == 0:
        raise UserWarning('Concatenation ' + key + ' has no science files')
    ww_ref = np.where(database.obs[key]['TYPE'] == 'REF')[0]
    if len(ww_ref) == 0:
        raise UserWarning('Concatenation ' + key + ' has no reference files')
    
    # Make ASN file.
    pro_id = '00000'
    asn_id = 'c0000'
    tar_id = 't000'
    asnfile = key + '_asn.json'
    asnpath = os.path.join(output_dir, asnfile)
    f = open(asnpath, 'w')
    f.write('{\n')
    f.write('    "asn_type": "coron3",\n')
    f.write('    "asn_rule": "candidate_Asn_Lv3Coron",\n')
    f.write('    "program": "' + pro_id + '",\n')
    f.write('    "asn_id": "' + asn_id + '",\n')
    f.write('    "target": "' + tar_id + '",\n')
    f.write('    "asn_pool": "jw' + pro_id + '_00000000t000000_pool.csv",\n')
    f.write('    "products": [\n')
    f.write('        {\n')
    name = key
    f.write('            "name": "' + name + '",\n')
    f.write('            "members": [\n')
    for i in ww_sci:
        f.write('                {\n')
        f.write('                    "expname": "' + os.path.abspath(database.obs[key]['FITSFILE'][i]) + '",\n')
        f.write('                    "exptype": "science"\n')
        f.write('                },\n')
    for i in ww_ref:
        f.write('                {\n')
        f.write('                    "expname": "' + os.path.abspath(database.obs[key]['FITSFILE'][i]) + '",\n')
        f.write('                    "exptype": "psf"\n')
        f.write('                },\n')
    f.write('            ]\n')
    f.write('        }\n')
    f.write('    ]\n')
    f.write('}\n')
    f.close()
    
    return asnpath
