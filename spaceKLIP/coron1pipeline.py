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

from jwst.datamodels import dqflags, RampModel
from jwst.pipeline import Detector1Pipeline, Image2Pipeline, Coron3Pipeline

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)


# =============================================================================
# MAIN
# =============================================================================

class Coron1Pipeline_spaceKLIP(Detector1Pipeline):
    """
    The spaceKLIP JWST stage 1 pipeline class.
    """
    
    def __init__(self,
                 **kwargs):
        
        # Initialize Detector1Pipeline class.
        super(Coron1Pipeline_spaceKLIP, self).__init__(**kwargs)
        
        # Set additional step parameters.
        self.saturation.grow_diagonal = False
        self.refpix.nlower = 4
        self.refpix.nupper = 4
        self.refpix.nleft = 4
        self.refpix.nright = 4
        self.refpix.nrow_off = 0
        self.refpix.ncol_off = 0
        self.ramp_fit.save_calibrated_ramp = False
        
        pass
    
    def process(self,
                input):
        """Process input file in stage 1.

        """
        
        # Open input as ramp model.
        input = RampModel(input)
        
        # Process MIR & NIR exposures differently.
        instrument = input.meta.instrument.name
        if instrument == 'MIRI':
            input = self.run_step(self.group_scale, input)
            input = self.run_step(self.dq_init, input)
            input = self.run_step(self.saturation, input)
            input = self.run_step(self.ipc, input)
            input = self.run_step(self.firstframe, input)
            input = self.run_step(self.lastframe, input)
            input = self.run_step(self.reset, input)
            input = self.run_step(self.linearity, input)
            input = self.run_step(self.rscd, input)
            input = self.run_step(self.dark_current, input)
            input = self.do_refpix(input)
        else:
            input = self.run_step(self.group_scale, input)
            input = self.run_step(self.dq_init, input)
            input = self.do_saturation(input)
            input = self.run_step(self.ipc, input)
            input = self.run_step(self.superbias, input)
            input = self.do_refpix(input)
            input = self.run_step(self.linearity, input)
            input = self.run_step(self.persistence, input)
            input = self.run_step(self.dark_current, input)
        input = self.run_step(self.jump, input)
        
        # Save calibrated ramp data.
        if self.ramp_fit.save_calibrated_ramp:
            self.save_model(input, suffix='ramp')
        
        # Run ramp fitting & gain scale correction.
        res = self.run_step(self.ramp_fit, input, save_results=False)
        input, ints_model = (res, None) if self.ramp_fit.skip else res
        if input is None:
            self.ramp_fit.log.info('NoneType returned from ramp fitting. Gain scale correction skipped')
        else:
            self.gain_scale.suffix = 'rate'
            input = self.run_step(self.gain_scale, input, save_results=False)
        if ints_model is not None:
            self.gain_scale.suffix = 'rateints'
            ints_model = self.run_step(self.gain_scale, ints_model, save_results=False)
            if self.save_results:
                self.save_model(ints_model, suffix='rateints')
        
        # Setup output file.
        self.setup_output(input)
        
        return input
    
    def run_step(self,
                 step_obj,
                 input,
                 save_results=None,
                 **kwargs):
        """ Run an individual step?

        # FIXME ought to document why this wrapper function exists - what purpose does it have?

        """
        
        # Check if results shall be saved.
        if step_obj.skip:
            really_save_results = False
        elif save_results is not None:
            really_save_results = save_results
        elif step_obj.save_results:
            really_save_results = True
        else:
            really_save_results = False
        
        # Run step. Don't save results yet.
        step_save_orig = step_obj.save_results
        step_obj.save_results = False
        res = step_obj(input)
        step_obj.save_results = step_save_orig
        
        # Check if group scale correction or gain scale correction were
        # skipped.
        if step_obj is self.group_scale:
            if res.meta.cal_step.group_scale == 'SKIPPED':
                really_save_results = False
        elif step_obj is self.gain_scale:
            if res.meta.cal_step.gain_scale == 'SKIPPED':
                really_save_results = False
        
        # Save results.
        if really_save_results:
            step_obj.output_dir = self.output_dir
            if isinstance(res, (tuple)):
                self.save_model(res[0], suffix=step_obj.suffix+'0')
                self.save_model(res[1], suffix=step_obj.suffix+'1')
            else:
                self.save_model(res, suffix=step_obj.suffix)
        
        return res
    
    def do_saturation(self,
                      input,
                      **kwargs):
        """ Wrapper to run modified version of saturation step

        # FIXME document the reasons for the added steps here

        """
        
        # Save original step parameter.
        npix_grow = self.saturation.n_pix_grow_sat
        
        # Run step with default settings.
        if self.saturation.grow_diagonal or npix_grow == 0:
            return self.run_step(self.saturation, input, **kwargs)
        
        # Run step with 1 fewer growth.
        self.saturation.n_pix_grow_sat = npix_grow - 1
        res = self.run_step(self.saturation, input, **kwargs)
        
        # Get saturated pixels and flag 4 neighbors.
        mask_sat = res.groupdq & dqflags.pixel['SATURATED'] == dqflags.pixel['SATURATED']
        mask_vp1 = np.roll(mask_sat, +1, axis=-2)
        mask_vm1 = np.roll(mask_sat, -1, axis=-2)
        mask_hp1 = np.roll(mask_sat, +1, axis=-1)
        mask_hm1 = np.roll(mask_sat, -1, axis=-1)
        
        # Ignore pixels which rolled around array border.
        mask_vp1[:, :,  0, :] = 0
        mask_vm1[:, :, -1, :] = 0
        mask_hp1[:, :, :,  0] = 0
        mask_hm1[:, :, :, -1] = 0
        
        # Combine saturation masks.
        mask_sat = mask_sat | mask_vp1 | mask_vm1 | mask_hp1 | mask_hm1
        
        # Flag saturated pixels.
        res.groupdq = res.groupdq | (mask_sat * dqflags.pixel['SATURATED'])
        
        # Delete unused masks.
        del mask_vp1, mask_vm1, mask_hp1, mask_hm1
        
        return res
    
    def do_refpix(self,
                  input,
                  **kwargs):
        
        # Check if full frame exposure.
        is_full_frame = 'FULL' in input.meta.subarray.name.upper()
        
        # Get number of custom reference pixel rows & columns.
        nlower = self.refpix.nlower
        nupper = self.refpix.nupper
        nleft = self.refpix.nleft
        nright = self.refpix.nright
        nall = nlower + nupper + nleft + nright
        
        # Run step with default settings if full frame exposure, otherwise run
        # step with custom reference pixel rows & columns.
        if is_full_frame or nall == 0:
            return self.run_step(self.refpix, input, **kwargs)
        else:
            return self.do_pseudo_refpix(input, **kwargs)
    
    def do_pseudo_refpix(self,
                         input,
                         **kwargs):
        
        # Get number of custom reference pixel rows & columns.
        nlower = self.refpix.nlower
        nupper = self.refpix.nupper
        nleft = self.refpix.nleft
        nright = self.refpix.nright
        nrow_off = self.refpix.nrow_off
        ncol_off = self.refpix.ncol_off
        nupper_off = -nrow_off if nrow_off != 0 else None
        nright_off = -ncol_off if ncol_off != 0 else None
        
        # Flag custom reference pixel rows & columns.
        self.refpix.log.info(f'Flagging [{nlower}, {nupper}] references rows at [bottom, top] of array')
        self.refpix.log.info(f'Flagging [{nleft}, {nright}] references columns at [left, right] of array')
        input.pixeldq[nrow_off:nrow_off + nlower, ncol_off:nright_off] = input.pixeldq[nrow_off:nrow_off + nlower, ncol_off:nright_off] | dqflags.pixel['REFERENCE_PIXEL']
        input.pixeldq[-nrow_off - nupper:nupper_off, ncol_off:nright_off] = input.pixeldq[-nrow_off - nupper:nupper_off, ncol_off:nright_off] | dqflags.pixel['REFERENCE_PIXEL']
        input.pixeldq[nrow_off:nupper_off, ncol_off:ncol_off + nleft] = input.pixeldq[nrow_off:nupper_off, ncol_off:ncol_off + nleft] | dqflags.pixel['REFERENCE_PIXEL']
        input.pixeldq[nrow_off:nupper_off, -ncol_off - nright:nright_off] = input.pixeldq[nrow_off:nupper_off, -ncol_off - nright:nright_off] | dqflags.pixel['REFERENCE_PIXEL']
        
        # Save original step parameter.
        use_side_orig = self.refpix.use_side_ref_pixels
        if nleft + nright == 0:
            self.refpix.use_side_ref_pixels = False
        else:
            self.refpix.use_side_ref_pixels = True
        
        # Run step with custom reference pixel rows & columns.
        res = self.run_step(self.refpix, input, **kwargs)
        
        # Apply original step parameter.
        self.refpix.use_side_ref_pixels = use_side_orig
        
        # Unflag custom reference pixel rows & columns.
        self.refpix.log.info('Removing custom reference pixel flags')
        res.pixeldq[nrow_off:nrow_off + nlower, ncol_off:nright_off] = res.pixeldq[nrow_off:nrow_off + nlower, ncol_off:nright_off] & ~dqflags.pixel['REFERENCE_PIXEL']
        res.pixeldq[-nrow_off - nupper:nupper_off, ncol_off:nright_off] = res.pixeldq[-nrow_off - nupper:nupper_off, ncol_off:nright_off] & ~dqflags.pixel['REFERENCE_PIXEL']
        res.pixeldq[nrow_off:nupper_off, ncol_off:ncol_off + nleft] = res.pixeldq[nrow_off:nupper_off, ncol_off:ncol_off + nleft] & ~dqflags.pixel['REFERENCE_PIXEL']
        res.pixeldq[nrow_off:nupper_off, -ncol_off - nright:nright_off] = res.pixeldq[nrow_off:nupper_off, -ncol_off - nright:nright_off] & ~dqflags.pixel['REFERENCE_PIXEL']
        
        return res


def run_obs(database,
            steps={},
            subdir='stage1'):
    """ Run all input files in a database through stage 1 processing

    Parameters
    ----------
    database :
    steps : dict
        Nested dict of parameter overrides per step, in the standard JWST pipeline format.
    subdir : str
        output subdirectory

    Returns None, but writes out files to subdir and updates the database

    """
    
    # Set output directory.
    output_dir = os.path.join(database.output_dir, subdir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Loop through concatenations.
    for i, key in enumerate(database.obs.keys()):
        log.info('--> Concatenation ' + key)
        
        # Loop through FITS files.
        nfitsfiles = len(database.obs[key])
        for j in range(nfitsfiles):
            
            # Skip non-stage 0 files.
            head, tail = os.path.split(database.obs[key]['FITSFILE'][j])
            if database.obs[key]['DATAMODL'][j] != 'STAGE0':
                log.info('  --> Coron1Pipeline: skipping non-stage 0 file ' + tail)
                continue
            log.info('  --> Coron1Pipeline: processing ' + tail)
            
            # Initialize Coron1Pipeline.
            pipeline = Coron1Pipeline_spaceKLIP(output_dir=output_dir)
            pipeline.save_results = True
            
            # Set step parameters.
            for key1 in steps.keys():
                for key2 in steps[key1].keys():
                    setattr(getattr(pipeline, key1), key2, steps[key1][key2])
            
            # Run Coron1Pipeline.
            fitspath = os.path.abspath(database.obs[key]['FITSFILE'][j])
            res = pipeline.run(fitspath)
            if isinstance(res, list):
                res = res[0]
            
            # Update spaceKLIP database.
            fitsfile = os.path.join(output_dir, res.meta.filename)
            if fitsfile.endswith('rate.fits'):
                if os.path.isfile(fitsfile.replace('rate.fits', 'rateints.fits')):
                    fitsfile = fitsfile.replace('rate.fits', 'rateints.fits')
            database.update_obs(key, j, fitsfile)
    
    pass
