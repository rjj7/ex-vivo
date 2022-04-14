#!/usr/bin/env python

import argparse

import os
import nibabel as nb
import numpy as np 

from dipy.core.gradients import gradient_table
from dipy.io.gradients import read_bvals_bvecs
from dipy.reconst.dti import (TensorModel, color_fa, fractional_anisotropy)

parser = argparse.ArgumentParser(
        description='Computes a nifti-compliant directional-encoded color' +
                    'FA image that can be loaded in TrackVis.')
parser.add_argument('dwi_data', metavar='input', help='Diffusion Imaging data')
parser.add_argument('bvals', metavar='bval file', help='bvals file, in FSL format.')
parser.add_argument('bvecs', metavar='bvec file', help='bvecs file, in FSL format.')
parser.add_argument('mask', metavar='mask',
	             help='Binary mask. Only data inside the mask will be'
                          'used for computations and reconstruction.')
parser.add_argument('rgb_trackvis', metavar='output', 
                    help='rgb output image in trackvis format. This is a 3D image where each voxel' +
	                 'contains a tuple of 3 elements, one for each value.')
parser.add_argument('--fa', dest='fa', metavar='file', default='',
        	    help='If specified, the FA map will be also saved with this name.')
parser.add_argument('--rgb', dest='rgb', metavar='file', default='',
        	    help='If specified, the RGB map will be also saved with this name.' +
                    'This is a 4D file where the fourth dimension encodes the RGB channels.')

args = parser.parse_args()

img = nb.load(args.dwi_data)
data = img.get_data()
affine = img.get_affine()

mask = nb.load(args.mask).get_data().astype(np.bool)

bvals, bvecs = read_bvals_bvecs(args.bvals, args.bvecs)
gtab = gradient_table(bvals, bvecs)

#WLS:weighted least squares
tenmodel = TensorModel(gtab, fit_method='WLS')
tenfit = tenmodel.fit(data, mask)

FA = fractional_anisotropy(tenfit.evals)
FA[np.isnan(FA)] = 0
FA = np.clip(FA, 0, 1)

if args.fa:
	fa_img = nb.Nifti1Image(FA.astype(np.float32), affine)
	nb.save(fa_img, args.fa)

RGB = color_fa(FA, tenfit.evecs)
RGB = np.array(255 * RGB, 'uint8')

if args.rgb:
        rgb_img = nb.Nifti1Image(RGB, affine)
        nb.save(rgb_img, args.rgb)

dest_dtype = np.dtype([('R', 'uint8'), ('G', 'uint8'), ('B', 'uint8')])
out_data = np.zeros(RGB.shape[:3], dtype=dest_dtype)

for ii in np.ndindex(RGB.shape[:3]):
    val = RGB[ii]
    out_data[ii] = (val[0], val[1], val[2])

new_hdr = rgb_img.header
new_hdr['dim'][4] = 1
new_hdr.set_intent(1001, name='Color FA')
new_hdr.set_data_dtype(dest_dtype)

nb.save(nb.Nifti1Image(out_data, affine, new_hdr),
            args.rgb_trackvis)

