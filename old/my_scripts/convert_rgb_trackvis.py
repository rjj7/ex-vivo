#!/usr/bin/env python

import argparse

import os
import nibabel as nb
import numpy as np 

parser = argparse.ArgumentParser(
        description='Creates a nifti-compliant directional-encoded color' +
                    'FA image that can be loaded in TrackVis.')

parser.add_argument('rgb_in', help='rgb input image. The 4th dimension.' +
				   'of this image contains 3 values')
parser.add_argument('rgb_out', help='rgb output image. This is a 3D' + 
				    'image where each voxel contains a' + 
				    'tuple of 3 elements, one for each value.')
args = parser.parse_args()

img_orig = nb.load(args.rgb_in)
rgb_in = img_orig.get_data()


dest_dtype = np.dtype([('R', 'uint8'), ('G', 'uint8'), ('B', 'uint8')])
out_data = np.zeros(rgb_in.shape[:3], dtype=dest_dtype)

for ii in np.ndindex(rgb_in.shape[:3]):
    val = rgb_in[ii]
    out_data[ii] = (val[0], val[1], val[2])

new_hdr = img_orig.header
new_hdr['dim'][4] = 1
new_hdr.set_intent(1001, name='Color FA')
new_hdr.set_data_dtype(dest_dtype)

nb.save(nb.Nifti1Image(out_data, img_orig.get_affine(), new_hdr),
            args.rgb_out)

