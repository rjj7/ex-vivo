#!/usr/bin/env python

import argparse

import nibabel as nb
import numpy as np

parser = argparse.ArgumentParser(
        description='Computes RGB map from OCT dyads.nii.')
parser.add_argument('dyads', metavar='input', help='Dyads from OCT')
parser.add_argument('retardance', metavar='retardance', help='Retardance from OCT')
parser.add_argument('rgb', metavar='output', 
                    help='rgb output image. This is a 4D image where where the' +
                    'fourth dimension encodes the RGB channels. Brightness is'  +
                    'modulated by the retardance')

args = parser.parse_args()

img = nb.load(args.dyads)
dyads = img.get_data()
affine = img.get_affine()

img = nb.load(args.retardance)
ret = img.get_data()

dyads=dyads.astype('float64')
for ii in np.ndindex (dyads.shape[:3]):
    dyads[ii][0]=0.001

rgb=np.abs(dyads)*np.clip(ret, 0, 1)[..., None]
rgb_img = nb.Nifti1Image(np.array(255*rgb, 'uint8'), affine)

rgb_img.header['qform_code']=1
rgb_img.header['sform_code']=1
rgb_img.header['xyzt_units']=18

nb.save(rgb_img, args.rgb)
