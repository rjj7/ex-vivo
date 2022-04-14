#!/usr/bin/env python

import argparse

import nibabel as nb
import numpy as np

from dipy.io.image import save_nifti

parser = argparse.ArgumentParser(
        description='Projects one of the vector components to the plane')
parser.add_argument('vector_file', metavar='input', help='4D orientation file')
parser.add_argument('plane', metavar='plane', help='Plane to project on (x, y, or z).')
parser.add_argument('output', metavar='output', 
                    help='Output image.')

args = parser.parse_args()

img = nb.load(args.vector_file)
data = img.get_data()
affine = img.get_affine()

if args.plane=='x':
    for ii in np.ndindex (data.shape[:3]):
        data[ii][0]=0.001
elif args.plane=='y':
    for ii in np.ndindex (data.shape[:3]):
        data[ii][1]=0.001
else:
    for ii in np.ndindex (data.shape[:3]):
        data[ii][2]=0.001

save_nifti(args.output, data, affine)
