#!/usr/bin/env python
import argparse

import numpy as np
import nibabel as nib
import os
import glob 

from nibabel.streamlines import save as save_trk
from nibabel.streamlines import Tractogram

from dipy.core.gradients import gradient_table, GradientTable
from dipy.core.ndindex import ndindex
from dipy.data import get_sphere, Sphere, HemiSphere
from dipy.direction import peaks_from_model, ProbabilisticDirectionGetter
from dipy.direction.peaks import peak_directions
from dipy.io import read_bvals_bvecs
from dipy.io.image import save_nifti
from dipy.direction.peaks import reshape_peaks_for_visualization
from dipy.segment.mask import median_otsu
from dipy.reconst.shm import sf_to_sh

parser = argparse.ArgumentParser(
        description='Computing Generalized Q-ball Imaging (GQI)')
parser.add_argument('dwi_data', help='DWI Data')
parser.add_argument('bvecs', help='Bvecs file, in FSL format')
parser.add_argument('bvals', help='Bvals file, in FSL format')
parser.add_argument('output_ODF', help='Output filename ODFs')

parser.add_argument('--mask', dest='mask', default=None, help='Binary mask')
parser.add_argument('--save_sh', dest='save_sh', default='',
                    help='Output filename ODF in spherical harmonic coefficients (default lmax=6')
parser.add_argument('--sh_order', help='Spherical harmonic order (Lmax term)', default=6)
parser.add_argument('--output_peaks', dest='output_peaks', default='', help='Output filename peaks reshaped for visualization in FSL and FS')
parser.add_argument('--output_peaks_values', help='Output filename peaks values')
parser.add_argument('--output_peaks_indices', help='Output filename peaks indices')
parser.add_argument('--output_peaks_dirs', help='Output filename peaks dirs')
parser.add_argument('diff_length', default='', help='')
parser.add_argument('--npeaks', type=int, dest='npeaks', default=3, help='Number of peaks to be selected, int. Default=3')
parser.add_argument('--rel_peak_thr', dest='rel_peak_thr', default=.5, help='Relative peak threshold, float. Default=.5', type=float)
parser.add_argument('--min_sep_angle', dest='min_sep_angle', default='25', help='Minimum separation angle, int. Default=25', type=int)

args = parser.parse_args()

if not args.diff_length:
	parser.error('You need to specify the diffusion length value')
else:
	diff_length = np.asarray(args.diff_length, dtype='float64')

if args.npeaks:
	npeaks = args.npeaks
else:
	npeaks = 3
if args.rel_peak_thr:
	relative_peak_threshold=args.rel_peak_thr
else:
	 relative_peak_threshold = 0.5
if args.min_sep_angle:
	min_separation_angle= args.min_sep_angle
else:
	min_separation_angle = 25

# Import DWI data, bvec, bval, and binary mask
print ("Imporing data")
fimg =  (args.dwi_data)
img = nib.load(fimg)
data = img.get_data()
print("Data Shape: " + str(data.shape))
 
# Bval and bvec file information import
bvals, bvecs = read_bvals_bvecs(args.bvals, args.bvecs)
gtab = gradient_table(bvals,bvecs)

# Mask 
if args.mask:
        mask = nib.load(args.mask).get_data().astype(np.bool)
else:
	parser.error('Please specify mask for now') #add median_otsu

#load sphere
sphere = get_sphere('symmetric362')

#GQI
print ("Start model fitting")
l_values = np.sqrt(gtab.bvals * 0.01506)
tmp=np.tile(l_values, (3,1))
gradsT = gtab.bvecs.T
b_vector = gradsT * tmp
b_vector = b_vector.T
gqi_vector = np.real(np.sinc(np.dot(b_vector, sphere.vertices.T)* diff_length/np.pi))
odfs = np.dot(data, gqi_vector)

#min/max normalization
ijk = np.ascontiguousarray(np.array(np.nonzero(mask)).T)
shape = data.shape[:-1]
odfs_norm = np.zeros ((shape + (len(sphere.vertices),)))

for (k, center) in enumerate(ijk):
    m = odfs[tuple(center.astype(np.int))].copy()
    m = m - (np.abs(m).min())
    m = m / (np.abs(m).max())
    odfs_norm[tuple(center.astype(np.int))] = m
    
peak_dirs = np.zeros(list(data.shape[0:3]) + [npeaks, 3])
peak_values = np.zeros ((shape + (npeaks,)))
peak_indices = np.zeros((shape + (npeaks,)), dtype='int')

odfs_norm =odfs_norm.astype(float)
print ("Model fitting done. Saving output...")
nib.save(nib.Nifti1Image(odfs_norm.astype(np.float32),
                                   img.affine, img.header), args.output_ODF)

print ("Extracting peak information...")
for idx in ndindex(shape):
	if mask[idx]:
		direction, pk, indices = peak_directions(odfs_norm[idx], sphere, relative_peak_threshold, min_separation_angle)
		n = min(npeaks, pk.shape[0])
		peak_dirs[idx][:n] = direction[:n]
		peak_values[idx][:n] = pk[:n]
		peak_indices[idx][:n] = indices[:n]
if args.output_peaks:
	print ("Saving peaks...")
	nib.save(nib.Nifti1Image(reshape_peaks_for_visualization(peak_dirs), img.affine), args.output_peaks)
if args.output_peaks_indices:
	print ("Saving peaks indices...")
	nib.save(nib.Nifti1Image(peak_indices, img.affine), args.output_peaks_indices)
if args.output_peaks_values:
        print ("Saving peaks values...")
        nib.save(nib.Nifti1Image(peak_values, img.affine), args.output_peaks_values)
if args.output_peaks_dirs:
        print ("Saving peaks dirs...")
        nib.save(nib.Nifti1Image(peak_dirs, img.affine), args.output_peaks_dirs)
print ("Saving ODFs in SH basis...")
if args.save_sh:
	if args.sh_order:
        	sh_order=args.sh_order
	else:
		sh_order=6
	odfs_sh = sf_to_sh(odfs_norm, sphere, sh_order=sh_order)
	nib.save(nib.Nifti1Image(odfs_sh.astype(np.float32),
                                  img.affine, img.header), args.save_sh)

print ("Done!")
