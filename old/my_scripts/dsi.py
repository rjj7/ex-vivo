#!/usr/bin/env python
import argparse

import numpy as np
import nibabel as nib
import os
import glob 

from dipy.core.gradients import gradient_table, GradientTable
from dipy.core.ndindex import ndindex
from dipy.data import get_sphere, Sphere, HemiSphere
from dipy.direction import peaks_from_model
from dipy.direction.peaks import peak_directions
from dipy.io import read_bvals_bvecs
from dipy.io.image import save_nifti
from dipy.reconst.dsi import DiffusionSpectrumModel
from dipy.direction.peaks import reshape_peaks_for_visualization
from dipy.segment.mask import median_otsu
from dipy.reconst.shm import sf_to_sh

parser = argparse.ArgumentParser(
        description='Computing Diffusion Spectrum Imaging (DSI)')
parser.add_argument('dwi_data', help='DWI Data')
parser.add_argument('bvecs', help='Bvecs file, in FSL format')
parser.add_argument('bvals', help='Bvals file, in FSL format')
parser.add_argument('output_ODF', help='Output filename ODFs')
parser.add_argument('--mask', dest='mask', default=None, help='Binary mask')
parser.add_argument('--gfa', dest='gfa', help='Save GFA map (Tuch 04)')
parser.add_argument('--rtop', dest='rtop', help='Save rtop-signal map')
parser.add_argument('--save_sh', dest='save_sh', default='',
                    help='Output filename ODF in spherical harmonic coefficients')
parser.add_argument('--sh_order', help='Order of the shperical harmonic (lmax term)', default=6)
parser.add_argument('--output_peaks', help='Output filename peaks reshaped for visualization in FSL or FS')
parser.add_argument('--output_peaks_dirs', help='Output filename peaks directions')
parser.add_argument('--output_peaks_indices', help='Output filename peaks indices')
parser.add_argument('--output_peaks_values', help='Output filenam peaks values')
parser.add_argument('filter_width', default='', help='Width of the Hannin Window')
parser.add_argument('--npeaks', dest='npeaks', default=3, action='store', help='Number of peaks to be selected, int. Default=3')
parser.add_argument('--rel_peak_thr', dest='rel_peak_thr', type=float, action='store', default=.5, help='Relative peak threshold, float. Default=.5')
parser.add_argument('--min_sep_angle', dest='min_sep_angle', default=25, type=int, action='store', help='Minimum separation angle, int. Default=25')
parser.add_argument('--r_start', dest='r_start', type=float, action='store', default=2.1, help='Radial start point of ODF Sampling. Default=2.1')
parser.add_argument('--r_end', dest='r_end', type=float, action='store', default=6.0,help='Radial end point of ODF sampling. Default=6')
parser.add_argument('--r_step', dest='r_step', type=float, action='store', default=0.1,help='Step size of ODF sampling. Default=.1')
parser.add_argument('--qgrid_size', dest='qgrid_size', type=int, action='store', help='Sets the size of the q_space grid. Default=17')
args = parser.parse_args()

if not args.filter_width:
	filter_width = 32
else:
	filter_width = np.asarray(args.filter_width, dtype='int')

if args.npeaks:
	npeaks = np.asarray(args.npeaks, dtype='int')
else:
	npeaks = 3
if args.rel_peak_thr:
	relative_peak_threshold=np.asarray(args.rel_peak_thr, dtype='float')
else:
	relative_peak_threshold = 0.5
if args.min_sep_angle:
	min_separation_angle=np.asarray(args.min_sep_angle, dtype='int')
else:
	min_separation_angle = 25
if args.qgrid_size:
	qgrid_size = np.asarray(args.qgrid_size, dtype='int')
else:
	qgrid_size = 17
if args.r_start:
	r_start = np.asarray(args.r_start, dtype='float')
else:
	r_start = 2.1
if args.r_end:
	r_end = np.asarray(args.r_end, dtype='float')
else:
	r_end = 6.0
if args.r_step:
	r_step = np.asarray(args.r_step, dtype='float')
else:
	r_step = 0.2


# Import DWI data, bvec, bval, and binary mask
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
sphere = get_sphere('symmetric724')

#DSI
print ("Starting model fitting...")
#qgrid_size = 17  #add as option
print(" - filter width: " + str(filter_width))
print(" - q grid size: " + str(qgrid_size))

dsmodel = DiffusionSpectrumModel(gtab, filter_width=filter_width, qgrid_size=qgrid_size)
dsfit = dsmodel.fit(data, mask)
odfs = dsfit.odf(sphere)
#GFA map
if args.gfa:
	from dipy.reconst.odf import gfa
	GFA = gfa(odfs)
	nib.save(nib.Nifti1Image(GFA, img.affine),
                 args.gfa)
#rtop signal
if args.rtop:
	data_norm = data / (data[..., 0, None]).astype(np.float)
	rtop = dsmodel.fit(data_norm).rtop_signal()
	nib.save(nib.Nifti1Image(rtop, img.affine),
                 args.rtop)

#save non normalized as a test
nib.save(nib.Nifti1Image(odfs.astype(np.float32),
                                   img.affine, img.header), 'odfs_dsi_non_norm.nii.gz')
#min/max normalization
ijk = np.ascontiguousarray(np.array(np.nonzero(mask)).T)
shape = data.shape[:-1]
odfs_norm = np.zeros ((shape + (len(sphere.vertices),)))

for (k, center) in enumerate(ijk):
    m = odfs[tuple(center.astype(np.int))].copy()
    m = m - (np.abs(m).min())
    m = m / (np.abs(m).max())
    odfs_norm[tuple(center.astype(np.int))] = m

odfs_norm =odfs_norm.astype(float)
print ("Saving ODFs...")
nib.save(nib.Nifti1Image(odfs_norm.astype(np.float32),
                                   img.affine, img.header), args.output_ODF)

print("Finding peaks...")
peak_dirs = np.zeros (list(data.shape[:-1]) + [npeaks.astype('int'), 3])
peak_values = np.zeros ((shape + (npeaks,)))
peak_indices = np.zeros((shape + (npeaks,)), dtype='int')

print(" - relative_peak_threshold: " + str(relative_peak_threshold))
print(" - min_separation_angle: " + str(min_separation_angle))

for idx in ndindex(shape):
	if mask[idx]:
		direction, pk, indices = peak_directions(odfs_norm[idx], sphere, relative_peak_threshold.astype('float'), min_separation_angle.astype('int'))
		n = min(npeaks, pk.shape[0])
		peak_dirs[idx][:n] = direction[:n]
		peak_values[idx][:n] = pk[:n]
		peak_indices[idx][:n] = indices[:n]
if args.output_peaks:
	print ("Saving peaks...")
	nib.save(nib.Nifti1Image(reshape_peaks_for_visualization(peak_dirs),
                                     img.affine),args.output_peaks)
if args.output_peaks_indices:
	nib.save(nib.Nifti1Image(peak_indices, img.affine),
                 args.output_peaks_indices)

if args.output_peaks_values:
        nib.save(nib.Nifti1Image(peak_values, img.affine),
                 args.output_peaks_values)

if args.output_peaks_dirs:
        nib.save(nib.Nifti1Image(peak_dirs, img.affine),
                 args.output_peaks_dirs)
if args.save_sh:
	print ("Saving ODFs in SH basis...")
	if args.sh_order:
        	sh_order=args.sh_order
	else:
		sh_order=6
	odf_sh = sf_to_sh(odfs_norm, sphere, sh_order=sh_order)
	nib.save(nib.Nifti1Image(odf_sh.astype(np.float32),
                                  img.affine, img.header), args.save_sh)


