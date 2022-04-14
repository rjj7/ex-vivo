#!/usr/bin/env python
import argparse

import numpy as np
import nibabel as nib
import os
import glob 

from dipy.core.gradients import gradient_table, GradientTable
from dipy.data import get_sphere
from dipy.direction import peaks_from_model
from dipy.direction.peaks import peak_directions
from dipy.io import read_bvals_bvecs
from dipy.io.image import save_nifti
from dipy.direction.peaks import reshape_peaks_for_visualization
from dipy.reconst.csdeconv import (ConstrainedSphericalDeconvModel, recursive_response) 

parser = argparse.ArgumentParser(
        description='Computing ODFs using constrained spherical deconvolution model (CSD) (Tournier et al., NeurImage 2007). Response function estimated recursively (Tax et al.,  )')

parser.add_argument('dwi_data', help='DWI Data')
parser.add_argument('bvecs', help='Bvecs file, in FSL format')
parser.add_argument('bvals', help='Bvals file, in FSL format')
parser.add_argument('output_ODF', help='Output filename ODFs, SH DIPY basis')
parser.add_argument('mask', help='Binary mask')
#add save response function + screenshot
#parser.add_argument('response', help='Txt file output respsonse function')
parser.add_argument('--peak_thr', type=float, default=0.01, help='How large the second peak can be relative to the first peak in order to call it a single fiber population')
parser.add_argument('--init_fa', dest='init_fa', type=float, default=0.08, help='FA of the initial fat response function (tensor)')
parser.add_argument('--init_trace', dest='init_trace', type=float, default=0.0021, help='trace of the initial fat response function (tensor)')
parser.add_argument('--iter', type=int, default=10, help='maximum nunber of iterations for calibration')
parser.add_argument('--convergence', dest='conv', default=0.001, type=float, help='maximum relative change of SH coefficients')

parser.add_argument('--lambda_', default=1, type=float, help='weight given to the constrained-positivity regularization part of the deconvolution equation')
parser.add_argument('--tau', default=0.1, type=float, help='threshold controlling the amplitude below which the corresponding fODF is assumed to be zero. Ideally, tau should be set to zero. However, to improve the stability of the algorithm, tau is set to tau*100%% of the mean fODF amplitude (here, 10%% by default)')
parser.add_argument('--save_sh_mrtrix', dest='save_sh_mrtrix', help='Output filename ODF SH in mrtrix basis')
parser.add_argument('--sh_order', type=int, help='Maximal order of the shperical harmonic (lmax term)', default=6)
parser.add_argument('--parallel', default=True, type=bool, help='use parallelization in calibration process')
parser.add_argument('--nbr_proc', type=int, default=4, help='number of processes to use')

parser.add_argument('--output_peaks', help='Output filename peaks reshaped for visualization in FSL or FS')
parser.add_argument('--output_peaks_dirs', help='Output filename peaks directions')
parser.add_argument('--output_peaks_indices', help='Output filename peaks indices')
parser.add_argument('--output_peaks_values', help='Output filenam peaks values')
parser.add_argument('--npeaks', dest='npeaks', type=int, default=3, action='store', help='Number of peaks to be selected, int. Default=3')
parser.add_argument('--rel_peak_thr', dest='rel_peak_thr', type=float, action='store', default=.5, help='Relative peak threshold, float. Default=.5')
parser.add_argument('--min_sep_angle', dest='min_sep_angle', default=25, type=int, action='store', help='Minimum separation angle, int. Default=25')

args = parser.parse_args()


if args.save_sh_mrtrix:
	sh_basis = 'tournier07'
else:
	sh_basis = 'descoteaux07'

print('Importing data...')
# Import DWI data, bvec, bval, and binary mask
img = nib.load(args.dwi_data)
data = img.get_data()
print("Data Shape: " + str(data.shape))
 
# Bval and bvec file information import
bvals, bvecs = read_bvals_bvecs(args.bvals, args.bvecs)
gtab = gradient_table(bvals,bvecs)

mask = nib.load(args.mask).get_data().astype(np.bool)

#load sphere
sphere = get_sphere('symmetric362')

#CSD
print ('Estimating response function...')
response = recursive_response(gtab, data, mask=mask, sh_order=args.sh_order,				peak_thr=args.peak_thr, init_fa=args.init_fa,
			 init_trace=args.init_trace, convergence=args.conv,
			iter=args.iter, parallel=args.parallel, 
			nbr_processes=args.nbr_proc)


print ("Starting model fitting...")
csd = ConstrainedSphericalDeconvModel(gtab, response,
		 sh_order=args.sh_order, lambda_=args.lambda_, tau=args.tau)
odfs = peaks_from_model(model=csd, data=data, sphere=sphere, mask=mask,
		 relative_peak_threshold=args.rel_peak_thr,
		min_separation_angle=args.min_sep_angle, return_sh=True, 
		sh_basis_type= sh_basis, return_odf=False)

print ("Saving ODFs...")
nib.save(nib.Nifti1Image(odfs.shm_coeff.astype(np.float32),
                         img.affine, img.header), args.output_ODF)

if args.output_peaks:
	print ("Saving peaks...")
	nib.save(nib.Nifti1Image(reshape_peaks_for_visualization(odfs.peak_dirs), img.affine), args.output_peaks)
if args.output_peaks_indices:
	print ("Saving peaks indices...")
	nib.save(nib.Nifti1Image(odfs.peak_indices, img.affine), args.output_peaks_indices)
if args.output_peaks_values:
	print ("Saving peaks values...")
	nib.save(nib.Nifti1Image(odfs.peak_values, img.affine), args.output_peaks_values)
if args.output_peaks_dirs:
	print ("Saving peaks dirs...")
	nib.save(nib.Nifti1Image(odfs.peak_dirs, img.affine), args.output_peaks_dirs)
print ('Done!')
