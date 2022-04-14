#!/bin/tcsh -ef

# This shell script extracts the fiber orientation vectors from dipy "dsi.peaks.nii.gz" files.
# That file has all 3 orientation vectors stacked in the 4th dim 
#     (size 9; 1-3 = peak 1, 4-6 = peak 2, 7-9 = peak 3)
# This splits the one file into 3 separate vecs niftis (3 is assuming you saved max 3 odf peaks per voxel)
#
# Input args = path to peaks nifti file, basename for output files

if ($#argv == 0) then
  echo "USAGE: $0 peaksfile outbase"
  exit 1
endif

set peaksfile = $1
set outbase = $2

fslroi $peaksfile "$outbase".v1.nii.gz 0 3
fslroi $peaksfile "$outbase".v2.nii.gz 3 3
fslroi $peaksfile "$outbase".v3.nii.gz 6 3

