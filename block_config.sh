# Configuration parameters for exvivo processing block
#
# Processing steps:
#

#change 0 to 1 for steps you wish to run

set dogetdata = 1       # Get orig data, bvec, bval
set doorient  = 1	# Fix orientation in DWI header?
set dodenoise = 0	# Do MC-PCA denoise of the data?
set dodegibb  = 0	# Do Gibbs ringing removal?
set doscale   = 0 	# Rescale image values
set dodrift   = 0	# Do temperature drift correction?
set domask    = 0	# Extract brain mask?
set gradcheck = 0 	# Check bvec orientation?
set doeddy    = 0	# Do eddy-current correction?
set doposteddy = 0 # Do post eddy cleanup? (copy bvals, rotate bvecs)
set dobiascor = 0	# Do bias correction?
set dotensor  = 0	# Do tensor fit (FSL)?
set dogqi     = 0	# Do GQI fit (DIPY)?
set dodsi     = 0	# Do DSI fit (DIPY)?
set doresamp  = 0	# Resamp grid to shell?
set domrtrix  = 0	# Do MRtrix processing?
set dotrack   = 0	# Do tractography?
set docsd     = 0 	# Do CSD fit (DIPY)?
set dodtk     = 0	# Do DTK processing?
set doprebed  = 0	# Do bedpostx pre-processing?
set dobed     = 0 	# Do bedpostx main processing?
set dopostbed = 0	# Do bedpostx post-processing?
set doprobtrk = 0	# Do seeded probabilistic tractography?
set doall     = 1
##################################################

# Manually set these for your scan + data:

############################################################
#### Raw data 
############################################################
set inputdir   = /autofs/space/pontus_001/users/data/4.7T/20220305_203422_I55_slab2_2022_03_05_1_1
set rawdata = dwi.nii.gz
set snum      = 12 #number of first b0 scan (folder number)
set enum      = 34 #number of last scan (folder number)
set realorien   = LAI #actual orientation of the sample
set dwidir   = /autofs/space/nyx_002/users/rjones/I55/slab2_20220305_test_proc #this directory is created to house (nearly) all of the output
# Resampling grad info
set bshell = /space/hemera/1/users/cmaffei/ex_vivo/ismrm_abstract/resampling_grad/bvals_multishell #path to shelled bvalues
set gshell = /space/hemera/1/users/cmaffei/ex_vivo/ismrm_abstract/resampling_grad/bvecs_multishell #path to shelled bvecs


############################################################
#### Eddy (FSL)
############################################################

#Eddy - job submission
set doEddyMLSC = 1  #run eddy on MLSC cluster with GPU
set doEddyPBS = 0  #run eddy on PBS/launchpad/nike/atalanta cluster (no GPU)
set doEddyLocal = 0 #run eddy locally with mulitple CPUs/cores

# vvv Have not been implemented in script yet!!! Just some ideas... although shouldn't need any for ex vivo data?
#Eddy - parameters
set eddyDataIsShelled = 1 #Assume, don't check, that data is shelled (default false)
set eddyNiters = "" #number of iterations (def 5)
set eddyReplaceOutliers = "" #detect + replace outlier slices (def false)
set eddyOutlierPos = "" #detect positive and negative outliers (def false)
set eddySepOffsMov = "" #do NOT attempt to separate field offset from subject movement (def false)
set eddyPostEddyAlign = "" #Do NOT perform a post-eddy alignment of shells (def false)

############################################################
#### Preprocessing parameters
############################################################

#dwidenoise parameters
set extent = 7 #size of the window for denoising. Default 7x7x7
#mrdegibbs parameters
set acq_plane = coronal #acquisition plane
#dwibiascorr parameters
set antsb = [100,3] #[scale(mm), spline], decrease scale if small sample; increase spline if very bad bias

############################################################
#### Mask parameters
############################################################

#binary brain mask thresholds
# (should be optimized depending on data (human block-v-monkey, 94T-v-47T, etc.))
set maskthreshlowb = 0.03
set maskthreshhighb = 0.085
set uselowbmask = 1 # use low b vol to create mask
set usehighbmask = 0 # use high b vol to create mask

#wmmask thresholds
set wmmask_fa_lthr = 0.1
set wmmask_fa_uthr = 0.7

############################################################
#### Modeling parameters
############################################################
# DTI
set dotensorfull = 1 #fit diff tensor using full DWI dataset
set dotensorlow = 0 # use only low bval to fit tensor
set tensorlowbThresh = 20000 #bval threshold for lowb-only data + DT fit
set cleantensorlow = 1 #remove DWI, bvecs, bvals from lowb-only data after dtifit
set dotensorkurt = 1 #output kurtosis maps from dtifit (fsl vers>=6)
set savetensor = 1 #save elements of diffusion tensor
set dotensorWLS = 0 #use weighted least squares for fitting

### vvv REMOVED GQI AND DSI FROM SCRIPT!
# GQI
set difflen = 0.4 #diffusion sampling length for GQI
# DSI (dipy)
set filter_width = 32 #Strength of the hanning filter
set filter_type = none
set qgrid_size = 75 #Sets the size of the q_space grid.
set r_start = 4 #ODF is sampled radially in the PDF. This parameters shows where the sampling should start
set r_end = 10 #Radial endpoint of ODF sampling
set r_step = 0.1 #Step size of the ODf sampling from r_start to r_end
set r_weight = 2
# ODFs (GQI, DSI)
set savesh = 1 #save gqi and dsi ODF in sherical harmonics
set outpeaks = 1 #save peaks info (dir, value, indices)
set savedsimaps = 1 #save rtop + GFA maps
set npeaks = 3 #number of peaks to select for voxel for GQI,CSD,DSI
set peak_thr = 0.1
set peak_sep = 25

# CSD
set init_fa = 0.06
set lmax = 6 #spherical harmoinic order CSD
set gm_vox = gm_voxels.nii.gz #filename for the manually selected gm voxels for manual csd response
set wm_vox = wm_voxels.nii.gz #filename for the manually selected wm voxels for manual csd response

############################################################
#### Tractography parameters
############################################################
set max_angle = 30 #max angle threshold
set cutoff = 0.1
set pmf_thr = 0.7
#set step_size =
set str_number = 5000 #total number of streamlines to retain
set seed_number = 5 #number of seeds per voxel
set fa_thr = 0.1 # Threshold to stop tractography
