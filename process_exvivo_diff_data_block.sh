#! /bin/tcsh -ef
#set echo=1

#This scripts assumes ANACONDA is being added to your $path
#The scripts used API from anaconda 3
#To add this ANACONDA version to yourpath you can run:
#setenv PATH /autofs/space/tinia_001/users/chiara/anaconda3_new/bin:$PATH

set noglob

if ($#argv == 0) then
  echo "USAGE: $0 configfile"
  exit 1
endif

echo "Sourcing configuration from $argv[1]"
source $argv[1]
echo "properly sourced configuration file"
if (! $?inputdir) then
  echo "ERROR: Must specify inputdir in configuration file"
endif
if (! $?snum) then 
  echo "ERROR: Must specify snum in configuration file"
endif
if (! $?enum) then
  echo "ERROR: Must specify enum in configuration file"
endif
if (! $?realorien) then
  echo "ERROR: Must specify realorien in configuration file"
endif
if (! $?dwidir) then
  echo "ERROR: Must specify dwidir in configuration file"
endif
if (! $?dogetdata) then
  set dogetdata = 0
endif
if (! $?dodenoise) then	
  set dodenoise = 0
endif
if (! $?dodrift) then
  set dodrift = 0
endif
if (! $?dodegibb) then 
  set dodegibb = 0
endif
if (! $?doscale) then 
  set doscale = 0
endif 
if (! $?doorient) then	
  set doorient = 0
endif
if (! $?domask) then	
  set domask = 0
endif
if (! $?doeddy) then	
  set doeddy = 0
endif
if (! $?dobiascor) then
  set dobiascor = 0
endif
if (! $?dotensor) then	
  set dotensor = 0
endif
if (! $?dogqi) then
  set dogqi = 0
endif
if (! $?dodsi) then
  set dodsi = 0
endif
if (! $?doresamp) then
  set doresamp = 0 
endif
if (! $?domrtrix) then
  set domrtrix = 0
endif 
if (! $?dodtk) then	
  set dodtk = 0
endif
if (! $?docsd) then 
  set doprobtrk = 0
endif
if (! $?maskthresh) then	
  set maskthresh = 0.01
endif
if (! $?extent) then 
  set extent = 7
endif
if (! $?antsb) then
    set antsb = [100,3]
endif
if (! $?difflen) then 
  set difflen = 0.4
endif
if (! $?filter_width) then
  set filter_width = 32
endif
if (! $?npeaks) then 
  set npeaks = 3
endif
if (! $?qgrid_size) then  
  set qgrid_size = 17
endif
if (! $?r_start) then
  set r_start = 2.1
endif
if (! $?r_end) then 
  set r_end = 6
endif
if (! $?r_step) then 
  set r_step = 0.2
endif
if (! $?str_number) then 
  set str_number = 5000
endif
if (! $?seed_number) then
  set seed_number = 5
endif
if (! $?pmf_thr) then 
set pmf_thr = 0.1
endif
if (! $?fa_init) then 
  set inita_fa = 0.06
endif
if (! $?peak_thr) then 
  set peak_thr = 0.01
endif
if (! $?fa_thr) then 
  set fa_thr = 0.1
endif

echo "Everything from the configuration file has been properly sourced."

# Create output directory
umask 002
mkdir -p $dwidir
set LF = $dwidir/log.txt

# Source Mrtrix
set mrtrix = /usr/pubsw/packages/mrtrix/current/bin/ 
source /usr/pubsw/packages/mrtrix/env.csh
# Set Anaconda
setenv PATH /autofs/space/tinia_001/users/chiara/anaconda3_new/bin:$PATH

set myscriptsdir = /autofs/space/hemera_001/users/cmaffei/my_scripts/ 

echo "Starting processing"
echo "----------------------------" >> $LF
echo "DATE `date`"                  >> $LF
echo "USER `whoami`"                >> $LF
echo "HOST `hostname`"              >> $LF

########################################################
if ($dogetdata) then 
  echo "Extracting bvecs and bvals from method files..."

  set arr = `seq $snum $enum`
  set kstart = 1
  set kend = $#arr

  @ k = $kstart
  printf '' > $dwidir/dwi_orig.bvals
  printf '' > $dwidir/dwi_orig.bvecs
  printf '' > $dwidir/params.txt

  while ($k <= $kend)
        set fn = $arr[$k]
        set bvalStr = `cat $inputdir/$fn/method | awk '$1 ~ /EffBval/ { getline; do {print $0; getline} while ($0 !~ /=/) }'`
        set bvecStr = `cat $inputdir/$fn/method | awk '$1 ~ /PVM_DwDir=/ {getline; do {print $0; getline} while ($0 !~ /=/) }'`
        set epiNEchoesStr = `cat $inputdir/$fn/method | awk '$1 ~ /PVM_EpiNEchoes=/ { do {print $0; getline} while ($0 !~ /=/) }'`
        set epiSpacingStr = `cat $inputdir/$fn/method | awk '$1 ~ /PVM_EpiEchoSpacing=/ { do {print $0; getline} while ($0 !~ /=/) }'`
        printf '%s %s\n' $epiNEchoesStr $epiSpacingStr >> $dwidir/params.txt
        printf '%s\n' $bvalStr >> $dwidir/dwi_orig.bvals
        printf '%s %s %s\n' $bvecStr >> $dwidir/dwi_orig.bvecs
        @ k = $k + 1
  end

  # Get values for first three columns of acqparams.txt file - based on phase encoding direction

  set phaseEncStr = `cat $inputdir/$enum/method | awk '$1 ~ /PVM_SPackArrReadOrient=/ {getline; do {print $0; getline} while ($0 !~ /=/) }'`

  if ($phaseEncStr == "L_R") then
        set phaseEncVal = "-1 0 0"
  else if ($phaseEncStr == "R_L") then
        set phaseEncVal = "1 0 0"
  else if ($phaseEncStr == "A_P") then
        set phaseEncVal = "0 -1 0"
  else if ($phaseEncStr == "P_A") then
        set phaseEncVal = "0 1 0"
  else if ($phaseEncStr == "H_F" || $phaseEncStr == "F_H") then
        set phaseEncVal = "0 0 -1"
  endif

  # Get epiNEchoes and epiEchoSpacing, use them to compute fourth column of acqparams.txt file
  set epiNEchoes = `cat $inputdir/$enum/method | awk '$1 ~ /PVM_EpiNEchoes=/ {print $0}'`
  set epiNEchoesVal = `printf "$epiNEchoes" | awk '{if ((nex=index($1,"=")) > 0) {nenum = substr($1,nex+1,length($1)); printf("%s", nenum)}}'`
  echo ".........."

  set epiSpacing = `cat $inputdir/$enum/method | awk '$1 ~ /PVM_EpiEchoSpacing=/ {print $0}'`
  set epiSpacingVal = `printf "$epiSpacing" | awk '{if ((nes=index($1,"=")) > 0) {esnum = substr($1,nes+1,length($1)); printf("%s", esnum)}}'`
  echo ".........."

  set acqpVal = `echo "($epiNEchoesVal - 1)*$epiSpacingVal*0.001" | bc -l`
  set acqpVal = `printf "%.7f" $acqpVal`
  echo ".........."

  # Get total number of lines/number of DW volumes
  set nbLines = `wc -l < $dwidir/dwi_orig.bvals`
  echo ".........."

  # Write to acqparams.txt file
  @ i = 1
  printf '' > $dwidir/acqparams.txt
  while ($i <= $nbLines)
    printf '%s %s %s %s\n' $phaseEncVal $acqpVal >> $dwidir/acqparams.txt
    @ i = $i + 1
  end
  echo ".........."
  # Write to index.txt file
     @ i = 1
     printf '' > $dwidir/index.txt
     while ($i <= $nbLines)
       printf '%s\n' $i >> $dwidir/index.txt
       @ i = $i + 1
     end

  # Copy data from from raw data directory
  set cmd = (mri_convert $inputdir/$rawdata $dwidir/dwi_orig.nii.gz)
  echo $cmd | tee -a $LF
  $cmd | tee -a $LF
endif 
##################################################

if ($doorient) then
  set cmd = mri_convert
  set cmd = ($cmd --in_orientation $realorien)
  set cmd = ($cmd $dwidir/dwi_orig.nii.gz)
  set cmd = ($cmd $dwidir/dwi_orig.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  #  Convert DWIs and gradient vectors to LAS orientation
  #  (Use this instead of flip4fsl, which inverts bvecs in x erroneously)
  #  Note, this will also produce the reoriented bvec file 
  #  NB: bvec file needs to be consistently named with corresponding dwi!!!
  set cmd = orientLAS
  set cmd = ($cmd $dwidir/dwi_orig.nii.gz)
  set cmd = ($cmd $dwidir/dwi_las.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
endif

if ($doall) then 
if ($dodenoise) then 
  #
  # Apply denoising Veraart 2016, MRtrix Implementation
  echo "Denoising data..."
  # 
  if (-e $dwidir/dwi_las.nii.gz) then
  set dwiname = dwi_las
  else
  set dwiname = dwi_orig
  endif
  set cmd = $mrtrix/dwidenoise
  set cmd = ($cmd -extent $extent)
  set cmd = ($cmd -noise $dwidir/noise_est.nii.gz)
  set cmd = ($cmd $dwidir/$dwiname.nii.gz)
  set cmd = ($cmd $dwidir/dwi_denoised.nii.gz)
  set cmd = ($cmd -force)
  echo $cmd | tee -a $LF
  $cmd |& tee -a $LF

  # Compute residuals
  set cmd = $mrtrix/mrcalc
  set cmd = ($cmd $dwidir/$dwiname.nii.gz)
  set cmd = ($cmd $dwidir/noise_est.nii.gz)
  set cmd = ($cmd -subtract)
  set cmd = ($cmd $dwidir/residuals.nii.gz)
  set cmd = ($cmd -quiet)
  set cmd = ($cmd -force)
  echo $cmd | tee -a $LF
  $cmd |& tee -a $LF
  # Compute RMS
  set cmd = $mrtrix/mrmath
  set cmd = ($cmd $dwidir/residuals.nii.gz)
  set cmd = ($cmd rms)
  set cmd = ($cmd -axis 3)
  set cmd = ($cmd $dwidir/rms_residuals.nii.gz)
  set cmd = ($cmd -quiet)
  set cmd = ($cmd -force)
  echo $cmd | tee -a $LF
  $cmd |& tee -a $LF
endif
endif

if ($doall) then
if ($dodegibb) then
  if (-e $dwidir/dwi_denoised.nii.gz) then
     set dwiname = dwi_denoised.nii.gz
  else if (-e $dwidir/dwi_las.nii.gz) then 
     set dwiname = dwi_las.nii.gz
  else
     set dwiname = dwi_orig.nii.gz
  endif
  #Specify acquisition plane in -axes option:
  # 0,1 axial 0,2 coronal 1,2 sagittal
  if ($?acq_plane) then
  if ($acq_plane == axial ) then 
     set axes = 0,1
  else if ($acq_plane == coronal) then 
     set axes = 0,2
  else 
     set exes = 1,2
  endif
  else 
  echo "Acq plane non specified. Coronal will be used."
  set axes = 0,2
  endif
  #
  # Attempts to remove Gibbs ringing artefacts from MR
  # Kellner et al 2016, MRtrix implementation
  # 
  set cmd = $mrtrix/mrdegibbs
  set cmd = ($cmd -axes $axes)
  set cmd = ($cmd $dwidir/$dwiname)
  set cmd = ($cmd $dwidir/dwi_degibb.nii.gz) 
  set cmd = ($cmd -force) 
  echo $cmd | tee -a $LF
  $cmd |& tee -a $LF
endif
endif 

if ($doall) then
if ($doscale) then 
  if (-e $dwidir/dwi_degibb.nii.gz) then 
     set dwiname = dwi_degibb.nii.gz 
  else if (-e $dwidir/dwi_denoised.nii.gz) then
     set dwiname = dwi_denoised.nii.gz
  else if  (-e $dwidir/dwi_las.nii.gz) then
     set dwiname = dwi_las.nii.gz
  else 
     set dwiname = dwi_orig.nii.gz
  endif
  set min = `fslstats $dwidir/$dwiname -R | awk '{ print $1 }'`
  set max = `fslstats $dwidir/$dwiname -R | awk '{ print $2 }'`
  set scaling = `echo "scale=5; 32767.0 / ( $max - $min )" | bc`
  fslmaths $dwidir/$dwiname -sub $min -mul $scaling $dwidir/$dwiname
endif
endif

if ($doall) then
if ($dodrift) then
  if (-e $dwidir/dwi_degibb.nii.gz) then
     set dwiname = dwi_degibb.nii.gz
  else if  (-e $dwidir/dwi_denoised.nii.gz) then
     set dwiname = dwi_denoised.nii.gz
  else if (-e $dwidir/dwi_las.nii.gz) then
     set dwiname = dwi_las.nii.gz
  else
     set dwiname = dwi_orig.nii.gz
  endif
  #
  # Normalize images to compensate for temperature drift
  # 
  set lowb = 1
  set cmd = "addpath $FREESURFER_HOME/matlab"
  set cmd = "addpath /autofs/space/turan_001/users/lzollei/dev/matlab"
  set cmd = "$cmd; fix_exvivo_dwi_drift("
  set cmd = "$cmd '$dwidir/dwi_drift.nii.gz', "
  set cmd = "$cmd '$dwidir/drift', "
  set cmd = "$cmd '$dwidir/$dwiname', "
  set cmd = "$cmd $#lowb);"
  echo $cmd | tee -a $LF
  echo $cmd | matlab -nosplash
endif
endif

if ($doall) then
if ($domask) then
  if (-e $dwidir/dwi_degibb.nii.gz) then
     set dwiname = dwi_degibb
  else if  (-e $dwidir/dwi_denoised.nii.gz) then
     set dwiname = dwi_denoised
  else if (-e $dwidir/dwi_las.nii.gz) then
     set dwiname = dwi_las
  else
     set dwiname = dwi_orig
  endif
  set cmd = fslroi
  set cmd = ($cmd $dwidir/$dwiname.nii.gz)
  set cmd = ($cmd $dwidir/${dwiname}_lowb.nii.gz)
  set cmd = ($cmd 0 1)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set dmax = `fslstats $dwidir/${dwiname}_lowb.nii.gz -R | awk '{print $2}'`
  set dthresh = `echo "$dmax * $maskthreshlowb" | bc`
    
  set cmd = mri_binarize
  set cmd = ($cmd --i $dwidir/${dwiname}_lowb.nii.gz)
  set cmd = ($cmd --min $dthresh)
  set cmd = ($cmd --o $dwidir/lowb_mask.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set cmd = $mrtrix/maskfilter
  set cmd = ($cmd $dwidir/lowb_mask.nii.gz)
  set cmd = ($cmd clean)
  set cmd = ($cmd -force)
  set cmd = ($cmd $dwidir/lowb_mask.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set cmd = $mrtrix/maskfilter
  set cmd = ($cmd $dwidir/lowb_mask.nii.gz)
  set cmd = ($cmd dilate)
  set cmd = ($cmd -npass 6 -force)
  set cmd = ($cmd $dwidir/lowb_mask.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
    
  set cmd = fslroi
  set cmd = ($cmd $dwidir/$dwiname.nii.gz)
  set cmd = ($cmd $dwidir/${dwiname}_highb.nii.gz)
  set cmd = ($cmd 1 2)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF 
    
  set dmax = `fslstats $dwidir/${dwiname}_highb.nii.gz -R | awk '{print $2}'`
  set dthresh = `echo "$dmax * $maskthreshhighb" | bc`
    
  set cmd = mri_binarize
  set cmd = ($cmd --i $dwidir/${dwiname}_highb.nii.gz)
  set cmd = ($cmd --min $dthresh)
  set cmd = ($cmd --o $dwidir/highb_mask.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
  set cmd = $mrtrix/maskfilter
  set cmd = ($cmd $dwidir/highb_mask.nii.gz)
  set cmd = ($cmd erode)
  set cmd = ($cmd -npass 9 -force)
  set cmd = ($cmd $dwidir/highb_mask.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
  set cmd = $mrtrix/maskfilter
  set cmd = ($cmd $dwidir/highb_mask.nii.gz)
  set cmd = ($cmd dilate)
  set cmd = ($cmd -npass 10 -force)
  set cmd = ($cmd $dwidir/highb_mask.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
endif
endif 

if ($doall) then
if ($gradcheck) then
  #checking gradients
  echo "Checking gradients..."
  echo "-------------------------------" >> $LF
  echo "------Checking Gradients-------" >> $LF
  if (-e $dwidir/dwi_las.nii.gz) then
     set dwiname = dwi_las
  else
     set dwiname = dwi_orig
  endif
  # Transpose bvec/bval files to make mrtrix happy 
  set cmd = /autofs/space/hemera_001/users/cmaffei/ex_vivo/proc.5.2020/code/transpose_gradients.py
  set cmd = (python $cmd)
  set cmd = ($cmd $dwidir/$dwiname.bvecs)
  set cmd = ($cmd $dwidir/${dwiname}_t.bvecs)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
  
  # Replace first bval with 0 to make Mrtrix happy
  set cmd = (sed -i "1s/.*/0/" $dwiname.bvals)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
  #
  set cmd = /autofs/space/hemera_001/users/cmaffei/ex_vivo/proc.5.2020/code/transpose_gradients.py
  set cmd = (python $cmd)
  set cmd = ($cmd $dwidir/$dwiname.bvals)
  set cmd = ($cmd $dwidir/${dwiname}_t.bvals)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
  
  set cmd = $mrtrix/dwigradcheck
  set cmd = ($cmd $dwidir/$dwiname.nii.gz)
  set cmd = ($cmd -mask $dwidir/lowb_mask.nii.gz)
  set cmd = ($cmd -fslgrad $dwidir/${dwiname}_t.bvecs $dwidir/${dwiname}_t.bvals)
  set cmd = ($cmd -export_grad_fsl $dwidir/${dwiname}_c.bvecs $dwidir/${dwiname}_c.bvals)
  set cmd = ($cmd -force)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  # Copying corrected bvecs bvals
  mv $dwidir/${dwiname}_c.bvecs $dwidir/${dwiname}_t.bvecs
  mv $dwidir/${dwiname}_c.bvals $dwidir/${dwiname}_t.bvals 
  # Transposing gradients back
  set cmd = /autofs/space/hemera_001/users/cmaffei/ex_vivo/proc.5.2020/code/transpose_gradients.py
  set cmd = (python $cmd)
  set cmd = ($cmd $dwidir/${dwiname}_t.bvals)
  set cmd = ($cmd $dwidir/$dwiname.bvals)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
  set cmd = /autofs/space/hemera_001/users/cmaffei/ex_vivo/proc.5.2020/code/transpose_gradients.py
  set cmd = (python $cmd)
  set cmd = ($cmd $dwidir/${dwiname}_t.bvecs)
  set cmd = ($cmd $dwidir/$dwiname.bvecs)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
endif
endif 

if ($doall) then
if ($doeddy) then
  echo "Running eddy..."
  #
  # Compensate for eddy-current distortions
  #
  if (-e $dwidir/dwi_drift.nii.gz) then
     set dwiname = dwi_drift
     set bname = dwi_las
  else if (-e $dwidir/dwi_degibb.nii.gz) then
     set dwiname = dwi_degibb
     set bname = dwi_las
  else if (-e $dwidir/dwi_denoised.nii.gz) then
     set dwiname = dwi_denoised
     set bname = dwi_las
  else if (-e $dwidir/dwi_las.nii.gz) then
     set dwiname = dwi_las
     set bname = dwi_las
  else
     set dwiname = dwi_orig
     set bname = dwi_orig
  endif
  set cmd = eddy
  set cmd = ($cmd --imain=$dwidir/$dwiname.nii.gz)
  if ($?uselowbmask) then
    set cmd = ($cmd --mask=$dwidir/lowb_mask.nii.gz)
  else if ($?usehighbmask) then
    set cmd = ($cmd --mask=$dwidir/highb_mask.nii.gz)
  endif
  set cmd = ($cmd --acqp=$dwidir/acqparams.txt)
  set cmd = ($cmd --index=$dwidir/index.txt)
  set cmd = ($cmd --bvecs=$dwidir/$bname.bvecs)
  set cmd = ($cmd --bvals=$dwidir/$bname.bvals)
  #set cmd = ($cmd --estimate_move_by_susceptibility) #added v 6.0.0
  set cmd = ($cmd --out=$dwidir/dwi)
  #set cmd = ($cmd --dont_peas) #Do NOT perform post-eddy alignment of shells
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set cmd = (cp $dwidir/dwi_las.bvals $dwidir/dwi.bvals)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
   
  set cmd = xfmrot
  set cmd = ($cmd $dwidir/dwi.eddy_parameters)
  set cmd = ($cmd $dwidir/$bname.bvecs)
  set cmd = ($cmd $dwidir/dwi.bvecs)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
endif
endif 

if ($doall) then
if ($dobiascor) then
  echo "Bias correction..."
  #
  # Apply bias correction ANTS
  #
  if (-e $dwidir/dwi.nii.gz) then
     set dwiname = dwi
     set bname = dwi
  else if (-e $dwidir/dwi_degibb.nii.gz) then
     set dwiname = dwi_degibb
     set bname = dwi_las
  else if (-e $dwidir/dwi_drift.nii.gz) then
     set dwiname = dwi_drift
     set bname = dwi_las
  else if (-e $dwidir/dwi_denoised.nii.gz) then
     set dwiname = dwi_denoised
     set bname = dwi_las
  else if (-e $dwidir/dwi_las.nii.gz) then
     set dwiname = dwi_las
     set bname = dwi_las
  else
     set dwiname = dwi_orig.nii.gz
     set bname = dwi_orig
  endif

  set cmd = $mrtrix/dwibiascorrect
  set cmd = ($cmd $dwidir/$dwiname.nii.gz)
  set cmd = ($cmd -bias $dwidir/bias_est.nii.gz)
  set cmd = ($cmd -ants)
  set cmd = ($cmd -fslgrad $dwidir/$bname.bvecs $dwidir/$bname.bvals)
  set cmd = ($cmd -ants.b $antsb)
  set cmd = ($cmd $dwidir/dwi_biascorr.nii.gz)
  if ($?uselowbmask) then
    set cmd = ($cmd -mask=$dwidir/lowb_mask.nii.gz)
  else if ($?usehighbmask) then
    set cmd = ($cmd -mask=$dwidir/highb_mask.nii.gz)
  endif
  set cmd = ($cmd -force)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF 
endif
endif 

if ($doall) then
if ($dotensor) then
  echo "Tensor fitting..."
  #
  # Fit tensors
  # We only want to apply the tensor on the lower bval volumes
  #
  if (-e $dwidir/dwi_biascorr.nii.gz) then
     set dwiname = dwi_biascorr
     set bname = dwi
     set maskname = dwi_las_mask
  else if (-e $dwidir/dwi.nii.gz) then
     set dwiname = dwi
     set bname = dwi
     set maskname = dwi_las_mask
  else if (-e $dwidir/dwi_degibb.nii.gz) then
     set dwiname = dwi_degibb
     set bname = dwi_las
     set maskname = dwi_las_mask
  else if (-e $dwidir/dwi_drift.nii.gz) then
     set dwiname = dwi_drift
     set bname = dwi_las
     set maskname = dwi_las_mask
  else if (-e $dwidir/dwi_denoised.nii.gz) then
     set dwiname = dwi_denoised
     set bname = dwi_las
     set maskname = dwi_las_mask
  else if (-e $dwidir/dwi_las) then
     set dwiname = dwi_las
     set bname = dwi_las
     set maskname = dwi_las_mask
  else
     set dwiname = dwi_orig
     set bname = dwi_orig
     set maskname = dwi_orig_mask
  endif

  if ($dotensorlow) then
    #extract lowbval volumes
    set bvals = `cat $dwidir/$bname.bvals`
    set bvecs = `cat $dwidir/$bname.bvecs`
    printf '' > $dwidir/${bname}_tensorbvecs
    printf '' > $dwidir/${bname}_tensorbvals
    @ k = 1
    while ($k <= $#bvals)
    set b = `echo "$bvals[$k]/1" | bc`
    if ($b <= 8000) then
    echo $bvals[$k] >> $dwidir/${bname}_tensorbvals
    awk 'NR=='$k' {printf "%.10f %.10f %.10f\n", $1, $2, $3} ' $dwidir/$bname.bvecs >> $dwidir/${bname}_tensorbvecs
    endif
    @ k = $k + 1
    end
  
    set tenbvals = `cat $dwidir/${bname}_tensorbvals`
    set cmd = fslroi
    set cmd = ($cmd $dwidir/$dwiname.nii.gz)
    set cmd = ($cmd $dwidir/${dwiname}_tensor.nii.gz)
    set cmd = ($cmd 0 $#tenbvals)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF

    set cmd = dtifit
    set cmd = ($cmd -k $dwidir/${dwiname}_tensor.nii.gz)
    set cmd = ($cmd -o $dwidir/dtifit_lowb)
    if ($?uselowbmask) then
      set cmd = ($cmd -m $dwidir/lowb_mask.nii.gz)
    else if ($?usehighbmask) then
      set cmd = ($cmd -m $dwidir/highb_mask.nii.gz)
    endif
    set cmd = ($cmd -r $dwidir/${bname}_tensorbvecs)
    set cmd = ($cmd -b $dwidir/${bname}_tensorbvals)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF
    set cmd = (rm $dwidir/${dwiname}_tensor.nii.gz)
    $cmd |& tee -a $LF
    set cmd = (rm $dwidir/${bname}_tensorbvecs)
    $cmd |& tee -a $LF
    set cmd = ($dwidir/${bname}_tensorbvals)
    $cmd |& tee -a $LF
  endif
  set cmd = dtifit
  set cmd = ($cmd -k $dwidir/$dwiname.nii.gz)
  set cmd = ($cmd -o $dwidir/dtifit)
  set cmd = ($cmd -m $dwidir/${maskname}.nii.gz)
  set cmd = ($cmd -r $dwidir/$bname.bvecs)
  set cmd = ($cmd -b $dwidir/$bname.bvals)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

#create rough wm mask
  echo " Creating wm mask.."
  echo " -------- Creating WM mask ---------" >> $LF 
  set cmd = fslmaths
  set cmd = ($cmd $dwidir/dtifit_FA.nii.gz)
  set cmd = ($cmd -thr 0.1)
  set cmd = ($cmd -uthr 0.7)
  set cmd = ($cmd -bin)
  set cmd = ($cmd $dwidir/wm_bin.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
  set cmd = $mrtrix/maskfilter
  set cmd = ($cmd $dwidir/wm_bin.nii.gz)
  set cmd = ($cmd erode)
  set cmd = ($cmd $dwidir/wm_bin.nii.gz)
  set cmd = ($cmd -force)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
  set cmd = maskfilter 
  set cmd = ($cmd $dwidir/wm_bin.nii.gz)
  set cmd = ($cmd clean)
  set cmd = ($cmd $dwidir/wm_bin.nii.gz)
  set cmd = ($cmd -force)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
  set cmd = maskfilter  
  set cmd = ($cmd $dwidir/wm_bin.nii.gz)
  set cmd = ($cmd dilate)
  set cmd = ($cmd $dwidir/wm_bin.nii.gz)
  set cmd = ($cmd -npass 2)
  set cmd = ($cmd -force)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  # Compute SNR --> move to wm mask
  set cmd = $mrtrix/mrstats
  set cmd = ($cmd -output mean)
  set cmd = ($cmd -mask $dwidir/wm_bin.nii.gz)
  set cmd = ($cmd $dwidir/noise_est.nii.gz)
  echo "Noise" >> snr.txt
  $cmd >> snr.txt
  set cmd = fslroi
  set cmd = ($cmd $dwidir/$dwiname.nii.gz)
  set cmd = ($cmd $dwidir/${dwiname}_lowb.nii.gz)
  set cmd = ($cmd 0 1)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
  set cmd = $mrtrix/mrstats
  set cmd = ($cmd -output mean)
  set cmd = ($cmd -mask $dwidir/wm_bin.nii.gz)
  set cmd = ($cmd $dwidir/${dwiname}_lowb.nii.gz)
  echo "Signal" >> snr.txt
  $cmd >> snr.txt
endif
endif

if ($dogqi) then
  echo "Fitting GQI model..."
  # More options available. $myscripts/gqi.py -h for info
  if (-e $dwidir/dwi_biascorr.nii.gz) then
   set dwiname = dwi_biascorr
   set bname = dwi
  else if (-e $dwidir/dwi.nii.gz) then
     set dwiname = dwi
     set bname = dwi
  else if (-e $dwidir/dwi_degibb.nii.gz) then
     set dwiname = dwi_degibb
     set bname = dwi_las
  else if (-e $dwidir/dwi_drift.nii.gz) then
     set dwiname = dwi_drift
     set bname = dwi_las
  else if (-e $dwidir/dwi_denoised.nii.gz) then
     set dwiname = dwi_denoised
     set bname = dwi_las
  else if (-e $dwidir/dwi_las) then
     set dwiname = dwi_las
     set bname = dwi_las
  else
     set dwiname = dwi_orig
     set bname = dwi_orig
  endif
  #######
  set cmd = $myscriptsdir/gqi.py 
  set cmd = ($cmd $dwiname.nii.gz)
  set cmd = ($cmd $dwidir/$bname.bvecs)
  set cmd = ($cmd $dwidir/$bname.bvals)
  set cmd = ($cmd $dwidir/$dwiname.gqi.nii.gz)
  set cmd = ($cmd $difflen)
  if ($?uselowbmask) then
    set cmd = ($cmd --mask $dwidir/lowb_mask.nii.gz)
  else if ($?usehighbmask) then
    set cmd = ($cmd --mask $dwidir/highb_mask.nii.gz)
  endif
  if ($?savesh) then 
   set cmd = ($cmd --save_sh $dwidir/$dwiname.gqi_sh.nii.gz)
  endif
  if ($?outpeaks) then 
   set cmd = ($cmd --output_peaks $dwidir/$dwiname.gqi.peaks.nii.gz)
   set cmd = ($cmd --output_peaks_values $dwidir/$dwiname.gqi.peaksval.nii.gz)
   set cmd = ($cmd --output_peaks_dirs $dwidir/$dwiname.gqi.peaksdir.nii.gz)
   set cmd = ($cmd --output_peaks_indices $dwidir/$dwiname.gqi.peaksind.nii.gz)
  endif
  set cmd = ($cmd --npeaks $npeaks)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  if ($dotrack) then
    if ( -e $dwidir/$dwiname.gqi_sh.nii.gz) then 
      set fodname = $dwiname.gqi_sh.nii.gz
    else 
      set fodname = $dwiname.gqi.nii.gz
    endif
    set cmd = $myscriptsdir/tractography.py
    set cmd = (python $cmd)
    set cmd = ($cmd $dwidir/$fodname)
    if ($?uselowbmask) then
      set cmd = ($cmd $dwidir/lowb_mask.nii.gz)
    else if ($?usehighbmask) then
      set cmd = ($cmd $dwidir/highb_mask.nii.gz)
    endif
    set cmd = ($cmd $dwidir/wm_bin.nii.gz)
    set cmd = ($cmd $dwidir/gqiprob.trk)
    set cmd = ($cmd --fa $dwidir/dtifit_FA.nii.gz)
    set cmd = ($cmd --fa_thr $fa_thr)
    set cmd = ($cmd --stopping_criteria threshold)
    set cmd = ($cmd --track_algo prob)
    set cmd = ($cmd --seeds_n $seed_number)
    if ($?step_size) then
      set cmd = ($cmd --step_size $step_size)
    endif
    set cmd = ($cmd --pmf_treshold $pmf_thr)
    set cmd = ($cmd --max_angle $max_angle)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF
  endif 
endif

if ($dodsi) then 
  echo "Ftting DSI model..."
  # More options available. $myscripts/dsi.py -h for info
  if (-e $dwidir/dwi_biascorr.nii.gz) then
    set dwiname = dwi_biascorr
    set bname = dwi
  else if (-e $dwidir/dwi.nii.gz) then
    set dwiname = dwi
    set bname = dwi
  else if (-e $dwidir/dwi_degibb.nii.gz) then
    set dwiname = dwi_degibb
    set bname = dwi_las
  else if (-e $dwidir/dwi_drift.nii.gz) then
    set dwiname = dwi_drift
    set bname = dwi_las
  else if (-e $dwidir/dwi_denoised.nii.gz) then
    set dwiname = dwi_denoised
    set bname = dwi_las
  else 
    set dwiname = dwi_las
    set bname = dwi_las
  endif
  ########
  set cmd = $myscriptsdir/dsi.py 
  set cmd = ($cmd $dwiname.nii.gz)
  set cmd = ($cmd $dwidir/$bname.bvecs)
  set cmd = ($cmd $dwidir/$bname.bvals)
  set cmd = ($cmd $dwidir/$dwiname.dsi.nii.gz)
  set cmd = ($cmd $filter_width)
  if ($?uselowbmask) then
    set cmd = ($cmd --mask $dwidir/lowb_mask.nii.gz)
  else if ($?usehighbmask) then
    set cmd = ($cmd --mask $dwidir/highb_mask.nii.gz)
  endif
  if ($?savesh) then 
   set cmd = ($cmd --save_sh $dwidir/$dwiname.dsi_sh.nii.gz)
  endif
  if ($?outpeaks) then 
   set cmd = ($cmd --output_peaks $dwidir/$dwiname.dsi.peaks.nii.gz)
   set cmd = ($cmd --output_peaks_values $dwidir/$dwiname.dsi.peaksval.nii.gz)
   set cmd = ($cmd --output_peaks_dirs $dwidir/$dwiname.dsi.peaksdir.nii.gz)
   set cmd = ($cmd --output_peaks_indices $dwidir/$dwiname.dsi.peaksind.nii.gz)
  endif
  set cmd = ($cmd --r_start $r_start)
  set cmd = ($cmd --r_end $r_end)
  set cmd = ($cmd --r_step $r_step)
  set cmd = ($cmd --npeaks $npeaks)
  set cmd = ($cmd --qgrid_size $qgrid_size)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
  
  if ($dotrack) then
    if ( -e $dwidir/$dwiname.dsi_sh.nii.gz) then 
      set fodname = $dwiname.dsi_sh.nii.gz
    else 
      set fodname = $dwiname.dsi.nii.gz
    endif
    set cmd = $myscriptsdir/tractography.py
    set cmd = ($cmd $dwidir/$fodname)
    if ($?uselowbmask) then
      set cmd = ($cmd --mask $dwidir/lowb_mask.nii.gz)
    else if ($?usehighbmask) then
      set cmd = ($cmd --mask $dwidir/highb_mask.nii.gz)
    endif
    set cmd = ($cmd $dwidir/wm_bin.nii.gz)
    set cmd = ($cmd $dwidir/dsiprob.trk)
    set cmd = ($cmd --fa $dwidir/dtifit_FA.nii.gz)
    set cmd = ($cmd --fa_thr $fa_thr)
    set cmd = ($cmd --stopping_criteria threshold)
    set cmd = ($cmd --track_algo prob)
    set cmd = ($cmd --seeds_n $seed_number)
    if ($?step_size) then
      set cmd = ($cmd --step_size $step_size)
    endif
    set cmd = ($cmd --pmf_treshold $pmf_thr)
    set cmd = ($cmd --max_angle $max_angle)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF
    endif
endif


if ($doresamp) then
  if (! $?bshell) then
    echo "ERROR: Must specify shelled bvals in configuration file"
  endif
  if (! $?gshell) then
    echo "ERROR: Must specify shelled bvecs in configuration file"
  endif
  set resampdir = $dwidir/shelled_data
  if (! -d $resampdir) then
    mkdir $resampdir
  endif
  if (-e $dwidir/dwi_las.nii.gz) then
   set dwiname = dwi_las
  else
   set dwiname = dwi_orig
  endif
  set cmd = "addpath /autofs/space/hemera_001/users/cmaffei/ex_vivo/ismrm_abstract/code"
  set cmd = "$cmd; addpath /autofs/space/hemera_001/users/cmaffei/ex_vivo/ismrm_abstract/code/irt/utilities"
  set cmd = "$cmd; addpath /autofs/space/hemera_001/users/cmaffei/ex_vivo/ismrm_abstract/code/irt/nufft"
  set cmd = "$cmd; addpath /autofs/space/hemera_001/users/cmaffei/ex_vivo/ismrm_abstract/code/irt/systems" 
  set cmd = "$cmd; qresample("
  set cmd = "$cmd '$resampdir/dwi_res.nii.gz',"
  set cmd = "$cmd '$dwidir/${dwiname}.nii.gz',"
  set cmd = "$cmd '$dwidir/${dwiname}.bvals',"
  set cmd = "$cmd '$dwidir/${dwiname}.bvecs',"
  set cmd = "$cmd '$bshell',"
  set cmd = "$cmd '$gshell',"
  set cmd = "$cmd '$dwidir/${dwiname}_mask.nii.gz');"
  #echo $cmd | tee -a $LF
  #echo $cmd | matlab -nosplash  
  
  if ($dodenoise) then   
    # Apply denoising Veraart 2016, MRtrix Implementation
    echo "Denoising data..."
    # 
    set cmd = $mrtrix/dwidenoise
    set cmd = ($cmd -extent $extent)
    set cmd = ($cmd -noise $resampdir/noise_est.nii.gz)
    set cmd = ($cmd $resampdir/dwi_res.nii.gz)
    set cmd = ($cmd $resampdir/dwi_denoised.nii.gz)
    set cmd = ($cmd -force)
    #echo $cmd | tee -a $LF
    #$cmd |& tee -a $LF

    # Compute residuals
    set cmd = $mrtrix/mrcalc
    set cmd = ($cmd $resampdir/dwi_res.nii.gz)
    set cmd = ($cmd $resampdir/noise_est.nii.gz)
    set cmd = ($cmd -subtract)
    set cmd = ($cmd $resampdir/residuals.nii.gz)
    set cmd = ($cmd -quiet)
    set cmd = ($cmd -force)
    echo $cmd | tee -a $LF
    $cmd |& tee -a $LF
    # Compute RMS
    set cmd = $mrtrix/mrmath
    set cmd = ($cmd $resampdir/residuals.nii.gz)
    set cmd = ($cmd rms)
    set cmd = ($cmd -axis 3)
    set cmd = ($cmd $resampdir/rms_residuals.nii.gz)
    set cmd = ($cmd -quiet)
    set cmd = ($cmd -force)
    echo $cmd | tee -a $LF
    $cmd |& tee -a $LF
    # Compute SNR
    set cmd = $mrtrix/mrstats
    set cmd = ($cmd -output mean)
    set cmd = ($cmd -mask $dwidir/wm_bin.nii.gz)
    set cmd = ($cmd $resampdir/noise_est.nii.gz)
    echo "Noise" >> $resampdir/snr.txt
    $cmd >> $resampdir/snr.txt
    set cmd = fslroi
    set cmd = ($cmd $resampdir/dwi_res.nii.gz)
    set cmd = ($cmd $resampdir/dwi_res_lowb.nii.gz)
    set cmd = ($cmd 0 1)
    echo $cmd | tee -a $LF
    $cmd |& tee -a $LF
    set cmd = $mrtrix/mrstats
    set cmd = ($cmd -output mean)
    set cmd = ($cmd -mask $dwidir/wm_bin.nii.gz)
    set cmd = ($cmd $resampdir/dwi_res_lowb.nii.gz)
    echo "Signal" >> $resampdir/snr.txt
    $cmd >> $resampdir/snr.txt
  endif

  if ($dodegibb) then
    if (-e $resampdir/dwi_denoised.nii.gz) then
      set dwiname = dwi_denoised.nii.gz
    else
       set dwiname = dwi_res.nii.gz
    endif
    #Specify acquisition plane in -axes option:
    # 0,1 axial 0,2 coronal 1,2 sagittal
    if ($?acq_plane) then
    if ($acq_plane == axial ) then 
       set axes = 0,1
    else if ($acq_plane == coronal) then 
       set axes = 0,2
    else 
       set exes = 1,2
    endif
    else 
    echo "Acq plane non specified. Coronal will be used."
    set axes = 0,2
    endif
    #
    # Attempts to remove Gibbs ringing artefacts from MR
    # Kellner et al 2016, MRtrix implementation
    # 
    set cmd = $mrtrix/mrdegibbs
    set cmd = ($cmd -axes $axes)
    set cmd = ($cmd $resampdir/$dwiname)
    set cmd = ($cmd $resampdir/dwi_degibb.nii.gz) 
    set cmd = ($cmd -force) 
    echo $cmd | tee -a $LF
    $cmd |& tee -a $LF
  endif

  if ($doscale) then 
    if (-e $resampdir/dwi_degibb.nii.gz) then 
      set dwiname = dwi_degibb.nii.gz 
    else if (-e $resampdir/dwi_denoised.nii.gz) then
      set dwiname = dwi_denoised.nii.gz
    else 
      set dwiname = dwi_res.nii.gz
    endif
    set min = `fslstats $resampdir/$dwiname -R | awk '{ print $1 }'`
    set max = `fslstats $resampdir/$dwiname -R | awk '{ print $2 }'`
    set scaling = `echo "scale=5; 32767.0 / ( $max - $min )" | bc`
    fslmaths $resampdir/$dwiname -sub $min -mul $scaling $resampdir/$dwiname
  endif

  if ($dodrift) then
    if (-e $resampdir/dwi_degibb.nii.gz) then
       set dwiname = dwi_degibb.nii.gz
    else if  (-e $resampdir/dwi_denoised.nii.gz) then
       set dwiname = dwi_denoised.nii.gz
    else
       set dwiname = dwi_res.nii.gz
    endif
    #
    # Normalize images to compensate for temperature drift
    # 
    set lowb = 1
    set cmd = "addpath $FREESURFER_HOME/matlab"
    set cmd = "addpath /autofs/space/turan_001/users/lzollei/dev/matlab"
    set cmd = "$cmd; fix_exvivo_dwi_drift("
    set cmd = "$cmd '$resampdir/dwi_drift.nii.gz', "
    set cmd = "$cmd '$resampdir/drift', "
    set cmd = "$cmd '$resampdir/$dwiname', "
    set cmd = "$cmd $#lowb);"
    echo $cmd | tee -a $LF
    echo $cmd | matlab -nosplash
  endif
  #
    
  if ($gradcheck) then
    #checking gradients
    echo "Checking gradients..."
    echo "-------------------------------" >> $LF
    echo "------Checking Gradients-------" >> $LF
   
    # Transpose bvec/bval files to make mrtrix happy 
    set cmd = /autofs/space/hemera_001/users/cmaffei/ex_vivo/proc.5.2020/code/transpose_gradients.py
    set cmd = (python $cmd)
    set cmd = ($cmd $gshell)
    set cmd = ($cmd $resampdir/bvecs)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF

    set cmd = /autofs/space/hemera_001/users/cmaffei/ex_vivo/proc.5.2020/code/transpose_gradients.py
    set cmd = (python $cmd)
    set cmd = ($cmd $bshell)
    set cmd = ($cmd $resampdir/bvals)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF
    
    set cmd = $mrtrix/dwigradcheck
    set cmd = ($cmd $resampdir/dwi_res.nii.gz)
    set cmd = ($cmd -mask $dwidir/lowb_mask.nii.gz)
    set cmd = ($cmd -fslgrad $resampdir/bvecs $resampdir/bvals)
    set cmd = ($cmd -export_grad_fsl $resampdir/dwi.bvecs $resampdir/dwi.bvals)
    set cmd = ($cmd -force)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF
  endif

  if ($doeddy) then
    echo "Running eddy..."
    #
    # Compensate for eddy-current distortions
    #
    if (-e $resampdir/dwi_drift.nii.gz) then
       set dwiname = dwi_drift
    else if (-e $dwidir/dwi_degibb.nii.gz) then
       set dwiname = dwi_degibb
    else if (-e $dwidir/dwi_denoised.nii.gz) then
       set dwiname = dwi_denoised
    else
       set dwiname = dwi_res
    endif

    awk 'NR==1{print }' $dwidir/acqparams.txt > $resampdir/acqparams.txt
    set vols = `fslnvols $resampdir/dwi_res.nii.gz`
    printf '' > $resampdir/index.txt
    @ k = 1
    while ($k <= $vols)
    echo 1 >> $resampdir/index.txt
    @ k = $k + 1
    end
    set cmd = eddy
    set cmd = ($cmd --imain=$resampdir/$dwiname.nii.gz)
    set cmd = ($cmd --mask=$dwidir/lowb_mask.nii.gz)
    set cmd = ($cmd --acqp=$resampdir/acqparams.txt)
    set cmd = ($cmd --index=$resampdir/index.txt)
    set cmd = ($cmd --bvecs=$resampdir/bvecs)
    set cmd = ($cmd --bvals=$resampdir/bvals)
    #set cmd = ($cmd --estimate_move_by_susceptibility) #added v 6.0.0
    set cmd = ($cmd --out=$resampdir/dwi)
    #set cmd = ($cmd --dont_peas) #Do NOT perform post-eddy alignment of shells
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF

    set cmd = xfmrot
    set cmd = ($cmd $resampdir/dwi.eddy_parameters)
    set cmd = ($cmd $resampdir/dwi.bvecs)
    set cmd = ($cmd $resampdir/dwi.bvecs)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF
  endif

  if ($dobiascor) then
    echo "Bias correction..."
    #
    # Apply bias correction ANTS
    #
    if (-e $resampdir/dwi.nii.gz) then
       set dwiname = dwi
    else if (-e $resampdir/dwi_drift.nii.gz) then
       set dwiname = dwi_drift
    else if (-e $resampdir/dwi_degibb.nii.gz) then
       set dwiname = dwi_degibb
    else if (-e $resampdir/dwi_denoised.nii.gz) then
       set dwiname = dwi_denoised
    else
       set dwiname = dwi_res
    endif

    set cmd = $mrtrix/dwibiascorrect
    set cmd = ($cmd $resampdir/$dwiname.nii.gz)
    set cmd = ($cmd -bias $resampdir/bias_est.nii.gz)
    set cmd = ($cmd -ants)
    set cmd = ($cmd -fslgrad $resampdir/dwi.bvecs $resampdir/dwi.bvals)
    set cmd = ($cmd -ants.b $antsb)
    set cmd = ($cmd $resampdir/dwi_biascorr.nii.gz)
    set cmd = ($cmd -mask $dwidir/lowb_mask.nii.gz)
    set cmd = ($cmd -force)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF 
  endif

  if ($dotensor) then
    echo "Tensor fitting..."
    #
    # Fit tensors
    # We only want to apply the tensor on the lower bval volumes
    #
    if (-e $resampdir/dwi_biascorr.nii.gz) then
       set dwiname = dwi_biascorr
    else if (-e $dwidir/dwi.nii.gz) then
       set dwiname = dwi
    else if (-e $dwidir/dwi_drift.nii.gz) then
       set dwiname = dwi_drift
    else if (-e $dwidir/dwi_degibb.nii.gz) then
       set dwiname = dwi_degibb
    else if (-e $dwidir/dwi_denoised.nii.gz) then
       set dwiname = dwi_denoised
    else
       set dwiname = dwi_res
    endif 

    #extract lowbval volumes
    set cmd = $mrtrix/dwiextract
    set cmd = ($cmd $resampdir/$dwiname.nii.gz)
    set cmd = ($cmd $resampdir/${dwiname}_tensor.nii.gz)
    set cmd = ($cmd -shells 0,4000,8000)
    set cmd = ($cmd -fslgrad $resampdir/dwi.bvecs $resampdir/dwi.bvals)
    set cmd = ($cmd -export_grad_fsl $resampdir/bvecs_tensor $resampdir/bvals_tensor)
    set cmd = ($cmd -force)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF

    set cmd = dtifit
    set cmd = ($cmd -k $resampdir/${dwiname}_tensor.nii.gz)
    set cmd = ($cmd -o $resampdir/dtifit)
    set cmd = ($cmd -m $dwidir/lowb_mask.nii.gz)
    set cmd = ($cmd -r $resampdir/bvecs_tensor)
    set cmd = ($cmd -b $resampdir/bvals_tensor)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF
endif
       
  if ($domrtrix) then
    set mrtrixdir = $resampdir/mrtrix
    if (! -d $mrtrixdir) then
      mkdir $resampdir/mrtrix
    endif 
    if (-e $resampdir/dwi_biascorr.nii.gz) then
       set dwiname = dwi_biascorr
    else if (-e $dwidir/dwi.nii.gz) then
       set dwiname = dwi
    else if (-e $dwidir/dwi_drift.nii.gz) then
       set dwiname = dwi_drift
    else if (-e $dwidir/dwi_degibb.nii.gz) then
       set dwiname = dwi_degibb
    else if (-e $dwidir/dwi_denoised.nii.gz) then
       set dwiname = dwi_denoised
    else
       set dwiname = dwi_res
    endif
	
    #convert nifit to mif mrtrix file format 
    set cmd = $mrtrix/mrconvert
    set cmd = ($cmd $resampdir/$dwiname.nii.gz)
    set cmd = ($cmd -fslgrad $resampdir/dwi.bvecs $resampdir/dwi.bvals)
    set cmd = ($cmd $mrtrixdir/$dwiname.mif)
    set cmd = ($cmd -force)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF
   
    echo " Running deterministic tractography..." 
    # whole brain tensor deterministic tractography 
    set cmd = $mrtrix/tckgen
    set cmd = ($cmd -algorithm Tensor_Det)
    set cmd = ($cmd $mrtrixdir/$dwiname.mif)
    set cmd = ($cmd $mrtrixdir/$dwiname.tensordet.tck)
    set cmd = ($cmd -select $str_number)
    set cmd = ($cmd -seed_image $dwidir/wm_bin.nii.gz)
    set cmd = ($cmd -mask $dwidir/lowb_mask.nii.gz)
    set cmd = ($cmd -force)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF
    set cmd = $myscriptsdir/tck2trk.py
    set cmd = ($cmd $mrtrixdir/$dwiname.tensordet.tck)
    set cmd = ($cmd $mrtrixdir/$dwiname.tensordet.trk)
    set cmd = ($cmd $dwidir/lowb_mask.nii.gz)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF

    # Compute response function for CSD (Tournier et al 2007)
    set cmd = $mrtrix/maskfilter
    set cmd = ($cmd $dwidir/wm_bin.nii.gz)
    set cmd = ($cmd erode)
    set cmd = ($cmd $mrtrixdir/wm_bin_er.nii.gz)
    set cmd = ($cmd -force)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF
    
    # Automatic tournier
    set cmd = $mrtrix/dwi2response
    set cmd = ($cmd tournier)
    set cmd = ($cmd -mask $mrtrixdir/wm_bin_er.nii.gz)
    set cmd = ($cmd -voxels $mrtrixdir/sfvox_tournier.nii.gz)
    set cmd = ($cmd $mrtrixdir/$dwiname.mif)
    set cmd = ($cmd $mrtrixdir/rf_tournier.txt)
    set cmd = ($cmd -force)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF
    # Compute fODFs
    set cmd = $mrtrix/dwi2fod
    set cmd = ($cmd csd)
    set cmd = ($cmd $mrtrixdir/$dwiname.mif)
    set cmd = ($cmd $mrtrixdir/rf_tournier.txt)
    set cmd = ($cmd $mrtrixdir/csd_fods_tournier.nii.gz)
    set cmd = ($cmd -mask $dwidir/lowb_mask.nii.gz) 
    set cmd = ($cmd -force)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF
    # Computing CSD Probabilistic Tractography 
    set cmd = $mrtrix/tckgen
    set cmd = ($cmd $mrtrixdir/csd_fods_tournier.nii.gz)
    set cmd = ($cmd $mrtrixdir/csdprob_tournier.tck)
    set cmd = ($cmd -seed_random_per_voxel $dwidir/wm_bin.nii.gz $seed_number)
    set cmd = ($cmd -mask $dwidir/lowb_mask.nii.gz)
    set cmd = ($cmd -angle $max_angle)
    set cmd = ($cmd -cutoff $cutoff)
    set cmd = ($cmd -force)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF
    # Coverting tck output to trk
    set cmd = $myscriptsdir/tck2trk.py
    set cmd = ($cmd $mrtrixdir/csdprob_tournier.tck)
    set cmd = ($cmd $mrtrixdir/csdprob_tournier.trk)
    set cmd = ($cmd $dwidir/lowb_mask.nii.gz)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF
    #
    
    # Manual 
    if ($?wm_vox && $?gm_vox) then 
      set cmd = $mrtrix/dwi2response
      set cmd = ($cmd manual)
      set cmd = ($cmd $mrtrixdir/$dwiname.mif)
      set cmd = ($cmd $dwidir/$wm_vox)
      set cmd = ($cmd $mrtrixdir/rf_man_wm.txt)
      set cmd = ($cmd -lmax 0,4,6,8)
      set cmd = ($cmd -force)
      echo $cmd |& tee -a $LF
      $cmd |& tee -a $LF
      set cmd = $mrtrix/dwi2response
      set cmd = ($cmd manual)
      set cmd = ($cmd $mrtrixdir/$dwiname.mif)
      set cmd = ($cmd $dwidir/$gm_vox)
      set cmd = ($cmd $mrtrixdir/rf_man_gm.txt)
      set cmd = ($cmd -lmax 0,0,0,0)
      set cmd = ($cmd -force)
      echo $cmd |& tee -a $LF
      $cmd |& tee -a $LF
      # Computing FODs
      set cmd = $mrtrix/dwi2fod
      set cmd = ($cmd msmt_csd)
      set cmd = ($cmd $mrtrixdir/$dwiname.mif)
      set cmd = ($cmd $mrtrixdir/rf_man_gm.txt)
      set cmd = ($cmd $mrtrixdir/csd_fods_gm.nii.gz)
      set cmd = ($cmd $mrtrixdir/rf_man_wm.txt)
      set cmd = ($cmd $mrtrixdir/csd_fods_wm.nii.gz)
      set cmd = ($cmd -mask $dwidir/lowb_mask.nii.gz)
      set cmd = ($cmd -force)
      echo $cmd |& tee -a $LF
      $cmd |& tee -a $LF
      # Computing CSD Probabilistic Tractography 
      set cmd = $mrtrix/tckgen
      set cmd = ($cmd $mrtrixdir/csd_fods_wm.nii.gz)
      set cmd = ($cmd $mrtrixdir/csdprob_manual.tck)
      set cmd = ($cmd -seed_random_per_voxel $dwidir/wm_bin.nii.gz $seed_number)
      set cmd = ($cmd -mask $dwidir/lowb_mask.nii.gz)
      set cmd = ($cmd -angle $max_angle)
      set cmd = ($cmd -cutoff $cutoff)
      set cmd = ($cmd -force)
      echo $cmd |& tee -a $LF
      $cmd |& tee -a $LF
      # Coverting tck output to trk
      set cmd = $myscriptsdir/tck2trk.py
      set cmd = ($cmd $mrtrixdir/csdprob_manual.tck)
      set cmd = ($cmd $mrtrixdir/csdprob_manual.trk)
      set cmd = ($cmd $dwidir/lowb_mask.nii.gz)
      echo $cmd |& tee -a $LF
      $cmd |& tee -a $LF
    endif

    # Compute FOD base WM mask
    set cmd = fslmaths
    set cmd = ($cmd $mrtrixdir/csd_fods_wm.nii.gz)
    set cmd = ($cmd -thr 0.11 -bin)
    set cmd = ($cmd $mrtrixdir/wm_fod_mask.nii.gz)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF
    set cmd = $mrtrix/maskfilter
    set cmd = ($cmd $mrtrixdir/wm_fod_mask.nii.gz)
    set cmd = ($cmd erode)
    set cmd = ($cmd $mrtrixdir/wm_fod_mask_er.nii.gz)
    set cmd = ($cmd -force)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF
  endif
  
  if ($docsd) then 
    # CSD DIPY, Tax et al recursive algorithm
    echo "Fitting CSD DIPY..."
    echo "-------------------------------" >> $LF
    echo "--------- CSD DIPY ------------" >> $LF
    set myscriptsdir = /autofs/space/hemera_001/users/cmaffei/my_scripts
    set dipydir = $resampdir/dipy
    set mrtrixdir = $resampdir/mrtrix
    if (! -d $dipydir) then
      mkdir $resampdir/dipy
    endif
    if (-e $dwidir/dwi_biascorr.nii.gz) then
      set dwiname = dwi_biascorr
    else if (-e $dwidir/dwi.nii.gz) then
      set dwiname = dwi
    else if (-e $dwidir/dwi_degibb.nii.gz) then
       set dwiname = dwi_degibb
    else if (-e $dwidir/dwi_drift.nii.gz) then
       set dwiname = dwi_drift
    else if (-e $dwidir/dwi_denoised.nii.gz) then
       set dwiname = dwi_denoised
    else
       set dwiname = dwi_res
    endif

    set cmd = $myscriptsdir/csd.py
    set cmd = ($cmd $resampdir/$dwiname.nii.gz)
    set cmd = ($cmd $resampdir/dwi.bvecs)
    set cmd = ($cmd $resampdir/dwi.bvals)
    set cmd = ($cmd $dipydir/csd_fods.nii.gz)
    set cmd = ($cmd $mrtrixdir/wm_bin_er.nii.gz)
    set cmd = ($cmd --init_fa $init_fa)
    set cmd = ($cmd --peak_thr $peak_thr)
    set cmd = ($cmd --save_sh_mrtrix $dipydir/csd_fods_mrtrixbasis.nii.gz)
    set cmd = ($cmd --sh_order $lmax)
    if ($?outpeaks) then
      set cmd = ($cmd --output_peaks $dipydir/csd.peaks.nii.gz)
      set cmd = ($cmd --output_peaks_values $dipydir/csd.peaksval.nii.gz)
      set cmd = ($cmd --output_peaks_dirs $dipydir/csd.peaksdir.nii.gz)
      set cmd = ($cmd --output_peaks_indices $dipydir/csd.peaksind.nii.gz)
    endif
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF
  
    if ($dotrack) then
      set cmd = $myscriptsdir/tractography.py
      set cmd = ($cmd $dipydir/csd_fods.nii.gz)
      if ($?uselowbmask) then
       set cmd = ($cmd  $dwidir/lowb_mask.nii.gz)
      else if ($?usehighbmask) then
       set cmd = ($cmd $dwidir/highb_mask.nii.gz)
      endif
      set cmd = ($cmd $dwidir/wm_bin.nii.gz)
      set cmd = ($cmd $dipydir/csdprob.trk)
      set cmd = ($cmd --fa $resampdir/dtifit_FA.nii.gz)
      set cmd = ($cmd --fa_thr $fa_thr)
      set cmd = ($cmd --stopping_criteria threshold)
      set cmd = ($cmd --track_algo prob)
      set cmd = ($cmd --seeds_n $seed_number)
      if ($?step_size) then
        set cmd = ($cmd --step_size $step_size)
      endif
      set cmd = ($cmd --pmf_treshold $pmf_thr)
      set cmd = ($cmd --max_angle $max_angle)
      echo $cmd |& tee -a $LF
      $cmd |& tee -a $LF
    endif
endif

set dtkdir = /usr/pubsw/packages/dtk/0.6.4.1_patched

if ($dodtk) then
  mkdir $dwidir/dtk
  set dtk = $dwidir/dtk
  #
  # Fit ODFs
  #
  if (-e $dwidir/dwi.nii.gz) then
    set dwiname = dwi
  else
    set dwiname = dwi_las
  endif
  
  # Move low-b images to the beginning of the series (for dtk)
  set cmd = mri_convert
  set cmd = ($cmd $dwidir/$dwiname.nii.gz)
  set cmd = ($cmd -f $lowblist $highblist)
  set cmd = ($cmd $dtk/tmp.$dwiname.nii)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  cp /dev/null $dtk/gradients.txt


  foreach k ($highblist)
    set k1 = `echo "$k+1" | bc`
    set xyz = `sed -n "$k1 p" $dwidir/$dwiname.bvecs`
    echo "$xyz[1], $xyz[2], $xyz[3]" >> $dtk/gradients.txt
  end

  set cmd = $dtkdir/hardi_mat
  set cmd = ($cmd $dtk/gradients.txt)
  set cmd = ($cmd $dtk/qbi_mat.dat)
  set cmd = ($cmd -ref $dtk/tmp.$dwiname.nii)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set nlow = $#lowblist
  set ndir = `echo "$#highblist + 1" | bc`

  setenv DSI_PATH $dtkdir/matrices

  set cmd = $dtkdir/odf_recon
  set cmd = ($cmd $dtk/tmp.$dwiname.nii)
  set cmd = ($cmd $ndir 181)
  set cmd = ($cmd $dtk/qbi)
  set cmd = ($cmd -b0 $nlow)
  set cmd = ($cmd -mat $dtk/qbi_mat.dat)
  set cmd = ($cmd -nt -p 3 -sn 1 -ot nii)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
  
  set cmd = $dtkdir/dti_recon
  set cmd = ($cmd $dwidir/dwi.nii.gz)
  set cmd = ($cmd $dtk/dti)
  set cmd = ($cmd -gm $dwidir/$dwiname.bvecs)
  set cmd = ($cmd -b 4080)
  set cmd = ($cmd -b0 $nlow)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF 

  # Clean up temporary files
  set cmd = (rm -f $dtk/tmp.$dwiname.nii $dtk/qbi_mat.dat)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  #
  # The foreach loop below run various inversions and swaps so you can check which version is correct
  # Run q-ball tractography
  #

  setenv DSI_PATH $dtkdir/matrices

  ##  Brian said no inversion, and yz swap is correct one, so following line should be used:
  ##  set trkfile = qbi.inv,no.swap,yz.trk  
  #
  # LZ (11/27/2017): the above is outdated - it is either (swap no, inv z) that works or (swap no, inv no) with a newer / correct version of DTK -- hard to track when that version got released and what version was used...
  #

  foreach inv (no x y z)
   foreach swap (no xy yz zx)
     set trkfile = qbi.inv,$inv.swap,$swap.trk
     set cmd = $dtkdir/odf_tracker
     set cmd = ($cmd $dwidir/qbi)
     set cmd = ($cmd $dwidir/tmp.qbi.trk)
     set cmd = ($cmd -at $angthresh)
     set cmd = ($cmd -m $dwidir/qbi_dwi.nii)
     set cmd = ($cmd -it nii)
     if ($inv != no)	set cmd = ($cmd -i$inv)
     if ($swap != no)	set cmd = ($cmd -s$swap)
     echo $cmd |& tee -a $LF
     $cmd |& tee -a $LF
	    
     set cmd = ($dtkdir/spline_filter $dwidir/tmp.qbi.trk 1 $dwidir/$trkfile)
     echo $cmd |& tee -a $LF
     $cmd |& tee -a $LF
	    
    # Clean up temporary files
     set cmd = (rm -f $dwidir/tmp.qbi.trk)
     echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF
    end
  end
endif


if ($doprebed) then
  #
  # Do bedpostx pre-processing
  # This is done separately because it needs to be run on a machine
  # with more memory than launchpad nodes
  #
  if (-e $dwidir/dwi.nii.gz) then
    set dwiname = dwi
    set maskname = ${mask}
  else
    set dwiname = dwi_las
    set maskname = ${mask}_las
  endif

  set cmd = (ln -sf $dwidir/$dwiname.nii.gz $dwidir/data.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set cmd = (ln -sf $dwidir/$dwiname.bvals $dwidir/bvals)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set cmd = (ln -sf $dwidir/$dwiname.bvecs $dwidir/bvecs)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set cmd = (ln -sf $dwidir/${maskname}_brain_mask.nii.gz)
  set cmd = ($cmd $dwidir/nodif_brain_mask.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  mkdir -p $dwidir.bedpostX
  mkdir -p $dwidir.bedpostX/diff_slices
  mkdir -p $dwidir.bedpostX/logs
  mkdir -p $dwidir.bedpostX/xfms

  set cmd = (bedpostx_preproc.sh $dwidir)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
endif

if ($dobed) then
  #
  # Run bedpostx main processing (slice by slice)
  #
  set nslice = `mri_info --nslices $dwidir/lowb.nii.gz`
  set cmdfile = $dwidir.bedpostX/commands.txt
  rm -f $cmdfile

  @ k = 0
  while ($k < $nslice)
    set cmd = bedpostx_single_slice.sh
    set cmd = ($cmd $dwidir $k)
    set cmd = ($cmd --nf=$nstick --fudge=1 --bi=1000 --nj=1250)
    set cmd = ($cmd --se=25 --model=1 --cnonlinear)
    echo $cmd >> $cmdfile
    @ k = $k + 1
  end

  set cmd = fsl_sub_mgh
  set cmd = ($cmd -l $dwidir.bedpostX/logs)
  set cmd = ($cmd -m a -N bedpostx)
  set cmd = ($cmd -t $cmdfile)
  set cmd = ($cmd -a ppn=1,vmem=32gb)
  echo "Run on launchpad:" |& tee -a $LF
  echo $cmd |& tee -a $LF
endif

if ($dopostbed) then
  #
  # Do bedpostx post-processing
  # This is done separately because it needs to be run on a machine
  # with more memory than launchpad nodes
  #
  # NOTE, cleanup for prebed happens here!
  #
  set cmd = (bedpostx_postproc.sh $dwidir)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
endif

if ($doprobtrk) then
  #
  # Probabilistic seed-based targeted tractography -- all seeds are used as targets as well (LZ: 11/27/2017)
  #
  set seedlist = ($seedroidir/*.nii)

  echo "Found $#seedlist seed ROIs in $seedroidir"

  set step = `mri_info --cres $seedlist[1]`
  set step = `echo "$step / 3" | bc -l`
  set step = `printf '%g' $step`

  foreach seed ($seedlist)
    set seedname = `basename $seed | sed 's/.nii//'`

    foreach target ($seedlist)
      set targetname = `basename $target | sed 's/.nii//'`
      if !($seedname == $targetname) then
        set outdir = $dwidir/dpath.uncorrected.targeted.$seedname.2.$targetname # -pd is not on
        mkdir -p $outdir
      
        if !(-e $dwidir/rois/target$targetname.txt) then
          echo $dwidir/rois/$targetname.nii >> $dwidir/rois/target$targetname.txt
        endif 
 
        set cmd = probtrackx2
        set cmd = ($cmd -x $seed)
        set cmd = ($cmd -s $dwidir.bedpostX/merged)
        set cmd = ($cmd -m $dwidir.bedpostX/nodif_brain_mask)
        set cmd = ($cmd --dir=$outdir)
        set cmd = ($cmd -l --onewaycondition)
        set cmd = ($cmd -c 0.2 -S 2000 -P 5000) # Note: use -c .5 curvature threshold (vs .2, the default) in order to have correspondance with the 60 degrees angle threshold of deterministic tractography
        set cmd = ($cmd --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0)
        set cmd = ($cmd --steplength=$step)
        set cmd = ($cmd --forcedir --opd) 
        set cmd = ($cmd --targetmasks=$dwidir/rois/target$targetname.txt --stop=$dwidir/rois/target$targetname.txt --waypoints=$dwidir/rois/target$targetname.txt --os2t --s2tastext)
        echo $cmd |& tee -a $LF
        # $cmd |& tee -a $LF
        pbsubmit -m `whoami` -q p30 -c "$cmd" -l nodes=1:ppn=4,vmem=28gb



        set outdir = $dwidir/dpath.corrected.targeted.$seedname.2.$targetname # -pd is off while dopd=0 ("Correct path distribution for the length of the pathways")
        mkdir -p $outdir
	if ($dopd) then
  	  set cmd = ($cmd --pd)
	endif
        echo $cmd |& tee -a $LF
        # $cmd |& tee -a $LF
        pbsubmit -m `whoami` -q p30 -c "$cmd" -l nodes=1:ppn=4,vmem=28gb
      endif
    end
  end
endif
