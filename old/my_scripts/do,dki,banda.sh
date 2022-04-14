#!/bin/tcsh
cd /space/erebus/1/users/data/trc/

set_anaconda
#DKI fit (DIPY 1.1.0)
foreach sub (BANDA*)
	set data =  $sub/dmri/dwi.nii.gz
	set bvals = $sub/dmri/dwi.bvals
	set bvecs = $sub/dmri/dwi.bvecs
	set mask = $sub/dmri/nodif_brain_mask.nii.gz
	#dipy_fit_dki $data $bvals $bvecs $mask --out_ak dkifit_ak.nii.gz --out_mk dkifit_mk.nii.gz --out_rk dkifit_rk.nii.gz --out_ga dkifit_ga.nii.gz --out_dir $sub/dmri/ --save_metrics ak mk ga rk
end

#DMRI PATHSTATS
foreach sub (BANDA*)
	echo $sub
	foreach tract ($sub/dpath/*avg16_syn_bbr) #new atlas
	echo $tract
	set intrc = $sub/dpath/${tract:t}
	set path_name = `echo ${tract:t} | awk '{print substr ( $0,1,5 ) }'`
	set out = $tract/dki_pathstats.overall.txt
	set out_voxel = $tract/dki_pathstats.byvoxel.txt
	set meas = ($sub/dmri/dkifit_mk.nii.gz $sub/dmri/dkifit_rk.nii.gz $sub/dmri/dkifit_ak.nii.gz $sub/dmri/dkifit_ga.nii.gz)
	set dtbase = $sub/dmri/dtifit
	
	set cmd = /usr/local/freesurfer/7.1-tracula-beta/bin/dmri_pathstats
	set cmd = ($cmd --intrc $intrc --dtbase $dtbase)
	set cmd = ($cmd --meas $meas --measname MK GA AK RK)
	set cmd = ($cmd --path $path_name --subj $sub)
	set cmd = ($cmd --out $out --outvox $out_voxel) 
	echo $cmd
	$cmd
	end
end
