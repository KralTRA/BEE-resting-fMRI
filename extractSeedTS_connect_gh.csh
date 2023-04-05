#!/bin/tcsh

#---
# Preprocessing script for BEE rs-fMRI data - starting with seed time-series extraction
#---
#---
# Author: R.M. Birn; Adapted by T.R.A. Kral
#---

#---------------------------
# Paths:
#---------------------------
set RootPath   = /path/
set ScriptPath = $RootPath/rsfc/scripts
set RoiPath    = $RootPath/UNCInfant012Atlas_20140325

#--- file containing list of seed region ---
set SeedFile = $ScriptPath/list.ROIs.txt

#---------------------------
# Parameters:
#---------------------------
set ign        = 10 
set dxyz       = 2
set fwhm       = 6
set mot_thresh = 0.2

set template   = /UNCInfant012Atlas_20140325/infant-neo-withCerebellum.nii

set SKIP_AHEAD     = 1
set BANDPASS       = 1
set CLEANUP        = 1
set DEBUG          = 1
set ALT_SKULLSTRIP = 0
set CONN_BOTH_RUNS = 1
set CENSOR_PREV    = 1

if ($CENSOR_PREV) then
   set cens_prev_label  = "_wPrev"
   set cens_prev_string = "-censor_prev_TR"
else
   set cens_prev_label  = ""
   set cens_prev_string = ""
endif

set file_info1 = "task-rest_dir-PA_run-02_bold_Unwarp"
set file_info2 = "task-rest_dir-AP_run-02_bold_Unwarp"
      
#---------------------------
# Loop:
#---------------------------
foreach dir ($argv) #(`cat $ScriptPath/dirs_good_0.4wP_excludeCoilIssue.txt`) #(`cat $ScriptPath/dirs_1_good.txt`) #(`cat $ScriptPath/dirs_scan2.txt`) #
   set subj = $dir

   set DataPath     = $RootPath/rsfc/processed-rmb/$subj/func
   set EpiPath		= /study/sjshort_ebds/rsfc/processed-trak/withBP
   set OrigDataPath = $DataPath
   set StrucPath	= $RootPath/rsfc/processed-rmb/$subj/struc
   
   if (! -e $EpiPath) then
   		mkdir $EpiPath
   endif
   
   set LOG_FILE     = $EpiPath/log.proc.txt
   
   set AnatFile     = $subj-T2_deob_nbs_uni.nii.gz
   set AnatFile_ns  = $subj-T2_deob_nbs_uni_ns.nii.gz
   set AnatFile_alignC   = $subj-T2_deob_nbs_uni_ns_alignC.nii.gz
   set AnatFileRS   = $subj-T2_deob_nbs_uni_ns_alignC_RS.nii.gz
   set AnatFile_atlas = $subj-T2_deob_nbs_uni_ns_alignC_RS.aw.nii
   
   echo "-------------------" |&  tee -a $LOG_FILE
   echo "Processing $dir ..." |&  tee -a $LOG_FILE
   
   cd $DataPath
   
    #--- find fMRI runs ---
	set fmri_runs = `/bin/ls $DataPath | grep _bold_Unwarp.nii.gz`
	
	#echo "fmri_runs: $fmri_runs"
	if ("$fmri_runs" == "") then
  		echo "---------- Could not find fmri runs" |& tee -a $LOG_FILE
   		exit()
	endif
	
# 	set basefile = sub-$subj'_task-rest_dir-AP_run-02_bold_Unwarp.deob.nii.gz['$ign']'
	set baseprefix = sub-$subj'_'$file_info1'.deob.nii.gz'
	set basesuffix = "[10]"
	
	foreach file_orig ($fmri_runs)
   
   		echo "file_orig: $file_orig"
   		set filename = $file_orig:r:r
      
#    		set NT = `3dinfo -nt $OrigDataPath/$file_orig`
   		set NT = `3dinfo -nt $DataPath/$file_orig`
   		if ($NT < 2) then
      		echo "---------- $file_orig has less than 2 time points ($NT)" |& tee -a $LOG_FILE
      		continue
   		else
   
      	set file0     = $filename.deob.rb10.aw.NgwcmdX${mot_thresh}.bp.s${fwhm}
      
      	set MaskTotal = Mask.Brain.B.${file0}.r.nii.gz
           
      	set motfile = mot.${file0}.1D
      
      if ($CENSOR_PREV) then
         set cens_prev_label  = "_wPrev"
         set cens_prev_string = "-censor_prev_TR"
      else
         set cens_prev_label  = ""
         set cens_prev_string = ""
      endif
      
      set CensorFilePrefix = motcensor_${mot_thresh}_${file0}${cens_prev_label}
      set CensorFile       = ${CensorFilePrefix}_censor.1D
      

      #====================================================================================
      CONNECTIVITY:
      #------------------------------------------------------------------------------------
      #--- Connectivity ---
     foreach seed (`cat $SeedFile`)

# 	 if (! -e $EpiPath/ts.$seed.$file0.1D ) then
	 echo "---Extracting data from ($seed)---"  |& tee -a $LOG_FILE
	 3dmaskave -mask $RoiPath/$seed.nii.gz -quiet $EpiPath/$file0.nii.gz > $EpiPath/tmp.ts.$seed.$file0.1D
	 1dnorm -demean $EpiPath/tmp.ts.$seed.$file0.1D - > $EpiPath/ts.$seed.$file0.1D
# 	 endif
	 
	 endif
 	end 
   end
   #====================================================================================
   CONNECTIVITY_BOTHRUNS:
   #------------------------------------------------------------------------------------
   if ($CONN_BOTH_RUNS) then
   
   cd $EpiPath
   
   set cor = aw.NgwcmdX${mot_thresh}.bp.s${fwhm}
   set file_1 = sub-$subj'_'$file_info1'.deob.rb10.'$cor
   set file_2 = sub-$subj'_'$file_info2'.deob.rb10.'$cor
   set file_cat = sub-$subj'_resting_cat.'$cor
   
   set CensorFile_1   = motcensor_${mot_thresh}_sub-$subj'_'$file_info1'.deob_wPrev_censor.1D'
   set CensorFile_2   = motcensor_${mot_thresh}_sub-$subj'_'$file_info2'.deob_wPrev_censor.1D'
   set CensorFile_cat = mot_censor.resting_cat.X${mot_thresh}_sub-$subj${cens_prev_label}_censor.1D
   
   set MaskTotal_1   = Mask.Brain.AE.sub-$subj'_'$file_info1'.deob.rb10.aw.nii.gz'
   set MaskTotal_2   = Mask.Brain.AE.sub-$subj'_'$file_info2'.deob.rb10.aw.nii.gz'
   set MaskTotal_cat = Mask.Brain.AE.sub-$subj'_resting_cat.aw.nii.gz'
   
#    if (! -e $EpiPath/$file_cat.nii.gz) then
      #--- Connectivity ---
     foreach seed (`cat $SeedFile`)
	 echo "---Connectivity ($seed) BOTH runs---"  |& tee -a $LOG_FILE
	 
	 #concatenate seed time series 
# 	 if (! -e $EpiPath/ts.$seed.$file_cat.1D) then
	 cat $EpiPath/ts.$seed.$file_1.1D $EpiPath/ts.$seed.$file_2.1D > $EpiPath/ts.$seed.$file_cat.1D
# 	 endif
	 
	 # 	 concatenate censor files 
	 cat $EpiPath/$CensorFile_1 $EpiPath/$CensorFile_2 > $EpiPath/$CensorFile_cat

	 #combine masks
	 if (! -e $EpiPath/$MaskTotal_cat) then
            3dcalc -a $EpiPath/$MaskTotal_1 -b $EpiPath/$MaskTotal_2 -expr "a*b" -prefix $EpiPath/$MaskTotal_cat
	 endif
	 
	 #combine scans
	 if (! -e $EpiPath/$file_cat.nii.gz) then
            3dTcat -prefix $EpiPath/$file_cat.nii.gz -session $EpiPath $EpiPath/$file_1.nii.gz $EpiPath/$file_2.nii.gz
	 endif 
	 
	 	 #determine concat points
	    set NN = `wc -l < $EpiPath/$CensorFile_1`
	 	echo "0 $NN" > $EpiPath/concat_runs.1D

	 
# 	 if (! -e $EpiPath/Fim.$seed.$file_cat.nii.gz) then
	 3dDeconvolve \
   	    -input $EpiPath/$file_1.nii.gz $EpiPath/$file_2.nii.gz \
   	    -polort 3 \
	    #-concat $EpiPath/concat_runs.1D \
   	    -mask $EpiPath/$MaskTotal_cat \
   	    -num_stimts 1 \
   	    -stim_file 1 $EpiPath/ts.$seed.$file_cat.1D -stim_label 1 seed.$seed \
   	    -fout -tout -rout \
   	    -bucket $EpiPath/Fim.$seed.$file_cat.nii.gz \
   	    -censor $EpiPath/$CensorFile_cat 
# 	 endif

      end 
#    endif #both exist
   
   endif #CONN_BOTH_RUNS
    
end #dir
