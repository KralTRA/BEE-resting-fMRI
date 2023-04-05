#!/bin/tcsh

#---
# Preprocessing script for BEE rs-fMRI data
#---
#---
# Author: R.M. Birn; adapted by T.R.A. Kral
#---

#---------------------------
# Paths:
#---------------------------
set RootPath   = /study/sjshort_ebds
set ScriptPath = $RootPath/rsfc/scripts
set MergePath  = $RootPath/rsfc/MERGE
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
set BANDPASS       = 0
set CLEANUP        = 0
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
# 1067 & 1168 have 4 runs; 1042, 1112, 1187 & 1131 have 3 runs (PA run4)
      
#---------------------------
# Loop:
#---------------------------
foreach dir ($argv) #(`cat $ScriptPath/dirs_good_0.4wP_excludeCoilIssue.txt`) #(`cat $ScriptPath/dirs_1_good.txt`) #(`cat $ScriptPath/dirs_scan2.txt`) #
   set subj = $dir

   set DataPath     = $RootPath/rsfc/processed-rmb/$subj/func
   set EpiPath		= /study/sjshort_ebds/rsfc/processed-trak/withBP/$subj
#    set EpiPath		= /scratch/sjshort_ebds/processed-trak/newTopup_ntrp_bandpass_clean/$subj
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
# 	set baseprefix = sub-$subj'_task-rest_dir-AP_run-02_bold.deob.nii.gz'
	set basesuffix = "[10]"
	
	foreach file_orig ($fmri_runs)
   
   		echo "file_orig: $file_orig"
   		set filename = $file_orig:r:r
   		set NT = `3dinfo -nt $DataPath/$file_orig`
   		if ($NT < 2) then
      		echo "---------- $file_orig has less than 2 time points ($NT)" |& tee -a $LOG_FILE
      		continue
   		else
   
      	set file0     = $filename
      
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
      
      #==================================
      # Jump ahead if files already exist
      #    (N.B. you can include multiple conditions using the syntax "A && B")
      if ($SKIP_AHEAD) then
      set goto_code = ( \
             aw.NgwcmdfX${mot_thresh}.s${fwhm}.mni 	CONNECTIVITY	$BANDPASS \
		     aw.NgwcmdX${mot_thresh}.s${fwhm}.mni 	CONNECTIVITY	1 \
		     aw.NgwcmdfX${mot_thresh}.s${fwhm} 		RESAMPLE_TO_MNI	$BANDPASS \
		     aw.NgwcmdX${mot_thresh}.s${fwhm} 		RESAMPLE_TO_MNI	1 \
		     aw.Ngwcmdf	 				SMOOTH		$BANDPASS \
		     aw.Ngwcmd    				SMOOTH		1 \
		     aw						SEGMENT		1 \
                   )

      set N = $#goto_code
      set n = 1
      while ($n < $N) 
	  set cor = $goto_code[$n]
	  @ n++
	  set g = $goto_code[$n]
	  @ n++
	  set condition = $goto_code[$n]
	  @ n++
	  if (-e $file0.$cor.nii.gz && $condition) then
	     set file = $file0.$cor
	     echo "$file0.nii.gz exists."
	     echo "jumping to $g"
	     goto $g
	  endif 
      end
      endif #SKIP_AHEAD
      
      #====================================================================================
      DEOBLIQUE:
      #------------------------------------------------------------------------------------
      #--- DeOblique ---
      set file_out = ${file0}.deob
#       if (`3dinfo -is_oblique $EpiPath/$file0.nii.gz`) then
	  if (! -e $EpiPath/$file_out.nii.gz) then
	 		if ($DEBUG) then
	 		echo "---------- Running Deoblique..."
	 		3dWarp -deoblique -prefix $EpiPath/$file_out.nii.gz $DataPath/$file0.nii.gz |& tee -a $LOG_FILE
      		endif
      else echo "Deoblique COMPLETED"
	  endif
      set file0 = $file_out
      
      #====================================================================================
	   ALIGN_CENTERS_RESAMPLE_ANAT:
	   #------------------------------------------------------------------------------------
		if (! -e $StrucPath/$AnatFile_alignC) then
			echo "------------- Aligning Center of Anatomical (T2) to EPI---"
# 			@Align_Centers -base $EpiPath/$file0.nii.gz -prefix $AnatFile_atlasC -dset $StrucPath/awpy/$AnatFile_atlas
			@Align_Centers -base $EpiPath/$file0.nii.gz -cm -prefix $AnatFile_alignC -dset $StrucPath/$AnatFile_ns # did we use the cm option?
		else echo "Align COMPLETED"
		endif

		if (! -e $StrucPath/$AnatFileRS) then
			echo "------------- Resample Anatomical (T2) to EPI---"
# 			3dresample -master $EpiPath/$file0.nii.gz -prefix $AnatFile_atlasRS -input $StrucPath/awpy/$AnatFile_atlasC
			3dresample -master $EpiPath/$file0.nii.gz -prefix $StrucPath/$AnatFileRS -input $StrucPath/$AnatFile_alignC
		else echo "Resample COMPLETED"
		endif
      
    #====================================================================================
      REGISTER_ANAT_TO_TEMPLATE:  
      #------------------------------------------------------------------------------------
      #--- Register anat to UNC template ---
      set anat = $AnatFileRS
      set anatfilename = $anat:r
      
      cd $StrucPath

	if (! -e $StrucPath/awpy/$AnatFile_atlas) then
		if (-e $StrucPath/awpy/anat.aff.nii) then
			rm -rf $StrucPath/awpy/
		endif
		echo "---------- Registering Anatomical (T2) to Template---"
		auto_warp.py -base $RootPath/UNCInfant012Atlas_20140325/infant-neo-withCerebellum_epiC-RS.nii \
		-input $anat -skull_strip_input no \
		-unifize_input no
	else
		echo "$dir T2 registration complete"
	endif

   
     #====================================================================================
      SEGMENT:
      #------------------------------------------------------------------------------------
      #--- Segment Anatomical ---
      set anat = $AnatFile_atlas
      set anatfilename = $anat:r
      if (-e $StrucPath/segment_${anatfilename}_seg_0.nii.gz || -e $StrucPath/Mask.CSF.nii.gz) then
	 echo "---------- FSL segmentation already exists."  |& tee -a $LOG_FILE
      else
	 echo "---------- starting FSL segmentation..."  |& tee -a $LOG_FILE
	 fast -t 1 -g -p -o $StrucPath/segment_$anat $StrucPath/awpy/$anat
	 echo "...done"
      endif
      
      if (! -e $StrucPath/Mask.CSF.nii.gz) then
      3dcopy $StrucPath/segment_${anatfilename}_seg_0.nii.gz $StrucPath/Mask.CSF.nii.gz
      3dcopy $StrucPath/segment_${anatfilename}_seg_1.nii.gz $StrucPath/Mask.GM.nii.gz
      3dcopy $StrucPath/segment_${anatfilename}_seg_2.nii.gz $StrucPath/Mask.WM.nii.gz
      endif

      #--- Resample ---
      if (! -e $StrucPath/Mask.CSF.${dxyz}mm.nii.gz) then
      3dresample -dxyz ${dxyz} ${dxyz} ${dxyz} -prefix $StrucPath/Mask.CSF.${dxyz}mm.nii.gz -input $StrucPath/Mask.CSF.nii.gz
      3dresample -dxyz ${dxyz} ${dxyz} ${dxyz} -prefix $StrucPath/Mask.GM.${dxyz}mm.nii.gz  -input $StrucPath/Mask.GM.nii.gz
      3dresample -dxyz ${dxyz} ${dxyz} ${dxyz} -prefix $StrucPath/Mask.WM.${dxyz}mm.nii.gz  -input $StrucPath/Mask.WM.nii.gz
      endif

      set MaskCSF = Mask.CSF.${dxyz}mm
      set MaskGM  = Mask.GM.${dxyz}mm
      set MaskWM  = Mask.WM.${dxyz}mm

      if (! -e $StrucPath/$MaskCSF.nii.gz) then
         echo "---------- Could not create segmentation."
	 continue
      endif
      
      endif
      end
      
	foreach file_orig ($fmri_runs)
   
   		echo "file_orig: $file_orig"
   		set filename = $file_orig:r:r

   		set NT = `3dinfo -nt $DataPath/$file_orig`
   		if ($NT < 2) then
      		echo "---------- $file_orig has less than 2 time points ($NT)" |& tee -a $LOG_FILE
      		continue
   		else
   
      	set file0     = $filename.deob
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
      
      #==================================
      # Jump ahead if files already exist
      #    (N.B. you can include multiple conditions using the syntax "A && B")
      if ($SKIP_AHEAD) then
      set goto_code = ( \
             aw.NgwcmdfX${mot_thresh}.s${fwhm}.mni 	CONNECTIVITY	$BANDPASS \
		     aw.NgwcmdX${mot_thresh}.s${fwhm}.mni 	CONNECTIVITY	1 \
		     aw.NgwcmdfX${mot_thresh}.s${fwhm} 		RESAMPLE_TO_MNI	$BANDPASS \
		     aw.NgwcmdX${mot_thresh}.s${fwhm} 		RESAMPLE_TO_MNI	1 \
		     aw.Ngwcmdf	 				SMOOTH		$BANDPASS \
		     aw.Ngwcmd    				SMOOTH		1 \
		     aw						SEGMENT		1 \
                   )

      set N = $#goto_code
      set n = 1
      while ($n < $N) 
	  set cor = $goto_code[$n]
	  @ n++
	  set g = $goto_code[$n]
	  @ n++
	  set condition = $goto_code[$n]
	  @ n++
	  if (-e $file0.$cor.nii.gz && $condition) then
	     set file = $file0.$cor
	     echo "$file0.nii.gz exists."
	     echo "jumping to $g"
	     goto $g
	  endif 
      end
      endif #SKIP_AHEAD

	  #====================================================================================
      VOLREG:
      #------------------------------------------------------------------------------------
      #--- motion correction ---
      set file_out = $file0.rb$ign
      if (! -e $EpiPath/$file_out.nii.gz) then
      3dvolreg -base "$EpiPath/$baseprefix$basesuffix" -1Dfile $EpiPath/$motfile -prefix $EpiPath/$file_out.nii.gz $EpiPath/$file0.nii.gz'['$ign'-$]' |& tee -a $LOG_FILE
      endif
      #--- if failed, try resampling to base or 3dAllineate (e.g. if grids don't match)---
      if (! -e $EpiPath/$file_out.nii.gz) then
#       3dAllineate -base $EpiPath/$basefile -master $EpiPath/$file0.nii.gz -warp shr -1Dparam_save $EpiPath/$motfile -prefix $EpiPat/$file_out.nii.gz -input $EpiPath/$file0.nii.gz'['$ign'-$]' 
      	3dresample -base $EpiPath/$baseprefix -input $EpiPath/$file0.nii.gz -prefix $EpiPath/$file0.rs.nii.gz
      	3dvolreg -base "$EpiPath/$baseprefix$basesuffix" -1Dfile $EpiPath/$motfile -prefix $EpiPath/$file_out.nii.gz $EpiPath/$file0.rs.nii.gz'['$ign'-$]' |& tee -a $LOG_FILE
      endif
      set file0 = $file_out
      
      #====================================================================================
      AUTOMASK:
      set file_out = $file0.am
      if (! -e $EpiPath/$file0.mask.nii.gz) then
      # create mask
      3dAutomask -prefix $EpiPath/$file0.mask.nii.gz $EpiPath/$file0.nii.gz
      endif
      
      if (! -e $EpiPath/$file_out.nii.gz) then
      echo "---------- applying brain mask to EPI"
      3dcalc -a $EpiPath/$file0.nii.gz -b $EpiPath/$file0.mask.nii.gz -expr 'a*b' -prefix $EpiPath/$file_out.nii.gz
      endif
      
      #====================================================================================
      COMPUTE_MOTION_CENSOR:
      #------------------------------------------------------------------------------------
      if (! -e $EpiPath/${CensorFilePrefix}_censor.1D) then
	 echo "---Motion Censoring---"  |& tee -a $LOG_FILE
	 1d_tool.py -infile $EpiPath/$motfile -censor_motion $mot_thresh $EpiPath/$CensorFilePrefix $cens_prev_string
      endif
      
      #====================================================================================
      ALIGN_EPI_ANAT:
      #------------------------------------------------------------------------------------
      #--- Align EPI to Anatomical ---
      set cor = a

      #--- Unifize EPI (bias field correction) --- this is written for T1 scans
      if (! -e $EpiPath/$file0.am.1u.nii.gz) then
		3dUnifize -input $EpiPath/$file0.am.nii.gz'['$ign']' -EPI -prefix $EpiPath/$file0.am.1u.nii.gz
      endif
      
      set fileU = $file0.am.1u

      #--- Align one volume ---
      if (! -e $EpiPath/xform.$file0-Anat.aff12.1D) then
 	 echo "---------- Aligning $fileU to $AnatFileRS ..."  |& tee -a $LOG_FILE
	    
	    3dAllineate  -base $StrucPath/$AnatFileRS \
	    -source $EpiPath/$fileU.nii.gz \
	    -prefix $EpiPath/$file0.${cor}.nii.gz \
	    -1Dmatrix_save $EpiPath/${fileU}.${cor}_mat.aff12.1D

	 #--- Save transforms ---
	 /bin/cp -f $EpiPath/${fileU}.${cor}_mat.aff12.1D $EpiPath/xform.$file0-Anat.aff12.1D
      endif
      
      #================================
      ALIGN_EPI_TEMPLATE:      
      #================================
      set file_out = ${file0}.aw
      if (! -e $EpiPath/$file_out.nii.gz) then
	 echo "---------- Applying nonlinear warping to $file0 (master: $AnatFile_atlas.nii.gz)..."  |& tee -a $LOG_FILE

	 #--- concatenate transformation matrixes ---
	 if (! -e $StrucPath/xform.Anat-Template.aff12.1D) then
	 	cat_matvec -ONELINE $StrucPath/awpy/anat.aff.Xat.1D > $StrucPath/xform.Anat-Template.aff12.1D
	 endif
	 
	 if (! -e $EpiPath/xform.$file0-Template2.aff12.1D) then
	 	cat_matvec -ONELINE $StrucPath/xform.Anat-Template.aff12.1D $EpiPath/xform.$file0-Anat.aff12.1D > $EpiPath/xform.$file0-Template.aff12.1D
	 endif

	 #--- Apply NL warp ---
	 3dNwarpApply \
            -nwarp $StrucPath/'awpy/anat.aff.qw_WARP.nii' $EpiPath'/xform.'$file0'-Template.aff12.1D' \
            -prefix $EpiPath/$file_out.nii.gz \
            -master $StrucPath/awpy/$AnatFile_atlas \
            -source $EpiPath/$file0.nii.gz
      endif
      set file0 = $file_out

      #====================================================================================
      MASK:
      #------------------------------------------------------------------------------------
      #--- Create Mask ---
      	 echo "---------- Creating brain masks..."  |& tee -a $LOG_FILE

      set MaskTotal = Mask.Brain.AE.$file0.nii.gz
      if (! -e $StrucPath/Mask.Brain.$AnatFile_atlas) then
      echo "3dAutomask -prefix $StrucPath/Mask.Brain.$AnatFile_atlas $StrucPath/awpy/$AnatFile_atlas"
      3dAutomask -prefix $StrucPath/Mask.Brain.$AnatFile_atlas $StrucPath/awpy/$AnatFile_atlas
      endif
      
      if (! -e $EpiPath/Mask.Brain.$file0.nii.gz) then
      echo "3dAutomask -prefix $EpiPath/Mask.Brain.$file0.nii.gz $EpiPath/$file0.nii.gz"
      3dAutomask -prefix $EpiPath/Mask.Brain.$file0.nii.gz $EpiPath/$file0.nii.gz
	  endif
	  
      if (! -e $EpiPath/$MaskTotal) then
      3dcalc -a $EpiPath/Mask.Brain.$file0.nii.gz -b $StrucPath/Mask.Brain.$AnatFile_atlas -expr "step(a+b)" -prefix $EpiPath/$MaskTotal
      endif
      
     #====================================================================================
      MEAN:
      #------------------------------------------------------------------------------------
      #--- Compute mean (to add back after 3dTproject) ---
      if (! -e $EpiPath/Mean.$file0.nii.gz) then
         3dTstat -mean -prefix $EpiPath/Mean.$file0.nii.gz $EpiPath/$file0.nii.gz
      endif
      
      set MeanFile = Mean.$file0.nii.gz
      
      #====================================================================================
      NUISANCE:
      #------------------------------------------------------------------------------------     
      #--- Nuisance regression, Temporal Filtering, Spatial Smoothing ---
      if($BANDPASS) then
	 set file_out = ${file0}.NgwcmdfX${mot_thresh}
	 set bp_cmd   = "-bandpass 0.01 0.1"
      else
	 set file_out = ${file0}.NgwcmdX${mot_thresh}
	 set bp_cmd   = ""
      endif

      if (! -e $EpiPath/$file_out.nii.gz) then

      #----------------------------------------
      # Motion
      #----------------------------------------
      if (! -e $EpiPath/mot.$file0.1x.1D) then
      1dnorm -demean -normx $EpiPath/$motfile'[0]' $EpiPath/mot.$file0.1x.1D
      1dnorm -demean -normx $EpiPath/$motfile'[1]' $EpiPath/mot.$file0.2x.1D
      1dnorm -demean -normx $EpiPath/$motfile'[2]' $EpiPath/mot.$file0.3x.1D
      1dnorm -demean -normx $EpiPath/$motfile'[3]' $EpiPath/mot.$file0.4x.1D 
      1dnorm -demean -normx $EpiPath/$motfile'[4]' $EpiPath/mot.$file0.5x.1D 
      1dnorm -demean -normx $EpiPath/$motfile'[5]' $EpiPath/mot.$file0.6x.1D
      endif
      
      #----------------------------------------
      # Derivative of Motion
      #----------------------------------------
      foreach mi (1 2 3 4 5 6)
         if (! -e $EpiPath/mot.$file0.${mi}x.d.1D && ! -z $EpiPath/mot.$file0.${mi}x.d.1D) then
         1d_tool.py -derivative -infile $EpiPath/mot.$file0.${mi}x.1D -write - > $EpiPath/tmp.mot.$file0.${mi}x.d.1D 
         1dnorm -demean -normx $EpiPath/tmp.mot.$file0.${mi}x.d.1D - > $EpiPath/mot.$file0.${mi}x.d.1D
#          rm -f tmp.mot.$file0.${mi}x.d.1D
         endif
      end
      
      cd $EpiPath

      #----------------------------------------
      # Global
      #----------------------------------------
      #--- Brain mask ---
      if ($DEBUG) echo "...global..."  |& tee -a $LOG_FILE
      if (! -e Mask.Brain.$file0.nii.gz) then
      3dAutomask -prefix $EpiPath/Mask.Brain.$file0.nii.gz $EpiPath/$file0.nii.gz
      endif
      
      if (! -e $EpiPath/ts.global.$file0.1D && ! -z $EpiPath/ts.global.$file0.1D) then
      3dROIstats -mask $EpiPath/Mask.Brain.$file0.nii.gz -mask_f2short -quiet $EpiPath/$file0.nii.gz > $EpiPath/tmp.ts.global.$file0.1D
      rm -f $EpiPath/ts.global.$file0.1D
      1dnorm -demean -normx $EpiPath/tmp.ts.global.$file0.1D $EpiPath/ts.global.$file0.1D
#       rm -f $EpiPath/tmp.ts.global.$file0.1D
      endif

      #---------------------------------------------------
      # Derivative of Global
      #---------------------------------------------------
      if ($DEBUG) echo "...global derivative..."  |& tee -a $LOG_FILE
      if (! -e $EpiPath/ts.global.d.$file0.1D && ! -z $EpiPath/ts.global.d.$file0.1D) then
      1d_tool.py -derivative -infile $EpiPath/ts.global.$file0.1D -write - > $EpiPath/tmp.ts.global.d.$file0.1D 
      rm -f $EpiPath/ts.global.d.$file0.1D
      1dnorm -demean -normx $EpiPath/tmp.ts.global.d.$file0.1D $EpiPath/ts.global.d.$file0.1D
#       rm -f $EpiPath/tmp.ts.global.d.$file0.1D
      endif

      #----------------------------------------
      # WM
      #----------------------------------------
      #--- Erode mask ---
      if ($DEBUG) echo "...WM..."  |& tee -a $LOG_FILE
      if (! -e $StrucPath/$MaskWM.erode.nii.gz) then
      3dcalc \
	 -datum short \
	 -a $StrucPath/$MaskWM.nii.gz \
	 -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
	 -expr "a*(1-amongst(0,b,c,d,e,f,g))" \
	 -prefix $StrucPath/$MaskWM.erode.nii.gz
      endif
      
      if (! -e $StrucPath/ts.WMe.$file0.1D && ! -z $StrucPath/ts.WMe.$file0.1D) then
      3dROIstats -mask $StrucPath/$MaskWM.erode.nii.gz -mask_f2short -quiet $EpiPath/$file0.nii.gz > $StrucPath/tmp.ts.WMe.$file0.1D
      1dnorm -demean -normx $StrucPath/tmp.ts.WMe.$file0.1D - > $StrucPath/ts.WMe.$file0.1D
#       rm -f tmp.ts.WMe.$file0.1D
      endif

      #---------------------------------------------------
      # Derivative of WM 
      #---------------------------------------------------
      if ($DEBUG) echo "...WM derivative..."  |& tee -a $LOG_FILE
      if (! -e $StrucPath/ts.WMe.d.$file0.1D && ! -z $StrucPath/ts.WMe.d.$file0.1D) then
      1d_tool.py -derivative -infile $StrucPath/ts.WMe.$file0.1D -write - > $StrucPath/tmp.ts.WMe.d.$file0.1D 
      1dnorm -demean -normx $StrucPath/tmp.ts.WMe.d.$file0.1D - > $StrucPath/ts.WMe.d.$file0.1D
#       rm -f $StrucPath/tmp.ts.WMe.d.$file0.1D
      endif

      #----------------------------------------
      # CSF 
      #----------------------------------------
      #--- Erode mask ---
      if ($DEBUG) echo "...CSF..."  |& tee -a $LOG_FILE
      if (! -e $StrucPath/$MaskCSF.erode.nii.gz) then
      3dcalc \
	 -datum short \
	 -a $StrucPath/$MaskCSF.nii.gz \
	 -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
	 -expr "a*(1-amongst(0,b,c,d,e,f,g))" \
	 -prefix $StrucPath/$MaskCSF.erode.nii.gz
      endif
      
      if (! -e $StrucPath/ts.CSFe.$file0.1D && ! -z $StrucPath/ts.CSFe.$file0.1D) then
      3dROIstats -mask $StrucPath/$MaskCSF.erode.nii.gz -mask_f2short -quiet $EpiPath/$file0.nii.gz > $StrucPath/tmp.ts.CSFe.$file0.1D
      1dnorm -demean -normx $StrucPath/tmp.ts.CSFe.$file0.1D - > $StrucPath/ts.CSFe.$file0.1D
#       rm -f $StrucPath/tmp.ts.CSFe.$file0.1D
      endif

      #----------------------------------------------------
      # Derivative of CSF 
      #----------------------------------------------------
      if ($DEBUG) echo "...CSF derivative..."  |& tee -a $LOG_FILE
      if (! -e $StrucPath/ts.CSFe.d.$file0.1D && ! -z $StrucPath/ts.CSFe.d.$file0.1D) then
      1d_tool.py -derivative -infile $StrucPath/ts.CSFe.$file0.1D -write - > $StrucPath/tmp.ts.CSFe.d.$file0.1D 
      1dnorm -demean -normx $StrucPath/tmp.ts.CSFe.d.$file0.1D - > $StrucPath/ts.CSFe.d.$file0.1D
#       rm -f $StrucPath/tmp.ts.CSFe.d.$file0.1D
      endif

      echo "---------- nuisance regression"  |& tee -a $LOG_FILE
      	 
      cd $EpiPath

      3dTproject \
	 -input $EpiPath/$file0.nii.gz \
	 -prefix $EpiPath/$file_out.nii.gz \
	 -mask $EpiPath/$MaskTotal \
	 -censor $EpiPath/$CensorFile \
	 -cenmode NTRP \
	 $bp_cmd \
	 -ort $EpiPath/ts.global.$file0.1D \
	 -ort $EpiPath/ts.global.d.$file0.1D \
	 -ort $StrucPath/ts.WMe.$file0.1D \
	 -ort $StrucPath/ts.WMe.d.$file0.1D \
	 -ort $StrucPath/ts.CSFe.$file0.1D \
	 -ort $StrucPath/ts.CSFe.d.$file0.1D \
	 -ort $EpiPath/mot.$file0.1x.1D \
	 -ort $EpiPath/mot.$file0.2x.1D \
	 -ort $EpiPath/mot.$file0.3x.1D \
	 -ort $EpiPath/mot.$file0.4x.1D \
	 -ort $EpiPath/mot.$file0.5x.1D \
	 -ort $EpiPath/mot.$file0.6x.1D \
	 -ort $EpiPath/mot.$file0.1x.d.1D \
	 -ort $EpiPath/mot.$file0.2x.d.1D \
	 -ort $EpiPath/mot.$file0.3x.d.1D \
	 -ort $EpiPath/mot.$file0.4x.d.1D \
	 -ort $EpiPath/mot.$file0.5x.d.1D \
	 -ort $EpiPath/mot.$file0.6x.d.1D 
	 	       
      endif
      set file0 = $file_out

      #====================================
      BANDPASS:      
      #====================================
      echo "---Bandpass Filtering---"  |& tee -a $LOG_FILE
      #--- bandpass filtering (optional; use if insufficient DOF with 3dTproject) ---
      set file_out = ${file0}.bp
      if (! -e $EpiPath/$file_out.nii.gz ) then
      3dBandpass -prefix $EpiPath/tmp.$file_out.nii.gz 0.01 0.1 $EpiPath/$file0.nii.gz
           
      #--- add mean back ---
      3dcalc -a $EpiPath/tmp.$file_out.nii.gz -b $EpiPath/$MeanFile -expr "a+b" -prefix $EpiPath/$file_out.nii.gz

      endif
      set file0 = $file_out

      #====================================
      SMOOTH:      
      #====================================
      echo "---Smoothing---"  |& tee -a $LOG_FILE
      #--- spatial smoothing (optional) ---
      set file_out = ${file0}.s${fwhm}
      if (! -e $EpiPath/$file_out.nii.gz) then
	 3dmerge -doall -1blur_fwhm $fwhm -prefix $EpiPath/$file_out.nii.gz $EpiPath/$file0.nii.gz
      endif
      set file0 = $file_out
      
      #====================================================================================
      CONNECTIVITY:
      #------------------------------------------------------------------------------------
      #--- Connectivity ---
     foreach seed (`cat $SeedFile`)

	 if (! -e $EpiPath/ts.$seed.$file0.1D ) then
	 echo "---Extracting data from ($seed)---"  |& tee -a $LOG_FILE
	 3dmaskave -mask $RoiPath/$seed.nii.gz -quiet $EpiPath/$file0.nii.gz > $EpiPath/tmp.ts.$seed.$file0.1D
	 1dnorm -demean $EpiPath/tmp.ts.$seed.$file0.1D - > $EpiPath/ts.$seed.$file0.1D
	 endif
	 
	 endif
 end # fmri runs ?
   
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
#    	-input $EpiPath/$file_cat.nii.gz \
   	    -input $EpiPath/$file_1.nii.gz $EpiPath/$file_2.nii.gz \
   	    -polort 3 \
	    #-concat $EpiPath/concat_runs.1D \
   	    -mask $EpiPath/$MaskTotal_cat \
   	    -num_stimts 1 \
   	    -stim_file 1 $EpiPath/ts.$seed.$file_cat.1D -stim_label 1 seed.$seed \
   	    -fout -tout -rout \
   	    -bucket $EpiPath/Fim.$seed.$file_cat.nii.gz \
   	    -censor $EpiPath/$CensorFile_cat # unless censor points removed with 3dTproject cenmode KILL
# 	 endif

      end 
#    endif #both exist
   
   endif #CONN_BOTH_RUNS
    
end #dir
