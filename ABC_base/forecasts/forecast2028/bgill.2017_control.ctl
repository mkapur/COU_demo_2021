#V3.20b
# star36.ctl, .dat - as star35, but with all junk code, comments, etc deleted ("clean") for the final document
# bgill.2017fixmat fixes the error in the maturity parameter slope- from 0.031 to 0.31 (units problem- my stupid/bad!!)


#C spawner-recruitment bias adjustment Not tuned For optimality
1 #_N_Growth_Patterns
1 #_N_Morphs_Within_GrowthPattern
#_Cond 1 #_Morph_between/within_stdev_ratio (no read if N_morphs=1)
#
2 #_Nblock_Patterns
#_Cond 
1 1 #_blocks_per_pattern
2000 2016 # begin and end years of blocks
2011 2016 #
#
0.5 #_fracfemale
0 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate
  #_no additional input for selected M option; read 1P per morph
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=not implemented; 4=not implemented
6 #_Growth_Age_for_L1
60 #_Growth_Age_for_L2 (999 to use as Linf)
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A)
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=read fec and wt from wtatage.ss
#_placeholder for empirical age-maturity by growth pattern
1 #_First_Mature_Age
1 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b
0 #_hermaphroditism option:  0=none; 1=age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
2 #_env/block/dev_adjust_method (1=standard; 2=logistic transform keeps in base parm bounds; 3=standard w/ no bound check)

#_growth_parms
#_LO 	HI 	INIT 	PRIOR 	PR_type SD 	PHASE env-var use_dev dev_minyr dev_maxyr dev_stddev Block Block_Fxn
0.02 	0.15 	0.063 	0.057 	0 	0.013 	-5 	0 0 0 0 0 0 0 # NatM_p_1_Fem_GP_1
2       32      12      13       0       99      -2      0 0 0 0 0.5 0 0 # F_Lmin
32      70      52      49      0       99     2      0 0 0 0 0.5 0 0 # F_Lmax
0.01    0.1     0.04    0.035    0       99      2      0 0 0 0 0.5 0 0 # F_VBK 
0.02    0.5     0.15     0.1     0       99      2      0 0 0 0 0.5 0 0 # F_CV-young
0.02    0.25    0.1     0.1     0       99       2      0 0 0 0 0.5 0 0 # F_CV-young
0.02 	0.25 	0.065 	0.058 	0 	0.013  -5 	0 0 0 0 0 0 0 # NatM_p_1_Mal_GP_1
2	45 	12 	9 	0 	99 	-3 	0 0 0 0 0 0 0 # L_at_Amin_Mal_GP_1
30 	60 	48.52 	43 	0 	99 	2 	0 0 0 0 0 0 0 # L_at_Amax_Mal_GP_1
0.02 	0.25 	0.046 	0.09 	0 	99 	2 	0 0 0 0 0 0 0 # VonBert_K_Mal_GP_1 
0.02 	0.75 	0.15 	0.1 	0 	99	2 	0 0 0 0 0 0 0 # CV_young_Mal_GP_1
0.02 	0.25 	0.1 	0.1 	0 	99 	2 	0 0 0 0 0 0 0 # CV_old_Mal_GP_1
-3 	3  1.132e-005 1.01e-005 	-1 	0.8 	-3 	0 0 0 0 0 0 0 # Wtlen_1_Fem
-3 	4 	3.1006 	3.12 	-1 	0.8 	-3 	0 0 0 0 0 0 0 # Wtlen_2_Fem
10 	60 	33.4 	32 	-1 	0.8 	-3 	0 0 0 0 0 0 0 # Mat50%_Fem
-3 	3 	-0.35 	-0.02 	-1 	0.8 	-3 	0 0 0 0 0 0 0 # Mat_slope_Fem
-3 	300 	159.473 	1 	-1 	0.8 	-3 	0 0 0 0	0 0 0 # Eggs/kg_inter_Fem
-3 	300 	95.00 	0 	-1 	0.8 	-3 	0 0 0 0 0 0 0 # Eggs/kg_slope_wt_Fem
-3 	3 1.132e-005 1.01e-005 	-1 	0.8 	-3 	0 0 0 0 0 0 0 # Wtlen_1_mal
-3 	4 	3.1006 	3.12 	-1 	0.8 	-3 	0 0 0 0 0 0 0 # Wtlen_2_mal 

# fecundity relationship 124637x + 74100

 0 0 0 0 -1 0 -4 0 0 0 0 0 0 0 # RecrDist_GP_1
 0 0 0 0 -1 0 -4 0 0 0 0 0 0 0 # RecrDist_Area_1
 0 0 0 0 -1 0 -4 0 0 0 0 0 0 0 # RecrDist_Seas_1
 0 0 0 0 -1 0 -4 0 0 0 0 0 0 0 # CohortGrowDev
#
#
#_seasonal_effects_on_biology_parms
 0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
#_Cond -4 #_MGparm_Dev_Phase
#
#_Spawner-Recruitment
3 #_SR_function: 1=B-H_flattop; 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=Shepard_3Parm
#_LO 	HI 	INIT 	PRIOR 	PR_type SD 	PHASE
6 	10 	8.1 	8.3 	-1 	10 	1  # SR_R0
0.2 	1.0 	0.718 	0.718 	2 	0.158 	-5  #_steepness
0 	2 	0.5 	0.5 	-1 	0.8 	-4 # SR_sigmaR
-5 	5 	0.1 	0 	-1 	1 	-3 # SR_envlink
-5 	5 	0 	0 	-1	1 	-4 # SR_R1_offset
0 	0 	0 	0 	-1 	0 	-99 # SR_autocorr
0 #_SR_env_link
0 #_SR_env_target_0=none;1=devs;_2=R0;_3=steepness
0 #do_recdev:  0=none; 1=devvector; 2=simple deviations
1970 # first year of main recr_devs; early devs can preceed this era
2016 # last year of main recr_devs; forecast devs start in following year
5 #_recdev phase
1 # (0/1) to read 13 advanced options
 0 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
 -4 #_recdev_early_phase
 0 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1 #_lambda for Fcast_recr_like occurring before endyr+1
 1900 #_last_early_yr_nobias_adj_in_MPD
 1970 #_first_yr_fullbias_adj_in_MPD
 2016 #_last_yr_fullbias_adj_in_MPD
 2016 #_first_recent_yr_nobias_adj_in_MPD
 1 #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)
 0 #_period of cycles in recruitment (N parms read below)
 -2 #min rec_dev
  2 #max rec_dev
  0 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
#Fishing Mortality info
0.3 # F ballpark for tuning early phases
-2001 # F ballpark year (neg value to disable)
3 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
2.9 # max F or harvest rate, depends on F_Method
# no additional F input needed for Fmethod 1
# if Fmethod=2; read overall start F value; overall phase; N detailed inputs to read
# if Fmethod=3; read N iterations for tuning for Fmethod 3
4  # N iterations for tuning F in hybrid method (recommend 3 to 7)
#
#_initial_F_parms
#_LO HI INIT PRIOR PR_type SD PHASE
 0 1 0 0.01 0 99 -1 # Impl_err_2002
 0 1 0 0.01 0 99 -1 # Impl_err_2002
 0 1 0 0.01 0 99 -1 # Impl_err_2002
#
#_Q_setup
 # Q_type options:  <0=mirror, 0/1=float, 2=parameter, 3=parm_w_random_dev, 4=parm_w_randwalk)
 #_Den-dep  env-var  extra_se  Q_type
 0 0 0 0 # 1 FISHERY1
 0 0 0 0 # 1 FISHERY2
 0 0 0 0 # 1 FISHERY3 
 0 0 0 0 # 2 SURVEY1
 0 0 0 0 # 3 SURVEY2
 0 0 0 0 # 3 SURVEY3
 0 0 0 0 # 3 SURVEY4
 0 0 0 0 # 3 SURVEY5
 0 0 0 0 # 3 SURVEY6
 0 0 0 0 # 3 SURVEY7

#
#_Cond 0 #_If q has random component, then 0=read one parm for each fleet with random q; 1=read a parm for each year of index
#_Q_parms(if_any)
# LO HI INIT PRIOR PR_type SD PHASE
#
#_size_selex_types24 is double normal
#_Pattern Discard Male Special
 24 0 0 0 # 1 FISHERY1
 24 0 0 0 # 1 FISHERY2
 24 0 0 0 # 1 FISHERY3
 1 0 0 0 # 2 SURVEY1
 5 0 0 4 # 3 SURVEY2
 1 0 0 0 # 3 SURVEY3
 5 0 0 1 # 3 SURVEY4
 5 0 0 2 # 3 SURVEY5
 5 0 0 3 # 3 SURVEY6
 5 0 0 6 # 3 SURVEY7
#
#_age_selex_types
#_Pattern ___ Male Special
 10 0 0 0 # 1 FISHERY1
 10 0 0 0 # 1 FISHERY 
 10 0 0 0 # 1 FISHERY
 10 0 0 0 # 2 SURVEY1
 10 0 0 0 # 3 SURVEY2
 10 0 0 0 # 1 SURVEY3
 10 0 0 0 # 1 SURVEY4
 10 0 0 0 # 1 SURVEY5
 10 0 0 0 # 1 SURVEY6
 10 0 0 0 # 1 SURVEY7
#_LO HI INIT PRIOR PR_type SD PHASE env-var use_dev dev_minyr dev_maxyr dev_stddev Block Block_Fxn
# size sel for south.fixed, double normal- but ascending only (pattern 24)
20	60	46	48	-1	10	3	0 	0 	0 	0 	0.5 	1 	0 # peak
-15	24	13	13	-1	10	-1	0 	0 	0 	0 	0.5 	0 	0 # init
-2	9	4	5	-1	10	4	0 	0 	0 	0 	0.5 	0 	0 # infl
-5 	20	11	5	-1	10	-2	0 	0 	0 	0 	0.5 	0 	0 # slope1
-20	1	-2	-5	-1	10	4	0 	0 	0 	0 	0.5 	0 	0 # final
-9	19	10	10	-1	10	-2	0 	0 	0 	0 	0.5 	0 	0 # final

# size sel for Central.fixed, double normal (pattern 24)
20	60	45	40	-1	10	3	0 	0 	0 	0 	0.5 	0 	0 # peak
-15	24	10	10	-1	10	-1	0 	0 	0 	0 	0.5 	0 	0 # init
-2	9	4	5	-1	10	4	0 	0 	0 	0 	0.5 	0 	0 # infl
-5 	20	11	5	-1	10	-2	0 	0 	0 	0 	0.5 	0 	0 # slope1
-20	1	-2	-5	-1	10	4	0 	0 	0 	0 	0.5 	0 	0 # final
-9	19	10	10	-1	10	-2	0 	0 	0 	0 	0.5 	0 	0 # final

# size sel for central trawl- double normal
20	60	45	40	-1	10	3	0 	0 	0 	0 	0.5 	2 	0 # peak
-15	24	10	10	-1	10	-1	0 	0 	0 	0 	0.5 	0 	0 # init
-2	9	4	5	-1	10	4	0 	0 	0 	0 	0.5 	0 	0 # infl
-5 	20	11	5	-1	10	-2	0 	0 	0 	0 	0.5 	0 	0 # slope1
-20	1	-2	-5	-1	10	4	0 	0 	0 	0 	0.5 	0 	0 # final
-9	19	10	10	-1	10	-2	0 	0 	0 	0 	0.5 	0 	0 # final

# triennial- logistic
20 	60      45      40      0       99        3      0       0       0       0       0.5     0       0       #peak
0.001   20     5.0	6.0     0       99        3      0       0       0       0       0.5     0       0       #init

# mirror sel. for NWFSC slope survey to triennial (same latitude range) 
0        20       1      1     0       99       -3      0       0       0       0       0.5     0       0       #peak
20       30       30 	 30    0       99        -2      0       0       0       0       0.5     0       0       #init

# size sel for NWC combined shelf/slope survey- logistic
16 	60      45      40      0       99        3      0       0       0       0       0.5     0       0       #peak
0.001   20     5.0	6.0     0       99        3      0       0       0       0       0.5     0       0       #init

# mirror sel. for ghost1 (south fixed)
0        20       1      1     0       99       -3      0       0       0       0       0.5     0       0       #peak
20       30       30 	 30    0       99        -2      0       0       0       0       0.5     0       0       #init
# mirror sel. for ghost2 (cen fixed)
0        20       1      1     0       99       -3      0       0       0       0       0.5     0       0       #peak
20       30       30 	 30    0       99        -2      0       0       0       0       0.5     0       0       #init
# mirror sel. for ghost1 (trawl fishery)
0        20       1      1     0       99       -3      0       0       0       0       0.5     0       0       #peak
20       30       30 	 30    0       99        -2      0       0       0       0       0.5     0       0       #init
# mirror sel. for ghost1 (combo survey)
0        20       1      1     0       99       -3      0       0       0       0       0.5     0       0       #peak
20       30       30 	 30    0       99        -2      0       0       0       0       0.5     0       0       #init

#_Cond 0 #_custom_sel-env_setup (0/1)
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no enviro fxns
#_Cond 
1 #_custom_sel-blk_setup (0/1)
-2 0 0 -0.1 0 99 4 #_placeholder when no block usage
-2 2 0 0.1 0 99 4 #_placeholder when no block usage
#-2 2 0 0.1 0 99 -4 #_placeholder when no block usage
#_Cond No selex parm trends
#_Cond 
# placeholder for selparm_Dev_Phase
#_Cond 
1 #_env/block/dev_adjust_method (1=standard; 2=logistic trans to keep in base parm bounds; 3=standard w/ no bound check)
#
# Tag loss and Tag reporting parameters go next
0  # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
1 #_Variance_adjustments_to_input_values - harmonic for LFs, CAALs
#fleet1	2	3	4	5	6	7	8	9	10	
0	0	0	0.101197	0	0	0	0	0	0	#_add_to_survey_CV
0	0	0	0	0	0	0	0	0	0	#_add_to_discard_stddev
0	0	0	0	0	0	0	0	0	0	#_add_to_bodywt_CV
0.067	0.090	0.206	0.354	1	0.125	1	1	1	1	#_mult_by_lencomp_N
0.036	0.022	0.223	1	1	0.0048	1	1	1	1	#_mult_by_agecomp
1	1	1	1	1	1	1	1	1	1	#_mult_by_size_age

#
5 #_maxlambdaphase
1 #_sd_offset
#
13 # number of changes to make to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch;
# 9=init_equ_catch; 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark
#like_comp fleet/survey  phase  value  sizefreq_method
 4 1 1 1 1
 1 5 1 1 1
 1 6 1 1 1
 4 7 1 0 1 
 4 8 1 0 1
 4 9 1 0 1
 5 7 1 0 1
 5 8 1 0 1
 5 9 1 0 1
 5 10 1 0 1
17 1 1 0 1
17 2 1 0 1
17 3 1 0 1
 
# 4 2 3 1 1
#
0 # (0/1) read specs for more stddev reporting
999
