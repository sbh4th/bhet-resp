////////////////////////////////////////////////////////////////////////////////

** Before running this code:
** (1) Ensure you have run all previous scripts to compile master data in accordance with BHET master data SOP and saved the resulting dataset to your desktop
** (2) Download "BHET_coal-ban-tracking_23Feb22", Lung_S1_to-merge", "Lung_S2_to-merge" and "BHET_PE_S4_29Mar23.dta" from OSF and save them to your computer
** (3) Set up your file paths accordingly: use "Search and replace" to find all instances of "/Users/leonasiaw/Desktop" within this do-file and edit them to reflect your file storage settings

////////////////////////////////////////////////////////////////////////////////

	use "/Users/leona/Desktop/BHET_master_data_2023May22.dta", clear

**	Merge in lung data from S1 and S2

*	Merge in lung data:
	merge m:1 ptc_id wave using "/Users/leona/Desktop/Lung_S1_to-merge.dta"
	rename _merge _mergeS1
	foreach v of varlist PEF_measure1-FVC_measure3{
	replace `v'=`v'_S1 if _mergeS1==3
	}
	merge m:1 ptc_id wave using "/Users/leona/Desktop/Lung_S2_to-merge.dta"
	rename _merge _mergeS2
	foreach v of varlist PEF_measure1-FVC_measure3{
	replace `v'=`v'_S2 if _mergeS2==3
	}
	
	drop if _mergeS1==2
	drop if _mergeS2==2
	drop dup _mergeS2 _mergeS1 *_S1 *_S2
		
*	Populate LF_trouble based on XY's 15/3/23 email:
	replace LF_trouble=0 if wave==1 & PEF_measure1!=.
	replace LF_trouble=0 if wave==2 & PEF_measure1!=.
	replace LF_trouble=1 if hh_id=="HH20180758" & wave==1
	replace LF_trouble=1 if hh_id=="HH20180550" & wave==1
	replace LF_trouble=1 if hh_id=="HH20180438" & wave==1
	replace LF_trouble=1 if hh_id=="HH20180783" & wave==1
	//	N.B.: Xiaoying said we'll wait to correct LF outliers
	
////////////////////////////////////////////////////////////////////////////////
	
**	Master data corrections (15Feb2023)

*	Merge in coal ban status:
	merge m:m ID_VILLAGE using "/Users/leona/Desktop/BHET_coal-ban-tracking_23Feb22.dta", keepusing(Ban_2018-Not_yet)
	drop _merge

*	Add S4 value label to "wave":
	label define wave 1 "Season 1 (2018-19)" 2 "Season 2 (2019-20)" 3 "Season 3 (2020-21)" 4 "Season 4 (2021-22)", replace

*	Fix typo:
	rename caffine caffeine

////////////////////////////////////////////////////////////////////////////////

**	Season 4 BP corrections

*	systolic_brachial1 is 10 in master, 100 in xml:
	replace systolic_brachial1=100 if ptc_id=="PT2018059401" & wave==4
*	diastolic_brachial4 is 8 in master, 81 in xml:
	replace diastolic_brachial4=81 if ptc_id=="PT2018068500" & wave==4
*	systolic_central1; diastolic_central1 both 0 in master, 116;84 in xml:
	replace systolic_central1=116 if ptc_id=="PT2018031800" & wave==4
	replace diastolic_central1=84 if ptc_id=="PT2018031800" & wave==4
*	systolic_central1; diastolic_central1 both 0 in master, 127;76 in xml:
	replace systolic_central1=127 if ptc_id=="PT2018031901" & wave==4
	replace diastolic_central1=76 if ptc_id=="PT2018031901" & wave==4
*	systolic_central1 is 724 in master, 124 in xml:
	replace systolic_central1=124 if ptc_id=="PT2018039001" & wave==4
*	systolic_central1 is 728 in master, 128 in xml:
	replace systolic_central1=128 if ptc_id=="PT2018119500" & wave==4
*	systolic_central1 is 773 in master, 113 in xml:
	replace systolic_central1=113 if ptc_id=="PT2018053004" & wave==4
*	systolic_central2 is 779 in master, 119 in xml:
	replace systolic_central2=113 if ptc_id=="PT2018034200" & wave==4
*	systolic_central3 is 21 in master, 121 in xml:
	replace systolic_central3=121 if ptc_id=="PT2019034400" & wave==4
*	systolic_central4 is 53 in master, 142 in xml:
	replace systolic_central4=142 if ptc_id=="PT2018113200" & wave==4
*	diastolic_central4 is 8 in master, 98 in xml:
	replace diastolic_central4=198 if ptc_id=="PT2018037600" & wave==4
*	diastolic_central4 is 790 in master, 90 in xml:
	replace diastolic_central4=90 if ptc_id=="PT2018086000" & wave==4
*	sd_bpID5 recorded incorrectly for PT2019022800 (fifth measurements correspond to sd_bpID 422):
	replace sd_bpID5=422 if ptc_id=="PT2019022800" & wave==4

////////////////////////////////////////////////////////////////////////////////

**	Recalculate average BP measurements (code from Brian)

// calculate average SB from last two readings
	* create row means of all combinations
	egen SB45 = rowmean(systolic_brachial4 systolic_brachial5) 
	egen SB34 = rowmean(systolic_brachial3 systolic_brachial4) 
	egen SB23 = rowmean(systolic_brachial2 systolic_brachial3) 
	egen SB12 = rowmean(systolic_brachial1 systolic_brachial2) 
	gen SB1 =  systolic_brachial1 

	* replace with correct rowmean value based on how many BP measurements taken:
	gen avg_SB = SB45
	replace avg_SB = SB34 if systolic_brachial5 == . 
	replace avg_SB = SB23 if systolic_brachial5 == . & systolic_brachial4 == . 
	replace avg_SB = SB12 if systolic_brachial5 == . & systolic_brachial4 == . & systolic_brachial3 == .  	
	replace avg_SB = SB1  if systolic_brachial5 == . & systolic_brachial4 == . & systolic_brachial3 == . & systolic_brachial2 == . 
										
// calculate average DB from last two readings
	* create row means of all combinations
	egen DB45 = rowmean(diastolic_brachial4 diastolic_brachial5) 
	egen DB34 = rowmean(diastolic_brachial3 diastolic_brachial4) 
	egen DB23 = rowmean(diastolic_brachial2 diastolic_brachial3) 
	egen DB12 = rowmean(diastolic_brachial1 diastolic_brachial2) 
	gen DB1 =  diastolic_brachial1 

	* replace with correct rowmean value based on how many BP measurements taken:
	gen avg_DB = DB45
	replace avg_DB = DB34 if diastolic_brachial5 == . 
	replace avg_DB = DB23 if diastolic_brachial5 == . & diastolic_brachial4 == . 
	replace avg_DB = DB12 if diastolic_brachial5 == . & diastolic_brachial4 == . & diastolic_brachial3 == .  	
	replace avg_DB = DB1  if diastolic_brachial5 == . & diastolic_brachial4 == . & diastolic_brachial3 == . & diastolic_brachial2 == . 

// calculate average SC from last two readings
	* create row means of all combinations
	egen SC45 = rowmean(systolic_central4 systolic_central5) 
	egen SC34 = rowmean(systolic_central3 systolic_central4) 
	egen SC23 = rowmean(systolic_central2 systolic_central3) 
	egen SC12 = rowmean(systolic_central1 systolic_central2) 
	gen SC1 =  systolic_central1 

	* replace with correct rowmean value based on how many BP measurements taken:
	gen avg_SC = SC45
	replace avg_SC = SC34 if systolic_central5 == . 
	replace avg_SC = SC23 if systolic_central5 == . & systolic_central4 == . 
	replace avg_SC = SC12 if systolic_central5 == . & systolic_central4 == . & systolic_central3 == .  	
	replace avg_SC = SC1  if systolic_central5 == . & systolic_central4 == . & systolic_central3 == . & systolic_central2 == . 
										
// calculate average DC from last two readings
	* create row means of all combinations
	egen DC45 = rowmean(diastolic_central4 diastolic_central5) 
	egen DC34 = rowmean(diastolic_central3 diastolic_central4) 
	egen DC23 = rowmean(diastolic_central2 diastolic_central3) 
	egen DC12 = rowmean(diastolic_central1 diastolic_central2) 
	gen DC1 =  diastolic_central1 

	* replace with correct rowmean value based on how many BP measurements taken:
	gen avg_DC = DC45
	replace avg_DC = DC34 if diastolic_central5 == . 
	replace avg_DC = DC23 if diastolic_central5 == . & diastolic_central4 == . 
	replace avg_DC = DC12 if diastolic_central5 == . & diastolic_central4 == . & diastolic_central3 == .  	
	replace avg_DC = DC1  if diastolic_central5 == . & diastolic_central4 == . & diastolic_central3 == . & diastolic_central2 == . 

////////////////////////////////////////////////////////////////////////////////

	*	Check against existing vars:
	order avg_SB, after(sys_brachial)
	order avg_SC, after(sys_central)
	order avg_DB, after(dia_brachial)
	order avg_DC, after(dia_central)

	*	Drop old vars and rename new:
	drop sys_central dia_central sys_brachial dia_brachial
	rename avg_SC sys_central
	rename avg_DC dia_central
	rename avg_SB sys_brachial
	rename avg_DB dia_brachial

**	Recalculate min/max BP measurements (Leona's code):

	gen min_db_new=min(diastolic_brachial1, diastolic_brachial2, diastolic_brachial3, diastolic_brachial4, diastolic_brachial5, diastolic_brachial6)
	gen min_dc_new=min(diastolic_central1, diastolic_central2, diastolic_central3, diastolic_central4, diastolic_central5, diastolic_central6)
	gen min_sb_new=min(systolic_brachial1, systolic_brachial2, systolic_brachial3, systolic_brachial4, systolic_brachial5, systolic_brachial6)
	gen min_sc_new=min(systolic_central1, systolic_central2, systolic_central3, systolic_central4, systolic_central5, systolic_central6)
	gen max_db_new=max(diastolic_brachial1, diastolic_brachial2, diastolic_brachial3, diastolic_brachial4, diastolic_brachial5, diastolic_brachial6)
	gen max_dc_new=max(diastolic_central1, diastolic_central2, diastolic_central3, diastolic_central4, diastolic_central5, diastolic_central6)
	gen max_sb_new=max(systolic_brachial1, systolic_brachial2, systolic_brachial3, systolic_brachial4, systolic_brachial5, systolic_brachial6)
	gen max_sc_new=max(systolic_central1, systolic_central2, systolic_central3, systolic_central4, systolic_central5, systolic_central6)

	*	Check against existing vars:
	order min_db_new, after(min_db)
	order min_dc_new, after(min_dc)
	order min_sb_new, after(min_sb)
	order min_sc_new, after(min_sc)
	order max_db_new, after(max_db)
	order max_dc_new, after(max_dc)
	order max_sb_new, after(max_sb)
	order max_sc_new, after(max_sc)

	*	Drop old vars and rename new:
	drop min_db min_dc min_sb min_sc max_db max_dc max_sb max_sc
	rename min*_new min*
	rename max*_new max*
	drop SB45-DC1
	
////////////////////////////////////////////////////////////////////////////////

**	Second round of edits (28Mar2023):

** 	Replace zero value for PT2018103900 in S1 BP data as per Talia's 24Mar23 email (use average of first 2 measurements):
	replace systolic_central3=((systolic_central1+systolic_central2)/2) if systolic_central3==0 & ptc_id=="PT2018103900"
	replace diastolic_central3=((diastolic_central1+diastolic_central2)/2) if diastolic_central3==0 & ptc_id=="PT2018103900"

**	Generate composite ban status variable:
	gen ban_status_composite=.
	replace ban_status_composite=1 if Ban_2019==1
	replace ban_status_composite=2 if Ban_2020==1
	replace ban_status_composite=3 if Ban_2021==1
	replace ban_status_composite=0 if Not_yet==1
	la define L_ban_status 0"No ban" 1"2019" 2"2020" 3"2021"
	la values ban_status_composite L_ban_status
	*** N.BBan status of 西王平 straddles 2019 and 2020 so 70 ptc are listed as having entered the ban in 2019 AND 2020XY/Talia's view is that this village shouldn't be considered treated until the later year we have indicated (so ban_status_composite==2020 for those ptc).
	order ban_status_composite, after(Not_yet)
	drop Ban_2018
	rename Ban_2019 ban_status_2019
	rename Ban_2020 ban_status_2020
	rename Ban_2021 ban_status_2021
	rename Not_yet ban_status_no
	
**	Generate comments var for each ptc:
	gen comments=""

**	Generate average grip strength vars for S4:
	egen avg_gripR4 = rowmean(grip_right_1 grip_right_2 grip_right_3) if !missing(grip_right_1, grip_right_2, grip_right_3)
	egen avg_gripL4 = rowmean(grip_left_1 grip_left_2 grip_left_3) if !missing(grip_left_1, grip_left_2, grip_left_3)
	replace avg_gripR=avg_gripR4
	replace avg_gripL=avg_gripL4

**	Generate rounded average grip strength vars:
	gen avg_gripRt4=.
	gen avg_gripLt4=.
	replace avg_gripRt4=round(avg_gripR4)
	replace avg_gripLt4=round(avg_gripL4)
	order avg_gripRt4, after(avg_gripRt) // compare against existing var
	order avg_gripLt4, after(avg_gripLt) // compare against existing var
	tostring avg_gripRt4, replace 
	tostring avg_gripLt4, replace 
	replace avg_gripRt=avg_gripRt4
	replace avg_gripLt=avg_gripLt4
	drop avg_gripRt4 avg_gripLt4 avg_gripR4 avg_gripL4

////////////////////////////////////////////////////////////////////////////////

///	Merge in PE data:

**	Merge and harmonize variables:
	merge m:1 ptc_id wave using "/Users/leona/Desktop/BHET_PE_S4_29Mar23.dta"
	replace DutyCycle_exposure=1 if wave==4 // per XY's 28Mar23 email 
	replace Volumem3="0.19024" if Volumem3=="0,19024" // in order to convert to numeric
	destring Volumem3, replace // convert to numeric
	replace Volumem3_exposure=Volumem3 if wave==4 // populate this var for S4
	replace AfterFlowrateLmin_exposure=Aftersampling if wave==4
	drop if _merge==2
	
**	Clean and destring:
	replace AfterFlowrateLmin_exposure="1.8" if AfterFlowrateLmin_exposure=="1，8"
	replace AfterFlowrateLmin_exposure="0.909" if AfterFlowrateLmin_exposure=="0,909"
	destring AfterFlowrateLmin_exposure, replace

**	Add LF and PE notes to comments var:
	replace comments="LF: " + lung_comments if lung_comments!=""
	replace comments=comments + "\ PE: "+ Notes if Notes!=""
	order comments, before(master_version)
	drop FilterID-_merge lung_comments
	

////////////////////////////////////////////////////////////////////////////////

	
**	MASTER DATA GENERAL HOUSEKEEPING (30May2023):

**	After checking distributions for all master data variables, I found outlying/weird values from previous seasons that I either resolved with Talia, Xiang, Jill, and Xiaoying or noted in the "comments" variable. While I was at it, I did some general housekeeping with regard to variable names, value labels, formatting, etc. 

*	Drop vars as per Talia's May 23 email:
	drop common* bedroom* kitchen* bathroom* storage* other*


*	Fix value labels:
	label define occupation 1 "Agriculture and related workers" 2 "Factory worker" 3 "Government worker" 4 "Professional/technical workers" 5 "Sales and service workers" 6 "Retired" 7 "House wife/husband" 8 "Self-employed" 9 "Unemployed" 10 "Other or not stated" 11 "Mining" 12 "Construction", replace

	label define freq_drink 1 "Never" 2 "1 or 2 times in the past year" 3 "3 to 11 times in the past year" 4 "Once a month" 5 "2 to 3 times a month" 6 "1-2 times per week" 7 "3-4 times per week" 8 "5-6 times per week" 9 "Every day", replace

	label define freq_exercising 1 "None" 2 "Occasionally (<once per week)" 3 "1-2 days per week" 4 "3-5 days per week" 5 "Everyday", replace

	label define freq_farming 1 "None" 2 "Occasionally (<1 once per week)" 3 "1-2 days per week" 4 "3-5 days per week" 5 "Everyday", replace

	label define snore 1 "No" 2 "Occasionally (<once per week)" 3 "1-2 days per week" 4 "3-5 days per week" 5 "Daily or almost every day" 6 "Unknown", replace

	label define noise_HomePump 1 "Not at all" 2 "Slightly" 3 "Moderately" 4 "Very" 5 "Extremely" 6 "Not applicable", replace

	label define noise_NeighborPump 1 "Not at all" 2 "Slightly" 3 "Moderately" 4 "Very" 5 "Extremely" 6 "Not applicable", replace

	label define cuff_size 1 "Small (17-26cm)" 2 "Medium (24-32cm)" 3 "Large (>32cm)", replace

	la define gender_health 1"Male" 2"Female"
	la values gender_health gender_health

	label define heated_time1 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define heated_time10 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours"4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define heated_time11 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours"4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define heated_time12 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours"4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define heated_time13 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define heated_time2 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define heated_time3 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define heated_time4 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define heated_time5 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define heated_time6 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define heated_time7 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define heated_time8 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define heated_time9 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define occupied_time1 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define occupied_time10 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define occupied_time11 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define occupied_time12 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define occupied_time13 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define occupied_time2 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define occupied_time3 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define occupied_time4 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define occupied_time5 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define occupied_time6 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define occupied_time7 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define occupied_time8 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	label define occupied_time9 1 "Seldom or only on special occasions" 2 "<=2 hours" 3 "2-4 hours" 4 "4-6 hours" 5 "6-8 hours" 6 "8-10 hours" 7 "10-12 hours" 8 "12-16 hours" 9 "all day and night", replace

	la define ID_survey_type 1 "Heat survey" 2 "Health survey" 6 "Nurse survey"
	la values ID_survey_type ID_survey_type
	
	la define NAs 1"Populated" 2"Unknown by participant" 3"Missing"

	la define grip_quality 1"Participant could not complete grip strength measurement" 2"Participant could not complete grip strength measurement (judged by participant)" 3"Participant could not complete grip strength measurement (judged by enumerator)"

	foreach v of varlist *quality{
		replace `v'="1" if `v'=="#"
		replace `v'="2" if `v'=="#n"
		replace `v'="2" if `v'=="#N"
		replace `v'="2" if `v'=="n"
		replace `v'="2" if `v'=="N"
		replace `v'="3" if `v'=="#p"
		replace `v'="3" if `v'=="#P"
		replace `v'="3" if `v'=="p"
		replace `v'="3" if `v'=="P"
		destring, replace
		la values `v' grip_quality
	}

*	Corrections as per May 24 meeting with Talia:
	replace FeNO="<5" if FeNO=="Ôºú5"

*	Corrections as per May 24 email exchange with Jill:
	replace hours_sleep=. if hours_sleep<1
	replace AGE_ASTHMA=. if AGE_ASTHMA==0
	replace AGE_CHRONIC_HEPATITIS=. if AGE_CHRONIC_HEPATITIS==0
	replace arm_circum_r=. if arm_circum_r<10
	replace arm_circum_r=. if arm_circum_r==275
	replace waist_circ=. if waist_circ==21.5
	replace grip_right_2=17 if grip_right_2==1817
	replace grip_right_3=38 if grip_right_3==380
	replace grip_left_3=32 if grip_left_3==320

*	Create and populate NA vars:
	gen hours_sleep_NAs=.
	replace hours_sleep_NAs=1 if hours_sleep!=.
	replace hours_sleep_NAs=2 if hours_sleep==99
	replace hours_sleep_NAs=3 if hours_sleep==.
	replace hours_sleep=. if hours_sleep==99
	order hours_sleep_NAs, after(hours_sleep)
	
	gen waist_circ_NAs=.
	replace waist_circ_NAs=1 if waist_circ!=.
	replace waist_circ_NAs=2 if waist_circ==99
	replace waist_circ_NAs=3 if waist_circ==.
	replace waist_circ=. if waist_circ==99
	order waist_circ_NAs, after(waist_circ)

	gen weight_NAs=.
	replace weight_NAs=1 if weight!=.
	replace weight_NAs=2 if weight==99
	replace weight_NAs=3 if weight==.
	replace weight=. if weight==99
	order weight_NAs, after(weight)

	gen Quantity_Briquettes_NAs=.
	replace Quantity_Briquettes_NAs=1 if Quantity_Briquettes!=.
	replace Quantity_Briquettes_NAs=2 if Quantity_Briquettes==999
	replace Quantity_Briquettes_NAs=3 if Quantity_Briquettes==.
	replace Quantity_Briquettes=. if Quantity_Briquettes==999
	order Quantity_Briquettes_NAs, after(Quantity_Briquettes)
	
	replace Price_Briquettes_NAs=2 if Price_Briquettes==999
	replace Price_Briquettes=. if Price_Briquettes==999

	gen Price_Honeycomb_coal_NAs=.
	replace Price_Honeycomb_coal_NAs=1 if Price_Honeycomb_coal!=.
	replace Price_Honeycomb_coal_NAs=2 if Price_Honeycomb_coal==999
	replace Price_Honeycomb_coal_NAs=2 if Price_Honeycomb_coal==800
	replace Price_Honeycomb_coal_NAs=3 if Price_Honeycomb_coal==.
	replace Price_Honeycomb_coal=. if Price_Honeycomb_coal==999
	order Price_Honeycomb_coal_NAs, after(Price_Honeycomb_coal)
	
	replace Price_Coal_NAs=2 if Price_coal==999
	replace Price_coal=. if Price_coal==999
	rename Price_coal Price_Coal

	replace Quantity_LPG_NAs=2 if Quantity_LPG==999
	replace Quantity_LPG=. if Quantity_LPG==999

	replace Price_LPG_NAs=2 if Price_LPG==999
	replace Price_LPG=. if Price_LPG==999

	gen Quantity_ELE_NAs=.
	replace Quantity_ELE_NAs=1 if Quantity_ELE!=.
	replace Quantity_ELE_NAs=2 if Quantity_ELE==999
	replace Quantity_ELE_NAs=3 if Quantity_ELE==.
	replace Quantity_ELE=. if Quantity_ELE==999
	order Quantity_ELE_NAs, after(Quantity_ELE)
	
	gen Price_ELE_NAs=.
	replace Price_ELE_NAs=1 if Price_ELE!=.
	replace Price_ELE_NAs=2 if Price_ELE==999
	replace Price_ELE_NAs=3 if Price_ELE==.
	replace Price_ELE=. if Price_ELE==999
	order Price_ELE_NAs, after(Price_ELE)

	replace Quantity_Wood_NAs=2 if Quantity_Wood==999
	replace Quantity_Wood=. if Quantity_Wood==999

	replace Price_Wood_NAs=2 if Price_Wood==999
	replace Price_Wood_NAs=2 if Price_Wood==9999
	replace Price_Wood=. if Price_Wood==999
	replace Price_Wood=. if Price_Wood==9999

	gen Price_Charcoal_NAs=.
	replace Price_Charcoal_NAs=1 if Price_Charcoal!=.
	replace Price_Charcoal_NAs=2 if Price_Charcoal==999
	replace Price_Charcoal_NAs=3 if Price_Charcoal==.
	replace Price_Charcoal=. if Price_Charcoal==999
	order Price_Charcoal_NAs, after(Price_Charcoal)
	
	gen Quantity_Other_NAs=.
	replace Quantity_Other_NAs=1 if Quantity_Other!=.
	replace Quantity_Other_NAs=2 if Quantity_Other==999
	replace Quantity_Other_NAs=3 if Quantity_Other==.
	replace Quantity_Other=. if Quantity_Other==999
	order Quantity_Other_NAs, after(Quantity_Other)

	gen Price_Other_NAs=.
	replace Price_Other_NAs=1 if Price_Other!=.
	replace Price_Other_NAs=2 if Price_Other==999
	replace Price_Other_NAs=3 if Price_Other==.
	replace Price_Other=. if Price_Other==999
	order Price_Other_NAs, after(Price_Other)

	la values *NAs NAs
	
*	Back-populate "measures" variable for S1 and S2:
	replace measures=6 if wave==1 & systolic_brachial1!=. & systolic_brachial2!=. & systolic_brachial3!=. & systolic_brachial4!=. & systolic_brachial5!=. & systolic_brachial6!=.
	replace measures=5 if wave==1 & systolic_brachial1!=. & systolic_brachial2!=. & systolic_brachial3!=. & systolic_brachial4!=. & systolic_brachial5!=. & systolic_brachial6==.
	replace measures=4 if wave==1 & systolic_brachial1!=. & systolic_brachial2!=. & systolic_brachial3!=. & systolic_brachial4!=. & systolic_brachial5==. & systolic_brachial6==.
	replace measures=3 if wave==1 & systolic_brachial1!=. & systolic_brachial2!=. & systolic_brachial3!=. & systolic_brachial4==. & systolic_brachial5==. & systolic_brachial6==.
	replace measures=2 if wave==1 & systolic_brachial1!=. & systolic_brachial2!=. & systolic_brachial3==. & systolic_brachial4==. & systolic_brachial5==. & systolic_brachial6==.
	replace measures=1 if wave==1 & systolic_brachial1!=. & systolic_brachial2==. & systolic_brachial3==. & systolic_brachial4==. & systolic_brachial5==. & systolic_brachial6==.
	replace measures=0 if wave==1 & systolic_brachial1==. & systolic_brachial2==. & systolic_brachial3==. & systolic_brachial4==. & systolic_brachial5==. & systolic_brachial6==.
	replace measures=6 if wave==2 & systolic_brachial1!=. & systolic_brachial2!=. & systolic_brachial3!=. & systolic_brachial4!=. & systolic_brachial5!=. & systolic_brachial6!=.
	replace measures=5 if wave==2 & systolic_brachial1!=. & systolic_brachial2!=. & systolic_brachial3!=. & systolic_brachial4!=. & systolic_brachial5!=. & systolic_brachial6==.
	replace measures=4 if wave==2 & systolic_brachial1!=. & systolic_brachial2!=. & systolic_brachial3!=. & systolic_brachial4!=. & systolic_brachial5==. & systolic_brachial6==.
	replace measures=3 if wave==2 & systolic_brachial1!=. & systolic_brachial2!=. & systolic_brachial3!=. & systolic_brachial4==. & systolic_brachial5==. & systolic_brachial6==.
	replace measures=2 if wave==2 & systolic_brachial1!=. & systolic_brachial2!=. & systolic_brachial3==. & systolic_brachial4==. & systolic_brachial5==. & systolic_brachial6==.
	replace measures=1 if wave==2 & systolic_brachial1!=. & systolic_brachial2==. & systolic_brachial3==. & systolic_brachial4==. & systolic_brachial5==. & systolic_brachial6==.
	replace measures=0 if wave==2 & systolic_brachial1==. & systolic_brachial2==. & systolic_brachial3==. & systolic_brachial4==. & systolic_brachial5==. & systolic_brachial6==.

*	Consolidate/reformat ELE* variables:
	rename ele_pay_photo ELE_pay_jpeg
	replace ELE_pay_photo=1 if ELE_pay_jpeg!=""	// populate binary photo var with S4
	replace ELE_id_photo_jpeg=ele_ID_photo if ele_ID_photo!=""
	replace ELE_ID_photo=1 if ELE_id_photo_jpeg!=""
	rename ELE_id_photo_jpeg ELE_ID_jpeg
	order ELE_ID_photo, before(ELE_ID_jpeg)
	//
	replace ELE_remaining_money_photo_jpeg=ele_remaining_money_photo if ele_remaining_money_photo!=""
	replace ELE_remaining_money_photo=1 if ELE_remaining_money_photo_jpeg!=""
	rename ELE_remaining_money_photo_jpeg ELE_remaining_money_jpeg
	order ELE_remaining_money_photo, before(ELE_remaining_money_jpeg)
	//
	replace ELE_last_month_photo_jpeg=ele_last_month_photo if ele_last_month_photo!=""
	replace ELE_last_month_photo=1 if ELE_last_month_photo_jpeg!=""
	rename ELE_last_month_photo_jpeg ELE_last_month_jpeg
	order ELE_last_month_photo, before(ELE_last_month_jpeg)
	//
	replace ELE_last_twomonths_photo_jpeg=ele_last_twomonths_photo if ele_last_twomonths_photo!=""
	replace ELE_last_twomonths_photo=1 if ELE_last_twomonths_photo_jpeg!=""
	rename ELE_last_twomonths_photo_jpeg ELE_last_twomonths_jpeg
	order ELE_last_twomonths_photo, before(ELE_last_twomonths_jpeg)
	//
	replace ele_id_photo_2_jpeg=ele_ID_photo_2 if ele_ID_photo_2!=""
	replace ele_id_photo_2=1 if ele_id_photo_2_jpeg!=""
	rename ele_id_photo_2_jpeg ELE_ID_2_jpeg
	replace ELE_ID_photo_2=ele_id_photo_2 if wave==4
	order ELE_ID_photo_2, before(ELE_ID_2_jpeg)
	rename ELE_ID_photo_2 ELE_ID_2_photo
	//
	replace ELE_pay_photo=ele_bill_photo if wave==1	// same var - populate the S4 version with S1
	drop ele_ID_photo_2 ele_last_twomonths_photo ele_last_month_photo ele_remaining_money_photo ele_ID_photo ele_id_photo_2 ele_bill_photo
	rename ele* ELE*
	rename ELE_remaining_money_photo_2 ELE_remaining_money_2_photo
	rename ELE_remaining_money_photo_2_jpeg ELE_remaining_money_2_jpeg
	rename ELE_last_month_photo_2 ELE_last_month_2_photo
	rename ELE_last_month_photo_2_jpeg ELE_last_month_2_jpeg
	rename ELE_last_twomonths_photo_2 ELE_last_twomonths_2_photo
	rename ELE_last_twomonths_photo_2_jpeg ELE_last_twomonths_2_jpeg

*	Streamline these var names: 
	rename ZEFLUOR_photo zefluor_photo
	rename ZEFLUOR_photo_jpeg zefluor_photo_jpeg
	rename QUARTZ_photo quartz_photo
	rename QUARTZ_photo_jpeg quartz_photo_jpeg
	rename marrital_stat_health marital_stat_health //misspelling

*	ELE_meterIDs and ELE_meterID appear to be the same var except ELE_meterIDs only exists for S2:
	replace ELE_meterID=ELE_meterIDs if ELE_meterID==. & ELE_meterIDs!=.
	drop ELE_meterIDs

*	As per June 1 email with Talia:
	drop TA_VIALS
	drop numBLD
	order Other_Fuel, before(FUEL_Other)
	la var Other_Fuel "Name of the other fuel used"
	la var FUEL_Other "Whether 'other' fuel type is used for cooking, heating, or heating water in winter"

*	As per June 1 email from Xiaoying:
	replace zefluor_id="NA" if zefluor_id=="9999999"
	replace zefluor_id="NA" if zefluor_id=="99999999"
	replace DurationHHMM_exposure="21:03" if DurationHHMM_exposure=="21.049"
	replace DurationHHMM_exposure="22:56" if DurationHHMM_exposure=="22.932"

*	As per June 6 meeting with Xiang:
	replace comments="cs_num_cigs: likely very occasional smoker (ptc smoked 1 per day in S4)" if cs_num_cigs==0
	replace comments="cs_age_smoking: very likely correct" if cs_num_cigs<8
	replace Quantity_Briquettes=5.5 if Quantity_Briquettes==550
	replace Quantity_Honeycomb_coal=800 if ptc_id=="PT2018113000" & wave==2
	replace Price_Honeycomb_coal=3 if ptc_id=="PT2018113000" & wave==2
	replace comments="Price_Coal==0: likely because not supposed to use coal but had some leftover from previous years" if Price_Coal==0
	replace Price_Biogas=2 if Price_Biogas==100
	replace comments="exp_heating_eqp: unsure as to why this value is so outsized" if exp_heating_eqp==231250 
	replace comments="hh_guests==750: they opened a hotel hosting 750 guests per month" if hh_guests==750 
	drop asset_complete
	// subtracted year they moved in from survey year to calculate how long they've lived there: 
	replace hh_residtime=61 if hh_residtime==1958
	replace hh_residtime=42 if hh_residtime==1980
	replace hh_residtime=6 if hh_residtime==2015
	// subtracted "years ago" value from survey year to calculate year house was built:
	replace house_builtyear=1978 if house_builtyear==43
	replace house_builtyear=1937 if house_builtyear==84
	replace house_builtyear=1934 if house_builtyear==87

////////////////////////////////////////////////////////////////////////////////

**	Xiaoying's corrections (emailed to Leona Jun 30, 2023)
	//use "/Users/leona/Desktop/BHET_master_data_7Jun2023.dta", clear
	replace occupation=1 if hh_id=="HH20180703" & wave==2
	replace occupation=1 if hh_id=="HH20180792" & wave==2
	replace occupation=7 if hh_id=="HH20180800" & wave==2
	replace occupation=1 if hh_id=="HH20180926" & wave==2
	replace occupation=1 if hh_id=="HH20181009" & wave==2
	replace marital_stat_health=1 if hh_id=="HH20180703" & wave==2
	replace marital_stat_health=4 if hh_id=="HH20180792" & wave==2
	replace marital_stat_health=1 if hh_id=="HH20180926" & wave==2
	replace marital_stat_health=1 if hh_id=="HH20181009" & wave==2
	replace house_builtyear=2013 if hh_id=="HH20190272" & wave==4
	replace house_builtyear_NAs=2 if hh_id=="HH20210353" & wave==4
	replace house_builtyear=2018 if hh_id=="HH20190272" & wave==4
	
////////////////////////////////////////////////////////////////////////////////

**	Harmonize NA vs. missing for heat_gpsStart and heat_gpsEnd (as per Talia's July 25 email)
	replace heat_gpsStart="" if heat_gpsStart=="NA,NA,NA"
	replace heat_gpsEnd="" if heat_gpsEnd=="NA,NA,NA"

////////////////////////////////////////////////////////////////////////////////

**	Generate master version variable:
	drop master_version
	gen master_version="26Jul2023 16:51:00"
	sort wave hh_id ptc_id
	save "/Users/leona/Desktop/BHET_master_data_26Jul2023.dta", replace	
	
////////////////////////////////////////////////////////////////////////////////
//////////////////////BACK-CORRECT PARTICIPANT DOB AND AGE//////////////////////
////////////////////////////////////////////////////////////////////////////////

//	(N.B.: These corrections are based on email correspondence between Talia and Leona and BHET team meetings in July, 2023)

**	Export datasheet containing ONLY S2 observations (ptc_id, survey date, DOB) and remerge with master data to replace S1 DOB with S2 DOB to harmonize DOB for participants that don't exist in S4:

	use "/Users/leona/Desktop/BHET_master_data_26Jul2023.dta", clear
	keep if wave==2
	keep hh_id ptc_id wave date dob_health
	rename dob_health dob_health_S2
	save "/Users/leona/Desktop/DOB_S2.dta"

**	Export datasheet containing ONLY S4 observations (ptc_id, survey date, DOB) and remerge with master data to replace S1 & 2 DOB with S4 to harmonize DOB for all remaining participants:

	use "/Users/leona/Desktop/BHET_master_data_26Jul2023.dta", clear
	keep if wave==4
	keep hh_id ptc_id wave date dob_health
	rename dob_health dob_health_S4
	save "/Users/leona/Desktop/DOB_S4.dta"

////////////////////////////////////////////////////////////////////////////////

**	Back-correct DOB:
	
	use "/Users/leona/Desktop/BHET_master_data_26Jul2023.dta", clear

*	Merge in exported S2 DOB:
	merge m:m ptc_id using "/Users/leona/Desktop/DOB_S2.dta", keepusing(dob_health_S2)
	drop _merge
	
*	Create new DOB var to house back-corrected DOB:
	gen dob_health_NEW=dob_health
	format %tc dob_health_NEW
	
*	Back-correct S1 DOB with S2 DOB:
	replace dob_health_NEW=dob_health_S2 if dob_health_S2!=.

*	Merge in exported S4 DOB:
	merge m:m ptc_id using "/Users/leona/Desktop/DOB_S4.dta", keepusing(dob_health_S4)
	drop _merge

*	Back-correct S1 and S2 DOB with S4 DOB:
	replace dob_health_NEW=dob_health_S4 if dob_health_S4!=.

////////////////////////////////////////////////////////////////////////////////

**	Convert date vars to recalculate age:

*	Convert survey datetime to date only:
	gen date2 = dofc(date)
	format date2 %td
	rename date datetime
	rename date2 date

*	Convert DOB to date only:
	gen dob_NEW = dofc(dob_health_NEW)
	format dob_NEW %td

*	Convert survey date to year only:
	gen survey_year=yofd(date)

*	Convert DOB to year only:
	gen dob_year_NEW=yofd(dob_NEW)

*	Create new age var based on these date vars:
	gen age_CORRECTED=survey_year-dob_year_NEW

*	Reorder and drop superfluous vars:
	rename dob_NEW dob_CORRECTED
	order dob_CORRECTED, after(dob_health)
	order age_CORRECTED, after(age_health)
	drop dob_health_S2-dob_year_NEW
	sort wave hh_id ptc_id

////////////////////////////////////////////////////////////////////////////////

**	Generate master version variable:
	drop master_version
	gen master_version="26Jul2023 16:54:00"
	sort wave hh_id ptc_id
	export delimited using "/Users/leona/Desktop/BHET_master_data_26Jul2023.csv", nolabel replace
	save "/Users/leona/Desktop/BHET_master_data_26Jul2023.dta", replace
	
////////////////////////////////////////////////////////////////////////////////

**	August 2023 revisions:

	use "/Users/leona/Desktop/BHET_master_data_26Jul2023.dta", clear
	
**	Drop S2 tester HH as per Xiang's Aug 2 email and Jill's Aug 19 response:
	drop if hh_id=="HH20190300" & wave==2
	drop if hh_id=="HH20190268" & wave==2
	drop if hh_id=="HH20190116" & wave==2
	drop if hh_id=="HH20190105" & wave==2
	drop if hh_id=="HH20190108" & wave==2
	drop if hh_id=="HH20190156" & wave==2
	drop if hh_id=="HH20190273" & wave==2
	drop if hh_id=="HH20190140" & wave==2
	drop if hh_id=="HH20190257" & wave==2

** Revised ptc gender corrections (emailed to Leona Sep 1, 2023):
	replace ptc_id="PT2018068101" if ptc_id=="PT2018068100" & wave==4
	replace ptc_id="PT2018096701" if ptc_id=="PT2018096700" & wave==4
	replace gender_health=2 if ptc_id=="PT2018096701" & wave==2
	replace ptc_id="PT2018094401" if ptc_id=="PT2018094400" & wave==4
	replace ptc_id="PT2018026701" if ptc_id=="PT2018026700" & wave==4
	replace ptc_id="PT2018029901" if ptc_id=="PT2018029900" & wave==4
	replace ptc_id="PT2018042701" if ptc_id=="PT2018042700" & wave==4
	replace ptc_id="PT2018052701" if ptc_id=="PT2018052700" & wave==2	
	replace ptc_id="PT2018052701" if ptc_id=="PT2018052700" & wave==4
	replace gender_health=2 if ptc_id=="PT2018059401" & wave==2
	replace gender_health=1 if ptc_id=="PT2018066301" & wave==2
	replace ptc_id="PT2018087201" if ptc_id=="PT2018087200" & wave==4
	replace ptc_id="PT2018088800" if ptc_id=="PT2018088801" & wave==4
	replace gender_health=2 if ptc_id=="PT2018089900" & wave==4
	replace ptc_id="PT2018094401" if ptc_id=="PT2018094400" & wave==4
	replace ptc_id="PT2018102401" if ptc_id=="PT2018102400" & wave==4
	replace gender_health=2 if ptc_id=="PT2018106801" & wave==2
	replace ptc_id="PT2018106801" if ptc_id=="PT2018106800" & wave==4

**	Generate master version variable:
	drop master_version
	gen master_version="6Sep2023 11:24:00"
	sort wave hh_id ptc_id
	export delimited using "/Users/leona/Desktop/BHET_master_data_6Sep2023.csv", nolabel replace
	save "/Users/leona/Desktop/BHET_master_data_6Sep2023.dta", replace
	
////////////////////////////////////////////////////////////////////////////////
