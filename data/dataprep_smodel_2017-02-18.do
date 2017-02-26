clear all
set more off
capture log close
set seed 1234
//set matsize 11000

cd "/home/lohia/projects/sex_sel/data/"


global date 151116

********************************
* SET NUMBER OF INCOME CLASSES *
********************************
* choose from 6, 8 or 10

global numincls = 10
********************************

*****************
* SET AGE GROUP *
*****************
* choose from "06" or "713"

global agegrp = "713"

******************

log using "Logs/treg_$date", replace

* Prepare data for bootstrap - running another do file
//do "Do files/BS_dataprep_0311"
/*
use "Census data/DataForReg.dta", clear

global PredictedIncome pred_income
//global ns = 5

keep if curr_hist_match==3 // villages for which we have both the historical and the current information

	gen d_06 = 0
	replace d_06 = 1 if d_girl6!=.
	bys HHid: egen HHunder6 = max(d_06)

***************
* First stage *
***************

	duplicates drop HHid, force

	* (1) Using hist_wealth1 "DD on account of land revenue per cultivated acre"
	xi: reg monthlyincome i.caste*hist_wealth1, cluster(A2)
	testparm hist_wealth1 _IcasXhist* // F-stat (298K) 20.40
	predict pred_income

	* count HH in each caste
	gen c1 = 1 if fs==2|fs==9
	egen numHH_caste_fs29 = count(c1), by(caste)
	lab var numHH_caste_fs29 "Number of fs2/9 HH in each caste"

	* count HH with atleast one under six and fs=2| fs=9
	gen c2 = 1 if HHunder6==1 & (fs==2|fs==9)
	egen numHH_caste_fs29_u6 = count(c2), by(caste)
	lab var numHH_caste_fs29_u6 "Number of fs2/9 HH in each caste with atleast one child under 6"

	keep HHid pred_* numHH_caste*

	merge 1:m HHid using "Census data/DataForReg.dta", keep(3) nogen

	//save "Temporary files/temp_2504", replace

	//use "Temporary files/temp_2504", clear

	* II * rank based on predicted income per OWN youngest generation/ OWN child +2 (not directly estimated)
			gen pred_inc_younggen2p =.
			replace pred_inc_younggen2p = $PredictedIncome /(total_yg_own+2) if LivesWithFather==1 & LivesWithMother==1 & d_youngestgen==1
			lab var pred_inc_younggen2p "Predicted HH montlyincome divided by no. in the youngest generation+ 2(parents)"

			gen log_pred_inc_younggen2p = log(pred_inc_younggen2p)

			sort caste log_pred_inc_younggen2p
			bys caste log_pred_inc_younggen2p : gen newid = 1 if _n==1
			by caste: gen rank_pred_inc_yg2p = sum(newid) if pred_inc_younggen2p!=.

			bys caste: egen num_ranks = max(rank_pred_inc_yg2p)

			gen ref_rank_pred_inc_yg2p = rank_pred_inc_yg2p/num_ranks
			lab var ref_rank_pred_inc_yg2p "Rank in within caste income distribution"

			drop newid
			drop num_ranks

	* X * rank based on actual/reported income divided by family size
			gen inc_famsize= monthlyincome/(total_yg_own+2) if LivesWithFather==1 & LivesWithMother==1 & d_youngestgen==1

			gen log_inc_famsize = log(inc_famsize)

			sort caste log_inc_famsize
			bys caste log_inc_famsize : gen newid = 1 if _n==1
			by caste: gen rank_inc_famsize = sum(newid) if inc_famsize!=.

			bys caste: egen num_ranks = max(rank_inc_famsize)

			gen ref_rank_inc_famsize = rank_inc_famsize/num_ranks
			lab var ref_rank_inc_famsize "Rank in within caste income distribution"

			drop newid
			drop num_ranks

*Prepare data for second stage
	keep if fs==2 | fs==9
	// lab def fs 2 "Pe and family", add
	// lab def fs 9 "Pe and siblings living with their parents", add

lab var pred_income "Predicted income"


* 2 * individual castes
keep if numHH_caste_fs29_u6>=500


* number of 0-6 year olds in a caste: to be reported as part of the graphs
		gen d_06 = 0
		replace d_06 = 1 if d_girl6!=.
		bys caste: egen num_06= total(d_06)


		* 7 to 13 year old girls
		gen d_girl713=.
		replace d_girl713 = 0 if C2_4>=7 & C2_4<=13
		replace d_girl713 = 1 if C2_4>=7 & C2_4<=13 & C2_3==2
		lab var d_girl713 "Dummy = 1 if girl between the age of 7-13 years, 0 for boys between the age of 7-13 years"

		gen d_713 = 0
		replace d_713 = 1 if d_girl713!=.
		bys caste: egen num_713= total(d_713)

save "Temporary files/readyforinc_cls", replace
*/

use "readyforinc_cls", clear

* use HH level data to calculate income classes by caste
duplicates drop HHid, force
global income pred_inc_younggen2p
gen inc_cls = .

* same as saying keep if numHH_caste_fs29_u6>=500
keep if caste==28| caste==38| caste==1| caste==14| caste==33 |caste==5|caste==10|caste==23|caste==7 ///
| caste==44 | caste==17 | caste==4



* Create income classes for castes
set seed 1234
		levelsof caste if numHH_caste_fs29_u6>=150, local(castenumber)
		foreach x of local castenumber{
			sort caste HHid $income
			xtile inc_cls_`x' = $income if caste==`x', nq($numincls)
			di "caste = `x'"
			//tab inc_cls_`x' if caste==`x'
			replace inc_cls = inc_cls_`x' if caste==`x'
		}

* drop values with no income class assigned
drop if inc_cls==.
gen numcls=.

* check how many income classes created per caste
* stata forces identical observations to be in the same quantile
	levelsof caste if numHH_caste_fs29_u6>=100, local(castenumber)
	foreach x of local castenumber{

		bys caste inc_cls: replace numcls = 1 if _n==1 & caste==`x'

	}
	by caste: replace numcls = sum(numcls)
	bys caste: egen tot_numcls = max(numcls)
	tab caste tot_numcls

* keep only the caste and income class data for each HHid
keep HHid caste inc_cls
keep if inc_cls!=.
save "Temp/HHinc_cls", replace

* merge back to individual level data to collapse
use "readyforinc_cls", clear
merge m:1 HHid using "Temp/HHinc_cls"

keep if _m==3

* our benchamrk regressions
	local ranklist pred_inc_yg2p
	local reglist d_girl6
			foreach reg of local reglist{
				foreach ranktype of local ranklist{
					* income per child *
					reg `reg' ref_rank_`ranktype', robust cluster(A2)

					areg `reg' ref_rank_`ranktype', robust cluster(A2) absorb(caste)
				}
			}

* collapse to caste - income class level
collapse (mean) d_girl6 $income (sum)d_$agegrp , by(caste inc_cls)
drop if inc_cls==.

* regress fraction of girls on average income
reg d_girl6 $income [fweight=d_$agegrp]

* regress fraction of girls on income class dummies
reg d_girl6 i.inc_cls

	* same regression graphically
	sort inc_cls caste
	by inc_cls: egen mean_frac = mean(d_girl6)
	tw line mean_frac inc_cls

sort caste inc_cls
ren pred_inc_younggen2p inc
gen f = d_girl6 * d_713
gen m = d_713 - f
ren d_713 t
gen csr = m/f
gen X = (m+f)/2
gen Hk = (csr-1)/(csr+1)
order inc m f t csr Hk X caste inc_cls d_girl6 mean_frac
levelsof caste, local(levels) 
foreach l of local levels {
	preserve
	keep if caste == `l'
	drop caste inc_cls d_girl6 mean_frac
	outsheet using caste`l'.csv, comma
	restore
}
