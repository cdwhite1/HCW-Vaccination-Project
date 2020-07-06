

clear

global PATH "DIRECTORY FOR ANALYSIS CODE"

set matsize 2000







************************************************************************************************
************************************************************************************************
****************************** TABLES
************************************************************************************************
************************************************************************************************







************************************************************************************************
***************Table 2:   Main Result; heterogeneity across seasons/age/law Types
***************Table S2:  Differential Effects for heterogeneity results in Table 2
***************Table S5:  Duplicates Table 2 with wild bootstrapped standard errors
************************************************************************************************


****************Main Result (for Tables 2, S2)

*Import analysis file
cd "$PATH\Influenza National\"
clear
use flu_hcw_national_analysis

*Drop flu-years with match rate<50%
drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 

*Covariates: temperature, humidity, precipitation, population shares, unemployment rate, poverty rate, indicator for long term care laws, indicator for childcare laws, state-specific linear trends
local controls  temp* hum* prcp* popshare* unemp pov offer_ltc req_child i.state#c.ym
*Set panel variable (for fixed effects) to be state
xtset state

*regression for main result
xtreg pi_rate offer_post  i.ym `controls' [aweight=meanpop] , fe cluster(state)

*report regression result (not necessary, but easier to read)
lincom offer_post

*Get wild cluster bootstrapped p-value and CI for Table 5
boottest offer_post=0, cluster(state)

*Calculate baseline mean of the outcome (adopting states, years prior to adoption)
sum pi_rate if max_offer==1 & yearstooffer<0 [aweight=meanpop], meanonly




****************Heterogeneity across seasons (for Tables 2, S2, S5)

*See annotations from above for additional comments

cd "$PATH\Influenza National\"
clear
use flu_hcw_national_analysis

*Generate interaction for law and peak season
gen offer_season = offer_post*season

drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 

*"sf_group" is an identifier for stateXseason groups. With the interaction model, we include stateXseason FE and stateXseason trends
local controls  temp* hum* prcp* popshare* unemp pov offer_ltc req_child i.sf_group#c.ym
xtset sf_group

*Fill in results matrix with means (peak and non-peak season)
sum pi_rate if max_offer==1 & yearstooffer<0 & season==1 [aweight=meanpop]
sum pi_rate if max_offer==1 & yearstooffer<0 & season==0 [aweight=meanpop]

*Regression for season interaction model
xtreg pi_rate offer_post offer_season i.ym `controls' [aweight=meanpop] , fe cluster(state)

*Results: Effect of laws during non-peak season (Table 2), plus bootstrapped version for Table 5
lincom offer_post
boottest offer_post=0, cluster(state)

*Results: Effect of laws during peak season (Table 2), plus bootstrapped version for Table 5
lincom offer_post + offer_season
boottest offer_post + offer_season=0, cluster(state)

*Results: Differential effect for peak vs. non-peak (Table S2), plus bootstrapped version for Table 5
lincom offer_season
boottest offer_season=0, cluster(state)




****************Heterogeneity across age (for Tables 2, S2, S5)

*See annotations from Table 2 for additional comments

cd "$PATH\Influenza National\"
clear
use flu_hcw_national_analysis

*To run an interacted model with age, we need to reshape the data to long format so that each observation is state-year-month-age; the next 4 lines of code do this reshape.
keep pi_rate_age_065 pi_rate_age_65 max_offer yearstooffer season meanpop temp* hum* prcp* popshare* unemp pov offer_ltc req_child sf ym state flu_year offer_post
rename pi_rate_age_065 pi_rate1
rename pi_rate_age_65 pi_rate2
reshape long pi_rate, i(state ym) j(age)

drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 

*Generate interaction with 65+
gen o65 = (age==2)
gen offer_age = offer_post*o65

*Generate new variables for stateXage and timeXage FE
egen sa_group = group(state age)
egen yma_group = group(ym age)

*Note controls include the stateXage time trends
local controls  temp* hum* prcp* popshare* unemp pov offer_ltc req_child i.sa_group#c.ym
xtset sa_group

*Calculate means for each group
sum pi_rate if max_offer==1 & yearstooffer<0 & o65==1 [aweight=meanpop]
sum pi_rate if max_offer==1 & yearstooffer<0 & o65==0 [aweight=meanpop]

*Age heterogeneity regression
xtreg pi_rate  offer_post offer_age i.yma `controls' [aweight=meanpop] , fe cluster(state)

*Results: effect of laws for <65 (Table 2)
lincom offer_post
boottest offer_post=0, cluster(state)

*Results: effect of laws for >65 (Table 2)
lincom offer_post + offer_age
boottest offer_post + offer_age=0, cluster(state)

*Results: Differential effect of laws for >65 vs. <65 (Table S2)
lincom offer_age
boottest offer_age=0, cluster(state)






****************Heterogeneity Law Types (for Tables 2, S2, S5)

*See annotations from Table 2 for additional comments

cd "$PATH\Influenza National\"
clear
use flu_hcw_national_analysis

*Define treat1, treat2, treat3 to represent: no declination, declination, mask-required laws
gen treat1 = 0
replace treat1=1 if offer_post==1 
gen treat2=0
replace treat2=1 if offer_dec_post==1
gen treat3=0
replace treat3=1 if require_post==1

drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 

*Same controls/FE as the baseline specification, since heterogeneity is just over different states.
local controls  temp* hum* prcp* popshare* unemp pov offer_ltc req_child i.state#c.ym
xtset state

*regression for heterogeneity by law type
xtreg pi_rate treat1 treat2 treat3 i.ym `controls'  [aweight=meanpop], fe cluster(state)

*calculate means outcome for each group (pre-treatment mean for states that adopt each type of law)
bysort state: egen max_treat1=max(treat1)
bysort state: egen max_treat2=max(treat2)
bysort state: egen max_treat3=max(treat3)
sum pi_rate [aw=meanpop] if max_treat1==1 & max_treat2==0 & max_treat3==0 & offer_post==0
sum pi_rate [aw=meanpop] if max_treat2==1 & max_treat3==0 & offer_post==0
sum pi_rate [aw=meanpop] if max_treat3==1 & offer_post==0

*Results: effect for law without any requirements, i.e., non-declination laws (Table 2)
lincom treat1
boottest treat1=0, cluster(state)

*Results: effect for law with declination requirement but no mask requirement (Table 2)
lincom treat1+treat2
boottest treat1+treat2=0, cluster(state)

*Results: effect for law with decliation plus mask requirement (Table 2)
lincom treat1+treat2+treat3
boottest treat1+treat2+treat3=0, cluster(state)

*Results: differential effect for declination law vs non-declination law
lincom treat2
boottest treat2=0, cluster(state)

*Results: differential effect for mask law vs non-declination law
lincom treat2+treat3
boottest treat2+treat3=0, cluster(state)

*Results: differential effect for mask law vs declination law
lincom treat3
boottest treat3=0, cluster(state)









************************************************************************************************
***************Table 3: Sensitivity Tests
************************************************************************************************

************************Varying Fixed Effects and Controls (Panel A)

cd "$PATH\Influenza National\"

clear
use flu_hcw_national_analysis

*Two sets of optional controls (with and without state-specific linear trends)
local covariates temp* hum* prcp* popshare* unemp pov offer_ltc req_child 

drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014

*excludes covariates, trends
xtset state
xtreg pi_rate offer_post   i.ym  [aweight=meanpop], fe cluster(state)
*report results (not necessary, but easier to read)
lincom offer_post

*excludes trends
xtset state
xtreg pi_rate offer_post `covariates'  i.ym  [aweight=meanpop], fe cluster(state)
lincom offer_post

*excludes covariates
xtset state
xtreg pi_rate offer_post i.state#c.ym  i.ym  [aweight=meanpop], fe cluster(state)
lincom offer_post


************************Varying Years (Panel B)

cd "$PATH\Influenza National\"

clear
use flu_hcw_national_analysis

local controls temp* hum* prcp* popshare* unemp pov offer_ltc req_child i.state#c.ym

*Include bad match years
xtset state
xtreg pi_rate offer_post  `controls' i.ym  [aweight=meanpop], fe cluster(state)
lincom offer_post

*Drop bad match years  for remaining regressions
drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014

*Drop H1N1 years
xtreg pi_rate offer_post  `controls' i.ym if flu_year!=2008 & flu_year != 2009   [aweight=meanpop], fe cluster(state)
lincom offer_post



************************Poisson (Panel C)

cd "$PATH\Influenza National\"

clear
use flu_hcw_national_analysis

*drop bad match years (<50%) -- note that the mean match rate is 23.74% in these 4 seasons (1998+), and 83.19% in the other 14 seasons.
drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 
xtset state

local controls  temp* hum* prcp* popshare* unemp pov offer_ltc req_child i.state#c.ym

*Estimate Poisson FE version of the main specification
/*
xtpoisson pi offer_post  i.ym `controls' , fe exposure(pop_tot) vce(robust)
disp exp(_b[offer_post])-1
*/

*Predicted values don't work properly with xtpoisson command, so use poisson and add in state dummies
poisson pi offer_post i.ym i.state `controls' , exposure(pop_tot) cluster(state)

*Display %change in P&I mortality from Poisson regressions
disp exp(_b[offer_post])-1

*Poisson coefficient is a %change, but all other estimates are presented as changes in levels. Below we get predicted changes in levels.
*get actual predicted number of P&I deaths
predict ed1

*get counterfactual predicted number of P&I deaths assuming no policy
gen offer_post2=offer_post
replace offer_post=0
predict ed2

*generate predicted mortality rates
gen ep1=100000*ed1/pop_tot
gen ep2=100000*ed2/pop_tot

*Weight by state population
egen mean_meanpop = mean(meanpop)
replace ep1 = ep1*meanpop/mean_meanpop
replace ep2 = ep2*meanpop/mean_meanpop

*test for the difference among treated states in the post-treatment period
ttest ep1=ep2 if offer_post2==1 

*Repeat the previous exercise to get the lower and upper bound of the 95% CI. To do this, we replace the coefficient in the b matrix with the lower (upper) bound

*Lower bound CI: 
drop ed1 ed2 ep1 ep2
replace offer_post=offer_post2
*"ereturn repost" will not run if it is not in an eclass program -- the program is totally useless other than to make the command run.
program eclassprogram, eclass
matrix A = e(b)
*Lower bound of the 95% CI from the Poisson regression is -.0400017
matrix A[1,1] = -.0400017
ereturn repost b = A
end
eclassprogram

predict ed1
replace offer_post=0
predict ed2
gen ep1=100000*ed1/pop_tot
gen ep2=100000*ed2/pop_tot
replace ep1 = ep1*meanpop/mean_meanpop
replace ep2 = ep2*meanpop/mean_meanpop
ttest ep1=ep2 if offer_post2==1

*Upper bound CI:
drop ed1 ed2 ep1 ep2
replace offer_post=offer_post2
program eclassprogram2, eclass
matrix A = e(b)
*Upper bound of the 95% CI from the Poisson regression is -.0019753
matrix A[1,1] = -.0019753
ereturn repost b = A
end
eclassprogram2

predict ed1
replace offer_post=0
predict ed2
gen ep1=100000*ed1/pop_tot
gen ep2=100000*ed2/pop_tot
replace ep1 = ep1*meanpop/mean_meanpop
replace ep2 = ep2*meanpop/mean_meanpop
ttest ep1=ep2 if offer_post2==1

*Report some key statistics:

* using XTPOISSON             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
*Poisson coefficient, etc :   |  -.0209885   .0096052    -2.19   0.029    -.0398144   -.0021626
	   
* using POISSON plus i.state  |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval] 
*Poisson coefficient, etc :   |  -.0209885   .0097008    -2.16   0.030    -.0400017   -.0019753

*Poisson % change (i.e., exp(coef)-1: -0.02077
*p-value: 0.029

*predicted difference in PI rate: -0.1322726 
*Lower bound 95% CI:  -.2497235
*Upper bound 95% CI:  -.0125673 








************************************************************************************************
***************Table S3: Effects of Laws on Adolescent and Infant Vaccination Rates (NIS Data)
************************************************************************************************


*************** NIS Teen

global root "DIRECTORY FOR NIS DATA"

* +++++++++++++++++++++++++++++++++++
	* Secondary Globals - Dont change this! - 
cd "$root"
global data     "$root/Data"
global output   "$root/Draft\NIS output"
global dofiles  "$root/Do Files"
global logs 	"$root/Logs"
* +++++++++++++++++++++++++++++++++++

**************************************
* First, rename state and year variables names in dataset of time varying controls, to match analysis file variable names
use "$data\UKCPR_controls.dta", clear

rename fips state
rename year flu_year

tempfile controls
save `controls'

**********************************
use "$data\NISteen_analysis file.dta", clear
gen flu_year = year-1

*merge in time-varying controls
merge m:1 state flu_year using `controls', keep(master matched)

//use dual-frame provider weights for 2011-2014 (note: cell phones added in 2011, & that is the only year with both a landline only & dual frame weight available)
	replace provwt= provwt_d if provwt==.
	replace rddwt=rddwt_d if rddwt==.

*generate state-specific linear time trends at the year level
gen trend = flu_year- 2006
	xi i.state*trend 	
	global linear "_IstaXtren_*"
	
*flag bad match years (<50%) -- note that the mean match rate is 23.74% in these 4 seasons (1998+), and 83.19% in the other 14 seasons.
gen bad_match = 0
	replace bad_match =1 if (flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 )
	
* code up the influenza laws	
do "$dofiles\NIS\code laws.do"

*define vector of time varying controls
global controls "offer_ltc req_child unemploymentrate povertyrate" //need to also add unemp and pov, maybe also MS & HS vaccine mandates?

*********regressions**********
	estimates clear


*All years, with state and influenza year FEs, including time varying controls	
	eststo: reg flu_pastyr  offer_post $controls i.flu_year i.state [aw=provwt],  vce(cluster state)
	estadd ysumm
	estadd local linear ""
	estadd local samp "all years"

*All years, with state and influenza year FEs, including time varying controls	ADDING LINEAR STATE TRENDS
	eststo: reg flu_pastyr  offer_post $controls i.flu_year i.state $linear [aw=provwt],  vce(cluster state)
	estadd ysumm
	estadd local linear "yes"
	estadd local samp "all years"
	
*Excluding bad match years, with state and influenza year FEs, including time varying controls	
	eststo: reg flu_pastyr  offer_post $controls i.flu_year i.state if bad_match==0 [aw=provwt],  vce(cluster state)
	estadd ysumm
	estadd local linear ""
	estadd local samp "exc. bad match yrs"

*Excluding bad match years, with state and influenza year FEs, including time varying controls	ADDING LINEAR STATE TRENDS
	eststo: reg flu_pastyr  offer_post $controls i.flu_year i.state $linear if bad_match==0 [aw=provwt],  vce(cluster state)
	estadd ysumm
	estadd local linear "yes"
	estadd local samp "exc. bad match yrs"
	
	
	esttab est* using "$output\HCW influenza laws_NIS teen.rtf", append keep(offer_post) /*
	*/ p label star(* 0.10 ** 0.05 *** 0.01) stats(N r2 ymean linear samp, label("N" "R-Squared" "Mean of Dependent" "linear trend" "Sample:")) title(DD Spec: Effects of HCW Influenza laws on teen vaccination rates, 2008-2017) /*
	*/ note("All columns are weighted, have state and flu year fixed effects, and include the following controls: $controls. SE are clustered at the state level.")

	

	*************** NIS Child

	
/*Policy notes:  Influenza vaccine was recommended for children aged 6 to 23 months starting with the
2004-05 season (CDC 2003) and for all children ≥ 6 months starting in 2010. */

* First, rename state and year variables names in dataset of time varying controls, to match analysis file variable names
use "$data\UKCPR_controls.dta", clear

rename fips state
rename year flu_year

tempfile controls
save `controls'

*****************************************************
	use "$data\nis_analysisfile.dta", clear

	gen flu_year = year-1

/*Merge in policy data */ 	
	merge m:1 state flu_year using `controls', keep(master matched)
	
	*flag bad match years (<50%) -- note that the mean match rate is 23.74% in these 4 seasons (1998+), and 83.19% in the other 14 seasons.
	gen bad_match = 0
	replace bad_match =1 if (flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 )
	
* code up the influenza laws	
	do "$dofiles\NIS\code laws.do"
	
//use dual weights for 2012-2017, and land-line weights for all other years
	replace provwt = PROVWT_D if year>=2012
	replace provwt = PROVWT_LL if year==2011
	label var provwt "LL weights 2003-2011, dual 2012-2017"
	
	gen provwt_2 = provwt
	replace provwt_2 = PROVWT_D if year==2011
	label var provwt_2 "Uses dual weights for 2011-2017"


*generate state-specific linear time trends at the year level
gen trend = flu_year- 2001
	
	xi i.state*trend 	
	global linear "_IstaXtren_*"

*define vector of time varying controls
global controls "offer_ltc req_child unemploymentrate povertyrate" 

*keep only the youngest age group
keep if agegrp==1

*************** regressions ***********************

*first do estimation for influenza vaccine, then for UTD PCV 
foreach var of varlist flu_any P_UTDPCV {

	estimates clear
*All years, with state and influenza year FEs, including time varying controls		
	eststo: reg `var' offer_post $controls i.flu_year i.state  [aw=provwt],  vce(cluster state)
	estadd ysumm
	estadd local linear ""
	estadd local samp "all years"

*All years, with state and influenza year FEs, including time varying controls	ADDING LINEAR STATE TRENDS
	eststo: reg `var' offer_post $controls  i.flu_year i.state $linear  [aw=provwt],  vce(cluster state)
	estadd ysumm
	estadd local linear "yes"
	estadd local samp "all years"
	
*Excluding bad match years, with state and influenza year FEs, including time varying controls	
	eststo: reg `var' offer_post $controls  i.flu_year i.state if bad_match==0  [aw=provwt],  vce(cluster state)
	estadd ysumm
	estadd local linear ""
	estadd local samp "exc. bad match yrs"
	
*Excluding bad match years, with state and influenza year FEs, including time varying controls	ADDING LINEAR STATE TRENDS
	eststo: reg `var' offer_post $controls   i.flu_year i.state $linear if bad_match==0  [aw=provwt],  vce(cluster state)
	estadd ysumm
	estadd local linear "yes"
	estadd local samp "exc. bad match yrs"
	
	
	esttab est* using "$output\HCW influenza laws_NIS child.rtf", append keep(offer_post) /*
	*/ p label star(* 0.10 ** 0.05 *** 0.01) stats(N r2 ymean linear samp, label("N" "R-Squared" "Mean of Dependent" "linear trend" "Sample:")) title(DD Spec: Effects of HCW Influenza laws on infant vaccination rates BY AGEGRP, 2003-2017) /*
	*/ note("All columns are weighted, have state and flu year fixed effects, and include the following controls: $controls. SE are clustered at the state level.")
}
	}

	*get baseline mean from ever adopters
	eststo: reg flu_any offer_post $controls   i.flu_year i.state $linear if bad_match==0  & agegrp==1 [aw=provwt],  vce(cluster state)
	
	sum flu_any if e(sample)==1 & yearstooffer<0 [aw=provwt] 
	
	eststo: reg P_UTDPCV offer_post $controls   i.flu_year i.state $linear if bad_match==0  & agegrp==1 [aw=provwt],  vce(cluster state)
	
	sum P_UTDPCV if e(sample)==1 & yearstooffer<0 [aw=provwt] 
	



	
	
	

************************************************************************************************
***************Table S4: Effects of Laws on Adult Vaccination Rates (BRFSS Data)
************************************************************************************************


cd "DIRECTORY FOR BRFSS DATA"
clear
use brfss_state_vax_hcw

rename flu flu_vaxrate

cd "$PATH\Influenza National\"
merge 1:1 state flu_year month using flu_hcw_national_analysis
drop if _m!=3
drop _m

collapse (mean) pvac* flu_vaxrate flu_u65 flu_o65 max_offer max_offer_all max_require yearstooffer offer_post meanpop popshare* unemp pov offer_ltc req_child , by(state flu_year)

*drop bad match years (<50%) -- note that the mean match rate is 23.74% in these 4 seasons (1998+), and 83.19% in the other 14 seasons.
drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 

local controls  unemp pov offer_ltc req_child 
xtset state

*Influenza vaccination rate regressions
xtreg flu_vaxrate offer_post  i.flu_year `controls' [aweight=meanpop] , fe cluster(state)

xtreg flu_u65 offer_post  i.flu_year `controls' [aweight=meanpop] , fe cluster(state)

xtreg flu_o65 offer_post  i.flu_year `controls' [aweight=meanpop] , fe cluster(state)

*Pneumococcal vaccination rate regressions
xtreg pvac offer_post  i.flu_year `controls' [aweight=meanpop] , fe cluster(state)

xtreg pvac_u65 offer_post  i.flu_year `controls' [aweight=meanpop] , fe cluster(state)

xtreg pvac_o65 offer_post  i.flu_year `controls' [aweight=meanpop] , fe cluster(state)

















************************************************************************************************
************************************************************************************************
****************************** FIGURES
************************************************************************************************
************************************************************************************************















**************************************************************************
*************** Figure 1
**************************************************************************

*************** Figure 1B

cd "$PATH\Influenza National\"
use flu_hcw_national_analysis.dta, clear

*Generate variable to represent number of people subject to laws in each year-month
gen pop_offer = pop_tot*offer_post
collapse (sum) pi pop_tot pop_offer, by(year month max_offer)

*Collapse to annual
gen flu_year=year
replace flu_year=year-1 if month<7
collapse (mean) pi pop_tot pop_offer, by(flu_year max_offer)

*Generate P&I mortality rate (mean monthly)
gen pi_rate = 100000*pi/pop_tot

*Generate variable for percent of population subject to laws in each year-month, note that this is defined the same for adopting and non-adopting states
*bysort flu_year max_offer: egen sumpop=sum(pop_tot)
bysort flu_year: egen sumpop=sum(pop_tot)
gen pct_pop = pop_offer/sumpop

*Create figure
label var pi_rate "P&I Mortality Rate"
label var flu_year "Influenza Year"

twoway bar pct_pop flu_year if max_offer==1, name(mort_trends, replace) ///
 xlabel(1995 "95/96" 1997 "97/98" 1999 "99/00" 2001 "01/02" 2003 "03/04" 2005 "05/06" 2007 "07/08" 2009 "09/10" 2011 "11/12" 2013 "13/14" 2015 "15/16", angle(20)) ///
 ytitle("Population Share", axis(2)) legend(rows(3))  ylabel(0 "0" .2 "0.2" .4 "0.4" .6 "0.6" .8 "0.8" ,angle(0) axis(2)) ///
 yaxis(2) fcolor(ebblue) fintensity(25) lcolor(ebblue*.25) graphregion(color(white)) ///
|| line pi_rate flu_year if max_offer==1,  lcolor(black) yaxis(1) ytitle("Monthly P&I Mortality Rate", axis(1)) ///
|| line pi_rate flu_year  if max_offer==0, lcolor(gs10) lpattern(shortdash) yaxis(1) ylabel(,angle(0) axis(1)) ///
 legend(cols(1) order( 1 "Share of Population Affected by Laws" 2 "P&I Mortality Rate in Treated States" 3 "P&I Mortality Rate in Non-Treated States")) ///
 xtitle("Influenza Year")
 
*************** Figure 1A

 *Pull in NHIS data
 use "$PATH\Influenza HCW National\Data\Figure 1 Data\nhis_00012.dta", clear
  
  *Keep only employed
  keep if empstat>=10 & empstat<30 
  
 *Generate variables for hospital/health workers
 gen ind_cat = 0
 replace ind_cat = 1 if ind1995==1910 //hospitals
 replace ind_cat = 2 if ind1995==1900|ind1995==1920 //other health services

 *Generate influenza year variable
gen flu_year = year-1 
replace flu_year = year if quarter==4
 
 recode vacflu12m (1=0) (2 3=1) (7 8 9=.)
 recode vacflush12m (1=0) (2=1) (7 8 9=.)
 
 collapse (mean) vacflu12m  vacflush12m [aw=perweight], by(flu_year ind_cat)
 
 merge m:1 flu_year using "$PATH\Influenza HCW National\Data\Figure 1 Data\CPS_collapsed.dta" 
 
 replace frac_treated=0 if ind_cat==0 & frac_treated!=.
 replace frac_treated=0 if ind_cat==2 & frac_treated!=.
   
 drop if flu_year<1997 | flu_year>2016
 label var flu_year "Influenza Year"
  
  twoway bar frac_treated flu_year if ind_cat==1,   title("",color(black)) name(vax_trends, replace) ///
 xlabel(1997 "97/98" 1999 "99/00" 2001 "01/02" 2003 "03/04" 2005 "05/06" 2007 "07/08" 2009 "09/10" 2011 "11/12" 2013 "13/14" 2015 "15/16", angle(20))  ///
 ytitle("Population Share")  ylabel(0 .2 "0.2" .4 "0.4" .6 "0.6" .8 "0.8", angle(horizontal)) ///
 fcolor(ebblue) fintensity(25) lcolor(ebblue*.25) graphregion(color(white)) ///
|| line vacflush12m  flu_year if ind_cat==1,  lcolor(black)  ///
|| line vacflush12m flu_year if ind_cat==0 , lcolor(gs10) lpattern(shortdash)  ///
 legend(cols(1) order( 1 "Share of Hospital Workers Affected by Laws" 2 "Share Vaccinated: Hospital Workers" 3 "Share Vaccinated: Employed, Non-Healthcare")) ///
 xtitle("Influenza Year")


***************Combine Figure 1A and 1B
graph combine vax_trends mort_trends, xsize(4) ysize(2) graphregion(color(white))




**************************************************************************
*************** Figure 2: Event Study
**************************************************************************

clear
use flu_hcw_national_analysis

*yearstooffer represents the number of years to/since law implementation. Note that we have the same
*number of observations from -7 to +3 if we exclude Alabama (they implemented a law in the first year of the sample).
*Because there is no pre-implementation data for Alabama, these estimates are nearly identical if we exclude Alabama.
tab yearstooffer
tab yearstooffer if state!=1

*drop bad match years (<50%) -- note that the mean match rate is 23.74% in these 4 seasons (1998+), and 83.19% in the other 14 seasons.
drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 

*Create a series of indicators for the number of years until law implementation.
forvalues i=2/8{
gen year_pre`i' = (yearstooffer==-`i')
}

*Indicator for -8 captures 8 or more years prior to implementation (included in regression so it is not part of the omitted group, but not displayed)
replace year_pre8 = 1 if yearstooffer<=-8

*Create a series of indicators for the number of years since law implementation.
forvalues i = 0/4{
gen year`i' = (yearstooffer==`i')
}

*Indicator for +4 captures 4 or more years since to implementation (included in regression so it is not part of the omitted group, but not displayed)
replace year4 = 1 if yearstooffer>=4 & !missing(yearstooffer)

xtset state

*****************All Month
*Estimate regression. Omit year_pre1 as reference group. Do not include controls/state-specific trends (event-studies are used to evaluate whether controls/trends are needed) 
xtreg pi_rate year_pre2-year_pre8 year0-year4  i.ym  [aweight=meanpop], fe cluster(state)
*Create and fill in matrix with coefficients/standard errors for figure
matrix D = J(11,3,.)
matrix D[1,1] = _b[year_pre7], _b[year_pre7] - ((invttail(e(df_r),.025))*_se[year_pre7]),_b[year_pre7] + ((invttail(e(df_r),.025))*_se[year_pre7])
matrix D[2,1] = _b[year_pre6], _b[year_pre6] - ((invttail(e(df_r),.025))*_se[year_pre6]),_b[year_pre6] + ((invttail(e(df_r),.025))*_se[year_pre6])
matrix D[3,1] = _b[year_pre5], _b[year_pre5] - ((invttail(e(df_r),.025))*_se[year_pre5]),_b[year_pre5] + ((invttail(e(df_r),.025))*_se[year_pre5])
matrix D[4,1] = _b[year_pre4], _b[year_pre4] - ((invttail(e(df_r),.025))*_se[year_pre4]),_b[year_pre4] + ((invttail(e(df_r),.025))*_se[year_pre4])
matrix D[5,1] = _b[year_pre3], _b[year_pre3] - ((invttail(e(df_r),.025))*_se[year_pre3]),_b[year_pre3] + ((invttail(e(df_r),.025))*_se[year_pre3])
matrix D[6,1] = _b[year_pre2], _b[year_pre2] - ((invttail(e(df_r),.025))*_se[year_pre2]),_b[year_pre2] + ((invttail(e(df_r),.025))*_se[year_pre2])
matrix D[7,1] = 0, 0, 0
matrix D[8,1] = _b[year0], _b[year0] - ((invttail(e(df_r),.025))*_se[year0]), _b[year0] +((invttail(e(df_r),.025))*_se[year0])
matrix D[9,1] = _b[year1], _b[year1] - ((invttail(e(df_r),.025))*_se[year1]), _b[year1] + ((invttail(e(df_r),.025))*_se[year1])
matrix D[10,1] = _b[year2], _b[year2] - ((invttail(e(df_r),.025))*_se[year2]), _b[year2] + ((invttail(e(df_r),.025))*_se[year2])
matrix D[11,1] = _b[year3], _b[year3] - ((invttail(e(df_r),.025))*_se[year3]), _b[year3] + ((invttail(e(df_r),.025))*_se[year3])

*****************Subsample for  Peak Months
xtreg pi_rate year_pre2-year_pre8 year0-year4 i.ym if season==1 [aweight=meanpop], fe cluster(state)
matrix DD = J(11,3,.)
matrix DD[1,1] = _b[year_pre7], _b[year_pre7] - ((invttail(e(df_r),.025))*_se[year_pre7]),_b[year_pre7] + ((invttail(e(df_r),.025))*_se[year_pre7])
matrix DD[2,1] = _b[year_pre6], _b[year_pre6] - ((invttail(e(df_r),.025))*_se[year_pre6]),_b[year_pre6] + ((invttail(e(df_r),.025))*_se[year_pre6])
matrix DD[3,1] = _b[year_pre5], _b[year_pre5] - ((invttail(e(df_r),.025))*_se[year_pre5]),_b[year_pre5] + ((invttail(e(df_r),.025))*_se[year_pre5])
matrix DD[4,1] = _b[year_pre4], _b[year_pre4] - ((invttail(e(df_r),.025))*_se[year_pre4]),_b[year_pre4] + ((invttail(e(df_r),.025))*_se[year_pre4])
matrix DD[5,1] = _b[year_pre3], _b[year_pre3] - ((invttail(e(df_r),.025))*_se[year_pre3]),_b[year_pre3] + ((invttail(e(df_r),.025))*_se[year_pre3])
matrix DD[6,1] = _b[year_pre2], _b[year_pre2] - ((invttail(e(df_r),.025))*_se[year_pre2]),_b[year_pre2] + ((invttail(e(df_r),.025))*_se[year_pre2])
matrix DD[7,1] = 0, 0, 0
matrix DD[8,1] = _b[year0], _b[year0] - ((invttail(e(df_r),.025))*_se[year0]), _b[year0] +((invttail(e(df_r),.025))*_se[year0])
matrix DD[9,1] = _b[year1], _b[year1] - ((invttail(e(df_r),.025))*_se[year1]), _b[year1] + ((invttail(e(df_r),.025))*_se[year1])
matrix DD[10,1] = _b[year2], _b[year2] - ((invttail(e(df_r),.025))*_se[year2]), _b[year2] + ((invttail(e(df_r),.025))*_se[year2])
matrix DD[11,1] = _b[year3], _b[year3] - ((invttail(e(df_r),.025))*_se[year3]), _b[year3] + ((invttail(e(df_r),.025))*_se[year3])

*****************Subsample for Non-Peak Season
xtreg pi_rate year_pre2-year_pre8 year0-year4 i.ym  if season==0 [aweight=meanpop], fe cluster(state)
matrix DDD = J(11,3,.)
matrix DDD[1,1] = _b[year_pre7], _b[year_pre7] - ((invttail(e(df_r),.025))*_se[year_pre7]),_b[year_pre7] + ((invttail(e(df_r),.025))*_se[year_pre7])
matrix DDD[2,1] = _b[year_pre6], _b[year_pre6] - ((invttail(e(df_r),.025))*_se[year_pre6]),_b[year_pre6] + ((invttail(e(df_r),.025))*_se[year_pre6])
matrix DDD[3,1] = _b[year_pre5], _b[year_pre5] - ((invttail(e(df_r),.025))*_se[year_pre5]),_b[year_pre5] + ((invttail(e(df_r),.025))*_se[year_pre5])
matrix DDD[4,1] = _b[year_pre4], _b[year_pre4] - ((invttail(e(df_r),.025))*_se[year_pre4]),_b[year_pre4] + ((invttail(e(df_r),.025))*_se[year_pre4])
matrix DDD[5,1] = _b[year_pre3], _b[year_pre3] - ((invttail(e(df_r),.025))*_se[year_pre3]),_b[year_pre3] + ((invttail(e(df_r),.025))*_se[year_pre3])
matrix DDD[6,1] = _b[year_pre2], _b[year_pre2] - ((invttail(e(df_r),.025))*_se[year_pre2]),_b[year_pre2] + ((invttail(e(df_r),.025))*_se[year_pre2])
matrix DDD[7,1] = 0, 0, 0
matrix DDD[8,1] = _b[year0], _b[year0] - ((invttail(e(df_r),.025))*_se[year0]), _b[year0] +((invttail(e(df_r),.025))*_se[year0])
matrix DDD[9,1] = _b[year1], _b[year1] - ((invttail(e(df_r),.025))*_se[year1]), _b[year1] + ((invttail(e(df_r),.025))*_se[year1])
matrix DDD[10,1] = _b[year2], _b[year2] - ((invttail(e(df_r),.025))*_se[year2]), _b[year2] + ((invttail(e(df_r),.025))*_se[year2])
matrix DDD[11,1] = _b[year3], _b[year3] - ((invttail(e(df_r),.025))*_se[year3]), _b[year3] + ((invttail(e(df_r),.025))*_se[year3])

*Make All Month, and Peak/Non-Peak plots separately

 coefplot (matrix(D[,1]),  ci((D[,2] D[,3])) lcolor(black) mcolor(black) mlcolor(black) msymbol(D) ylabel(-1 "-1" -.5 "-0.5" 0 "0" .5 "0.5" 1 "1", angle(0)) ///
 mlwidth(thin) connect(direct) ciopts(recast(rcap) lcolor(black))) ,  xlabel(1 "-7" 2 "-6" 3 "-5" 4 "-4" 5 "-3" 6 "-2" 7 "-1" 8 "0" 9 "1" 10 "2" 11 "3" ) ///
 leg(on)  xtitle("Influenza-Years Relative to Law") ytitle("Change in Monthly P&I Deaths per 100,000 Population") graphregion(color(white)) vertical  yline(0, lcolor(black) lwidth(medthin)) xline(8, lcolor(gs10) lwidth(thin)) name(allmonth, replace)

 coefplot (matrix(DD[,1]), ci((DD[,2] DD[,3])) lcolor(orange) mcolor(orange) mlcolor(orange) msymbol(O) msize(medlarge) offset(.1) ylabel(-1 "-1" -.5 "-0.5" 0 "0" .5 "0.5" 1 "1", angle(0))  ///
 mlwidth(med) connect(direct) ciopts(recast(rcap) lcolor(orange))) ///
  (matrix(DDD[,1]), ci((DDD[,2] DDD[,3])) lcolor(midblue) mcolor(midblue) mlcolor(midblue) msymbol(T) offset(-.1) ///
 mlwidth(med) connect(direct) ciopts(recast(rcap) lcolor(midblue))) ///
  ,  xlabel(1 "-7" 2 "-6" 3 "-5" 4 "-4" 5 "-3" 6 "-2" 7 "-1" 8 "0" 9 "1" 10 "2" 11 "3") ///
  xtitle("Influenza-Years Relative to Law") ytitle("Change in Monthly P&I Deaths per 100,000 Population") graphregion(color(white)) vertical  yline(0, lcolor(black) lwidth(medthin)) xline(8, lcolor(gs10) lwidth(thin)) name(byseason, replace)

*Combine plots, edit legend manually.
 graph combine allmonth byseason, ycomm graphregion(color(white)) row(1) ysize(2) xsize(4)
  

  
  
  
**************************************************************************
*************** Figure S1: Vaccination laws over time
**************************************************************************

clear 
set obs 56
gen state_fip =_n
*offer means any offer law or requirement
gen offeryear = . 
replace offeryear =1994 if state==1 /*AL -- not exactly an offer law*/
replace offeryear =2007 if state==6 /*CA - hospital*/
replace offeryear =2012 if state==8 /*CO - hospital, ASC, LTCF - requirement*/
replace offeryear =2008 if state==11 /*DC - All HCF - supposed requirement*/
replace offeryear =2010 if state==13 /*GA - Hospital*/
replace offeryear =2010 if state==17 /*IL - All HCF*/
replace offeryear =2002 if state==23 /*ME - All HCF*/
replace offeryear =2008 if state==24 /*MD - Hospital*/
replace offeryear =2009 if state==25 /*MA - Hospital*/
replace offeryear =2011 if state==31 /*NE - Hospital -- language is super weak*/
replace offeryear =2005 if state==33 /*NH - Hospital, residential care, adult day care, assisted living -- some evidence suggests it wasn't well implemented*/
replace offeryear =2013 if state==36 /*NY - All HCF - seems like a requirement*/
replace offeryear =2009 if state==40 /*OK - Hospital*/
replace offeryear =2012 if state==44 /*RI - All HCF -- seems like a requirement*/
replace offeryear =2007 if state==47 /*TN - Hospital*/

keep if offeryear!=.

gen adoptorder = . 
replace adoptorder =1 if state==1 /*AL -- not exactly an offer law*/
replace adoptorder =2 if state==23 /*ME - All HCF*/
replace adoptorder =3 if state==33 /*NH - Hospital, residential care, adult day care, assisted living -- some evidence suggests it wasn't well implemented*/
replace adoptorder =4 if state==6 /*CA - hospital*/
replace adoptorder =5 if state==47 /*TN - Hospital*/
replace adoptorder =6 if state==11 /*DC - All HCF - supposed requirement*/
replace adoptorder =7 if state==24 /*MD - Hospital*/
replace adoptorder =8 if state==25 /*MA - Hospital*/
replace adoptorder =9 if state==40 /*OK - Hospital*/
replace adoptorder =10 if state==13 /*GA - Hospital*/
replace adoptorder =11 if state==17 /*IL - All HCF*/
replace adoptorder =12 if state==31 /*NE - Hospital -- language is super weak*/
replace adoptorder =13 if state==8 /*CO - hospital, ASC, LTCF - requirement*/
replace adoptorder =14 if state==44 /*RI - All HCF -- seems like a requirement*/
replace adoptorder =15 if state==36 /*NY - All HCF - seems like a requirement*/

gen endyear=2018

*reverse order, so 1st adopting state gets plotted 1st
gsort -adoptorder
gen adoptorder2= _n

global state_labels2 "1 "NEW YORK" 2 "RHODE ISLAND" 3 "COLORADO" 4 "NEBRASKA" 5 "ILLINOIS" 6 "GEORGIA" 7 "OKLAHOMA" 8 "MASSACHUSETTS" 9 "MARYLAND" 10 "D.C." 11 "TENNESSEE" 12 "CALIFORNIA" 13 "NEW HAMPSHIRE" 14 "MAINE" 15 "ALABAMA""

twoway pcspike adoptorder2 offeryear adoptorder2 endyear, mlabels(state) || ///
	pci 0 1995 16 1995  "Start", lpattern(dash) lcolor(dkgreen) lwidth(medthick) || ///
	pci 0 2017	16 2017	 "End",	lpattern(dash) lcolor(red) lwidth(thin) ///
	ylabel($state_labels2 , labsize(vsmall) angle(horizontal)) ///
	xlabel(1994 "1994/95" 1998 "1998/99 " 2002 "2002/03" 2006 "2006/07" 2010 "2010/11" 2014 "2014/15" 2018 "2018/19", angle(20))  ///
 	graphr(color(white)) lwidth(medthick) ///
	ytitle("") xtitle("Influenza Year") ///
	legend(cols(3) order(1 2 3 ) label(1 "Law in effect" ) label(2 "Start of Mortality Data") label(3 "End of Mortality Data") size(vsmall)) ///
	subtitle("Hospital Worker Influenza Vaccination Laws, by State")
	


**************************************************************************
*************** Figure S2: Predicted Trends
**************************************************************************

*Import analysis file
cd "$PATH\Influenza National\"
clear
use flu_hcw_national_analysis

*Preserve data that includes the full sample (for graphs), but run regression only on the analysis sample
preserve
drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 
local controls  temp* hum* prcp* popshare* unemp pov offer_ltc req_child i.state#c.ym
xtset state
reg pi_rate offer_post  i.ym i.state `controls' [aweight=meanpop] ,  cluster(state)
restore

*Get actual predicted values; don't make predictions for the years that are dropped from the sample 
predict predicted_pi if (flu_year!=1997 & flu_year!=2003 & flu_year!=2004 & flu_year !=2007 & flu_year!=2014 )

*Get counterfactual predicted values
gen offer_post2=offer_post
replace offer_post=0
predict predicted_pi_counterfactual if (flu_year!=1997 & flu_year!=2003 & flu_year!=2004 & flu_year !=2007 & flu_year!=2014 )
ttest predicted_pi=predicted_pi_counterfactual if offer_post2==1

*Generate a variable for the percent of the population affected
gen pop_offer = pop_tot*offer_post2
collapse (mean) predicted* (rawsum) pi pop_tot pop_offer [aw=meanpop], by(year month max_offer)

gen flu_year=year
replace flu_year=year-1 if month<7
collapse (mean) pi pop_tot pop_offer predicted*, by(flu_year max_offer)

*Generate variable for percent of population subject to laws in each year-month, note that this is defined the same for adopting and non-adopting states
*bysort flu_year max_offer: egen sumpop=sum(pop_tot)
bysort flu_year: egen sumpop=sum(pop_tot)
gen pct_pop = pop_offer/sumpop

*Create figure
label var predicted_pi "P&I Mortality Rate"
label var flu_year "Influenza Year"

 twoway bar pct_pop flu_year if max_offer==1, name(mort_trends, replace) ///
 xlabel(1995 "95/96" 1997 "97/98" 1999 "99/00" 2001 "01/02" 2003 "03/04" 2005 "05/06" 2007 "07/08" 2009 "09/10" 2011 "11/12" 2013 "13/14" 2015 "15/16", angle(20)) ///
 ytitle("Population Share", axis(2)) legend(rows(3))  ylabel(0 "0" .2 "0.2" .4 "0.4" .6 "0.6" .8 "0.8" ,angle(0) axis(2)) ///
 yaxis(2) fcolor(ebblue) fintensity(25) lcolor(ebblue*.25) graphregion(color(white)) ///
|| line predicted_pi flu_year if max_offer==1,  lcolor(black) yaxis(1) ytitle("Monthly P&I Mortality Rate", axis(1)) ///
|| line predicted_pi_counterfactual flu_year  if max_offer==1, lcolor(black) lpattern(shortdash) yaxis(1) ylabel(,angle(0) axis(1)) ///
 legend(cols(1) rows(4) order( 1 "Share of Population Affected by Laws" 2 "Predicted Outcome in Treated States" 3 "Counterfactual Predicted Outcome in Treated States")) ///
 xtitle("Influenza Year")
 
 

 
 
 
 
 
 
**************************************************************************
*************** Figure S3: Synthetic Control Trends by State
**************************************************************************

*Begin with synthetic control analysis

*first, generate mean state population, which will be used later
cd "$PATH\Data\Population\1969-2017 SEER"
use Pop_Age_State_SEER_1969_2017_monthly, clear
collapse (mean) pop_tot, by(state)
save mean_state_pop.dta, replace

cd "$PATH\Influenza National\"

clear
use flu_hcw_national_analysis

*Collapse to annual level
collapse (sum) pi pop_tot offer_post max_offer, by(state flu_year)

*Re-generate offyear variable
gen offeryear = . 
replace offeryear =1994 if state==1 /*AL -- not exactly an offer law*/
replace offeryear =2007 if state==6 /*CA - hospital*/
replace offeryear =2012 if state==8 /*CO - hospital, ASC, LTCF - requirement*/
replace offeryear =2008 if state==11 /*DC - All HCF - supposed requirement*/
replace offeryear =2010 if state==13 /*GA - Hospital*/
replace offeryear =2010 if state==17 /*IL - All HCF*/
replace offeryear =2002 if state==23 /*ME - All HCF*/
replace offeryear =2008 if state==24 /*MD - Hospital*/
replace offeryear =2009 if state==25 /*MA - Hospital*/
replace offeryear =2011 if state==31 /*NE - Hospital -- language is super weak*/
replace offeryear =2005 if state==33 /*NH - Hospital, residential care, adult day care, assisted living -- some evidence suggests it wasn't well implemented*/
replace offeryear =2013 if state==36 /*NY - All HCF - seems like a requirement*/
replace offeryear =2009 if state==40 /*OK - Hospital*/
replace offeryear =2012 if state==44 /*RI - All HCF -- seems like a requirement*/
replace offeryear =2007 if state==47 /*TN - Hospital*/

*Re-generate P&I mortality rate
gen pi_rate = 100000*pi/pop_tot

*Set panel and time variables (for synth command)
tsset state flu_year

*Set up file to store regression estimates for each state
cap postclose synth_estimates
postfile synth_estimates beta ste state str20 stateabr using synth_estimates.dta, replace


*Each state is treated separately, CA is first and is commented, the other states follow the same outline

*CA
preserve
*Keep only CA and non-adopting states
keep if (state==6 | max_offer==0)

*Run synth command, match on the outcome in each pre-treatment year (not average of the outcome over all pre-treatment years).
synth pi_rate pi_rate(2006) pi_rate(2005) pi_rate(2004) pi_rate(2003) pi_rate(2002) pi_rate(2001) pi_rate(2000) pi_rate(1999) ///
pi_rate(1998) pi_rate(1997) pi_rate(1996) pi_rate(1995), trunit(6) trperiod(2007) keep(synth_6, replace)

*Pull in main data to run synthetic control regressions
clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_6 
drop _m
replace _W_Weight = 1 if state==6
*Simple regression
xtset state
xtreg pi_rate offer_post  i.ym  [aweight=_W_Weight] , fe cluster(state)
post synth_estimates (_b[offer_post]) (_se[offer_post]) (6) ("CA")
*Use file created by synth command to manually plot synthetic control figure
clear
use synth_6
twoway line _Y_t _time, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) /// 
xtitle("")  title("California", color(black)) name(CA, replace) ///
||  line _Y_s _time, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(2007,  lcolor(gs10) lwidth(medthin))
restore 

*CO
preserve
keep if (state==8 | max_offer==0)
synth pi_rate pi_rate(2011) pi_rate(2010) pi_rate(2009) pi_rate(2008) pi_rate(2007) pi_rate(2006) pi_rate(2005) ///
pi_rate(2004) pi_rate(2003) pi_rate(2002) pi_rate(2001) pi_rate(2000) pi_rate(1999) ///
pi_rate(1998) pi_rate(1997) pi_rate(1996) pi_rate(1995) , trunit(8) trperiod(2012) keep(synth_8, replace)

clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_8 
drop _m
replace _W_Weight = 1 if state==8
xtset state
xtreg pi_rate offer_post  i.ym  [aweight=_W_Weight] , fe cluster(state)
post synth_estimates (_b[offer_post]) (_se[offer_post]) (8) ("CO")

clear
use synth_8
twoway line _Y_t _time, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) /// 
xtitle("")  title("Colorado", color(black)) name(CO, replace) ///
||  line _Y_s _time, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(2012,  lcolor(gs10) lwidth(medthin))
restore 

*DC
preserve
keep if (state==11 | max_offer==0)
synth pi_rate pi_rate(2007) pi_rate(2006) pi_rate(2005) ///
pi_rate(2004) pi_rate(2003) pi_rate(2002) pi_rate(2001) pi_rate(2000) pi_rate(1999) ///
pi_rate(1998) pi_rate(1997) pi_rate(1996) pi_rate(1995)  , trunit(11) trperiod(2008) keep(synth_11, replace)

clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_11
drop _m
replace _W_Weight = 1 if state==11
xtset state
xtreg pi_rate offer_post  i.ym  [aweight=_W_Weight] , fe cluster(state)
post synth_estimates (_b[offer_post]) (_se[offer_post]) (11) ("DC")

clear
use synth_11
twoway line _Y_t _time, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) /// 
xtitle("")  title("DC", color(black)) name(DC, replace) ///
||  line _Y_s _time, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(2008,  lcolor(gs10) lwidth(medthin))
restore 

*GA
preserve
keep if (state==13 | max_offer==0)
synth pi_rate pi_rate(2009) pi_rate(2008) pi_rate(2007) pi_rate(2006) pi_rate(2005) ///
pi_rate(2004) pi_rate(2003) pi_rate(2002) pi_rate(2001) pi_rate(2000) pi_rate(1999) ///
pi_rate(1998) pi_rate(1997) pi_rate(1996) pi_rate(1995) , trunit(13) trperiod(2010) keep(synth_13, replace)

clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_13
drop _m
replace _W_Weight = 1 if state==13
xtset state
xtreg pi_rate offer_post  i.ym  [aweight=_W_Weight] , fe cluster(state)
post synth_estimates (_b[offer_post]) (_se[offer_post]) (13) ("GA")

clear
use synth_13
twoway line _Y_t _time, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) /// 
xtitle("")  title("Georgia", color(black)) name(GA, replace) ///
||  line _Y_s _time, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(2010,  lcolor(gs10) lwidth(medthin))
restore 

*IL
preserve
keep if (state==17 | max_offer==0)
synth pi_rate pi_rate(2009) pi_rate(2008) pi_rate(2007) pi_rate(2006) pi_rate(2005) ///
pi_rate(2004) pi_rate(2003) pi_rate(2002) pi_rate(2001) pi_rate(2000) pi_rate(1999) ///
pi_rate(1998) pi_rate(1997) pi_rate(1996) pi_rate(1995) , trunit(17) trperiod(2010) keep(synth_17, replace)

clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_17
drop _m
replace _W_Weight = 1 if state==17
xtset state
xtreg pi_rate offer_post  i.ym  [aweight=_W_Weight] , fe cluster(state)
post synth_estimates (_b[offer_post]) (_se[offer_post]) (17) ("IL")

clear
use synth_17
twoway line _Y_t _time, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) /// 
xtitle("")  title("Illinois", color(black)) name(IL, replace) ///
||  line _Y_s _time, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(2010,  lcolor(gs10) lwidth(medthin))
restore 

*ME
preserve
keep if (state==23 | max_offer==0)
synth pi_rate pi_rate(2001) pi_rate(2000) pi_rate(1999) ///
pi_rate(1998) pi_rate(1997) pi_rate(1996) pi_rate(1995) , trunit(23) trperiod(2002) keep(synth_23, replace)

clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_23
drop _m
replace _W_Weight = 1 if state==23
xtset state
xtreg pi_rate offer_post  i.ym  [aweight=_W_Weight] , fe cluster(state)
post synth_estimates (_b[offer_post]) (_se[offer_post]) (23) ("ME")

clear
use synth_23
twoway line _Y_t _time, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) /// 
xtitle("")  title("Maine", color(black)) name(ME, replace) ///
||  line _Y_s _time, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(2002,  lcolor(gs10) lwidth(medthin))
restore 

*MD
preserve
keep if (state==24 | max_offer==0)
synth pi_rate pi_rate(2007) pi_rate(2006) pi_rate(2005) ///
pi_rate(2004) pi_rate(2003) pi_rate(2002) pi_rate(2001) pi_rate(2000) pi_rate(1999) ///
pi_rate(1998) pi_rate(1997) pi_rate(1996) pi_rate(1995) , trunit(24) trperiod(2008) keep(synth_24, replace)

clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_24
drop _m
replace _W_Weight = 1 if state==24
xtset state
xtreg pi_rate offer_post  i.ym  [aweight=_W_Weight] , fe cluster(state)
post synth_estimates (_b[offer_post]) (_se[offer_post]) (24) ("MD")

clear
use synth_24
twoway line _Y_t _time, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) /// 
xtitle("")  title("Maryland", color(black)) name(MD, replace) ///
||  line _Y_s _time, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(2008,  lcolor(gs10) lwidth(medthin))
restore 

*MA
preserve
keep if (state==25 | max_offer==0)
synth pi_rate pi_rate(2008) pi_rate(2007) pi_rate(2006) pi_rate(2005) ///
pi_rate(2004) pi_rate(2003) pi_rate(2002) pi_rate(2001) pi_rate(2000) pi_rate(1999) ///
pi_rate(1998) pi_rate(1997) pi_rate(1996) pi_rate(1995) , trunit(25) trperiod(2009) keep(synth_25, replace)

clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_25
drop _m
replace _W_Weight = 1 if state==25
xtset state
xtreg pi_rate offer_post  i.ym  [aweight=_W_Weight] , fe cluster(state)
post synth_estimates (_b[offer_post]) (_se[offer_post]) (25) ("MA")

clear
use synth_25
twoway line _Y_t _time, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) /// 
xtitle("")  title("Massachusetts", color(black)) name(MA, replace) ///
||  line _Y_s _time, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(2009,  lcolor(gs10) lwidth(medthin))
restore 

*NE
preserve
keep if (state==31 | max_offer==0)
synth pi_rate  pi_rate(2010) pi_rate(2009) pi_rate(2008) pi_rate(2007) pi_rate(2006) pi_rate(2005) ///
pi_rate(2004) pi_rate(2003) pi_rate(2002) pi_rate(2001) pi_rate(2000) pi_rate(1999) ///
pi_rate(1998) pi_rate(1997) pi_rate(1996) pi_rate(1995) , trunit(31) trperiod(2011) keep(synth_31, replace)

clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_31
drop _m
replace _W_Weight = 1 if state==31
xtset state
xtreg pi_rate offer_post  i.ym  [aweight=_W_Weight] , fe cluster(state)
post synth_estimates (_b[offer_post]) (_se[offer_post]) (31) ("NE")

clear
use synth_31
twoway line _Y_t _time, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) /// 
xtitle("")  title("Nebraska", color(black)) name(NE, replace) ///
||  line _Y_s _time, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(2011,  lcolor(gs10) lwidth(medthin))
restore 

*NH
preserve
keep if (state==33 | max_offer==0)
synth pi_rate pi_rate(2004) pi_rate(2003) pi_rate(2002) pi_rate(2001) pi_rate(2000) pi_rate(1999) ///
pi_rate(1998) pi_rate(1997) pi_rate(1996) pi_rate(1995) , trunit(33) trperiod(2005) keep(synth_33, replace)

clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_33
drop _m
replace _W_Weight = 1 if state==33
xtset state
xtreg pi_rate offer_post  i.ym  [aweight=_W_Weight] , fe cluster(state)
post synth_estimates (_b[offer_post]) (_se[offer_post]) (33) ("NH")

clear
use synth_33
twoway line _Y_t _time, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) /// 
xtitle("")  title("New Hampshire", color(black)) name(NH, replace) ///
||  line _Y_s _time, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(2005,  lcolor(gs10) lwidth(medthin))
restore 

*NY
preserve
keep if (state==36 | max_offer==0)
synth pi_rate pi_rate(2012) pi_rate(2011) pi_rate(2010) pi_rate(2009) pi_rate(2008) pi_rate(2007) pi_rate(2006) pi_rate(2005) ///
pi_rate(2004) pi_rate(2003) pi_rate(2002) pi_rate(2001) pi_rate(2000) pi_rate(1999) ///
pi_rate(1998) pi_rate(1997) pi_rate(1996) pi_rate(1995) , trunit(36) trperiod(2013) keep(synth_36, replace)

clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_36
drop _m
replace _W_Weight = 1 if state==36
xtset state
xtreg pi_rate offer_post  i.ym  [aweight=_W_Weight] , fe cluster(state)
post synth_estimates (_b[offer_post]) (_se[offer_post]) (36) ("NY")

clear
use synth_36
twoway line _Y_t _time, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) /// 
 xtitle("") title("New York", color(black)) name(NY, replace) ///
||  line _Y_s _time, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(2013,  lcolor(gs10) lwidth(medthin))
restore 

*OK
preserve
keep if (state==40 | max_offer==0)
synth pi_rate pi_rate(2008) pi_rate(2007) pi_rate(2006) pi_rate(2005) ///
pi_rate(2004) pi_rate(2003) pi_rate(2002) pi_rate(2001) pi_rate(2000) pi_rate(1999) ///
pi_rate(1998) pi_rate(1997) pi_rate(1996) pi_rate(1995) , trunit(40) trperiod(2009) keep(synth_40, replace)

clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_40
drop _m
replace _W_Weight = 1 if state==40
xtset state
xtreg pi_rate offer_post  i.ym  [aweight=_W_Weight] , fe cluster(state)
post synth_estimates (_b[offer_post]) (_se[offer_post]) (40) ("OK")

clear
use synth_40
twoway line _Y_t _time, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) /// 
xtitle("")  title("Oklahoma", color(black)) name(OK, replace) ///
||  line _Y_s _time, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(2009,  lcolor(gs10) lwidth(medthin))
restore 

*RI
preserve
keep if (state==44 | max_offer==0)
synth pi_rate pi_rate(2011) pi_rate(2010) pi_rate(2009) pi_rate(2008) pi_rate(2007) pi_rate(2006) pi_rate(2005) ///
pi_rate(2004) pi_rate(2003) pi_rate(2002) pi_rate(2001) pi_rate(2000) pi_rate(1999) ///
pi_rate(1998) pi_rate(1997) pi_rate(1996) pi_rate(1995) , trunit(44) trperiod(2012) keep(synth_44, replace)

clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_44
drop _m
replace _W_Weight = 1 if state==44
xtset state
xtreg pi_rate offer_post  i.ym  [aweight=_W_Weight] , fe cluster(state)
post synth_estimates (_b[offer_post]) (_se[offer_post]) (44) ("RI")

clear
use synth_44
twoway line _Y_t _time, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) /// 
xtitle("")  title("Rhode Island", color(black)) name(RI, replace) ///
||  line _Y_s _time, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(2012,  lcolor(gs10) lwidth(medthin))
restore 

*TN
preserve
keep if (state==47 | max_offer==0)
synth pi_rate pi_rate(2006) pi_rate(2005) ///
pi_rate(2004) pi_rate(2003) pi_rate(2002) pi_rate(2001) pi_rate(2000) pi_rate(1999) ///
pi_rate(1998) pi_rate(1997) pi_rate(1996) pi_rate(1995) , trunit(47) trperiod(2007) keep(synth_47, replace)

clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_47
drop _m
replace _W_Weight = 1 if state==47
xtset state
xtreg pi_rate offer_post  i.ym  [aweight=_W_Weight] , fe cluster(state)
post synth_estimates (_b[offer_post]) (_se[offer_post]) (47) ("TN")

clear
use synth_47
twoway line _Y_t _time, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) /// 
xtitle("")  title("Tennessee", color(black)) name(TN, replace) ///
||  line _Y_s _time, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(2007,  lcolor(gs10) lwidth(medthin))
postclose synth_estimates


*Plot all synthetic control plots on same graph, locate plots and edit legend manually
grc1leg2 ME NH CA TN MD DC MA OK GA IL NE CO RI NY, graphregion(color(white)) 




**************************************************************************
*************** Figure S4: Synthetic Control Treatment Effects
**************************************************************************

*Open postfile of regression estimates
clear
use synth_estimates

*Merge in mean population data for weighted estimates
merge m:1 state using  "$PATH\Data\Population\1969-2017 SEER\mean_state_pop"
keep if _m==3

*Get average (unweighted and weighted) regression estimates
sum beta 
sum beta [aw=pop_tot]

*Calculate average estimates for each law type
gen lawcat = 0
replace lawcat = 1 if stateab=="NH" | stateab=="DC" | stateab=="GA" | stateab=="NE" 
replace lawcat = 2 if stateab=="ME" | stateab=="CA" | stateab=="TN" | stateab=="MD" | stateab=="MA"  | stateab=="OK"  | stateab=="IL"  
replace lawcat = 3 if stateab=="CO" | stateab=="RI" | stateab=="NY" 
sum beta if lawcat==1 [aw=pop_tot]
sum beta if lawcat==2 [aw=pop_tot] 
sum beta if lawcat==3 [aw=pop_tot]

*Plot effect sizes, delete y-axis and add text manually
gsort  - beta
gen n= _n

*Plot effect sizes with different colors for different law types
gen betalabel = beta
replace betalabel = beta-.15 if beta<0

*Color code by law category, edit legend and add numbers manually.
twoway bar beta n if lawcat==1 , graphregion(color(white)) legend(row(4)) xtitle("") xlab("") xscale(lstyle(none)) ///
 ytitle("Synthetic Control Coefficient Estimate") yline(-0.464, lcolor(gs6) lpattern(dash))  yline(0, lcolor(black)) ///
  color(ebblue) fintensity(25) lcolor(ebblue*.25) barw(.65) ylabel(-2 "-2" -1.5 "-1.5" -1 "-1" -.5 "-0.5" 0 "0",angle(0)) ///
 || bar beta n if lawcat==2 , color(ebblue) fintensity(75) lcolor(ebblue*.75) barw(.65) ///
 || bar beta n if lawcat==3 , color(ebblue) fintensity(150) lcolor(ebblue*1.5) barw(.65) ///
 || scat betalabel n , mlab(stateab) msym(none) mlabpos(12) mlabcolor(black) mlabsize(medsmall)


 
 
 
 
**************************************************************************
*************** Figure S5: Match Rate Figure
**************************************************************************

*Import analysis file
clear
cd "$PATH\Influenza National\"
use flu_hcw_national_analysis.dta

*Create bar graph -- need to make low match bars red manually
 twoway (bar matchrate flu_year, barw(.9)  xlabel(1995 "95/96" 1997 "97/98" 1999 "99/00" 2001 "01/02" 2003 "03/04" 2005 "05/06" 2007 "07/08" 2009 "09/10" 2011 "11/12" 2013 "13/14" 2015 "15/16", angle(20)) ///
 fcolor(ebblue) fintensity(25) lcolor(ebblue) graphregion(color(white)) ylabel(0 "0" .2 "0.2" .4 "0.4" .6 "0.6" .8 "0.8" 1 "1", angle(0)) xtitle("Influenza Season")) 

 
 
 
 
 
 **************************************************************************
*************** Figure S6: P&I Mortality Trends by Year of Implementation
**************************************************************************
 
 cd "$PATH\Influenza National\"
use flu_hcw_national_analysis.dta, clear

gen yearbin = 0
*2002-2005
replace yearbin = 1 if state==23 |  state==33 
*2006-2009 
replace yearbin = 2 if state==6 | state==47 | state==24 | state==11 | state==25 | state==40 
*2010-2013
replace yearbin = 3 if state==13 | state==17 | state==31 |  state==8 | state==44 | state==36 

*Generate variable to represent number of people subject to laws in each year-month
gen pop_offer = pop_tot*offer_post
collapse (sum) pi pop_tot pop_offer, by(year month yearbin)

*Collapse to annual
gen flu_year=year
replace flu_year=year-1 if month<7
collapse (mean) pi pop_tot pop_offer, by(flu_year yearbin)

*Generate P&I mortality rate (mean monthly)
gen pi_rate = 100000*pi/pop_tot

*Generate variable for percent of population subject to laws in each year-month, note that this is defined the same for adopting and non-adopting states
bysort flu_year yearbin: egen sumpop=sum(pop_tot)
gen pct_pop = pop_offer/sumpop

*Create figure
label var pi_rate "P&I Mortality Rate"
label var flu_year "Influenza Year"

*Edit legend manually
twoway line pi_rate flu_year if yearbin==0,  ///
 xlabel(1995 "95/96" 1997 "97/98" 1999 "99/00" 2001 "01/02" 2003 "03/04" 2005 "05/06" 2007 "07/08" 2009 "09/10" 2011 "11/12" 2013 "13/14" 2015 "15/16", angle(20)) ///
 ytitle("Monthly P&I Mortality Rate") ylabel(,angle(0)) legend(rows(3))  lcolor(black) graphregion(color(white)) lwidth(thick) ///
|| connected pi_rate flu_year  if yearbin==1, lcolor(gs10) mcolor(blue) msize(medsmall) msymb(Oh) ///
|| connected pi_rate flu_year  if yearbin==2, lcolor(gs10) mcolor(red) msize(small)  msymb(T) ///
|| connected pi_rate flu_year  if yearbin==3, lcolor(gs10) mcolor(midgreen) msize(small)  msymb(S) 





**************************************************************************
****************Figure S6 and S7: Estimates by 39 causes of death.
**************************************************************************

clear
use flu_hcw_national_analysis

*Drop pre-1999/00 (UCR39 not defined pre-99)
drop if flu_year<1999

drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 

*Create matrix to fill in with results
mat X = J(39,3,.)

*Use postfile command to create a new file that will be filled in with estimates
cap postclose mtpvalues
postfile mtpvalues beta ste str20 estm using mtpvalues.dta, replace

*Loop over all 39 causes of death (note that Cause 24 is P&I)
local i =1
foreach var in 24 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39{

*Note the UCR_PI_24 variable is the number of deaths with any primary/secondary diagnosis for P&I
*UCR_PI_* for other categories represents deaths with another primary category and no secondary diagosis for P&I
replace UCR_PI_`var' = 100000*UCR_PI_`var'/pop_tot

sum UCR_PI_`var', meanonly
local UCR_mean = r(mean)

local controls temp* hum* prcp* popshare* unemp pov offer_ltc req_child i.state#c.ym
xtset state

*Run regression, for peak months only where largest effect for P&I is expected
xtreg UCR_PI_`var' offer_post  i.ym `controls' if season==1 [aweight=meanpop], fe cluster(state)
lincom offer_post
mat X[`i',1] = r(estimate), r(estimate) -(invttail(e(df_r),.025))*r(se), r(estimate) + (invttail(e(df_r),.025))*r(se)
post mtpvalues (r(estimate)) (r(se)) ("est_`var'")

local i = `i'+1
}

postclose mtpvalues

*Open file with estimates, create p-values for each estimate. Create adjusted p-values correcting for multiple hypothesis testing
*Use a variety of adjustments. Note that bonferroni is the most conservative and simes is the most liberal among these.
clear 
use mtpvalues
gen z = beta/ste
gen p = 2*normal(-abs(z))
qqvalue p, method(holm) qvalue(holm)
qqvalue p, method(bonferroni) qvalue(bonf)
qqvalue p, method(sidak) qvalue(sidak)
qqvalue p, method(holland) qvalue(holl)
qqvalue p, method(hochberg) qvalue(hochberg)
qqvalue p, method(simes) qvalue(simes)
qqvalue p, method(yekutieli) qvalue(yekutieli)
browse

*Create figure for Bonferroni p-values. Make P&I bar red and add text manually.
hist bonf, fcolor(ebblue) fintensity(25) lcolor(ebblue) lwidth(vvthin) width(.045) freq graphregion(color(white)) xtitle("Bonferroni-Corrected P-Value")  xlabel(0(.1)1) name(ucr39_hist_bonf, replace)
hist beta, fcolor(ebblue) fintensity(25) lcolor(ebblue) lwidth(vvthin) width(.015) freq graphregion(color(white)) xtitle("Coefficient Estimate")  name(ucr39_hist_beta, replace)
graph combine ucr39_hist_beta ucr39_hist_bonf, graphregion(color(white)) ycomm

*Create figure for Simes p-values. Make P&I bar red and add text manually.
hist simes, fcolor(ebblue) fintensity(25) lcolor(ebblue) lwidth(vvthin) width(.045) freq graphregion(color(white)) xtitle("Simes-Corrected P-Value")  xlabel(0(.1)1) name(ucr39_hist_simes, replace)
hist beta, fcolor(ebblue) fintensity(25) lcolor(ebblue) lwidth(vvthin) width(.015) freq graphregion(color(white)) xtitle("Coefficient Estimate")  name(ucr39_hist_beta, replace)
graph combine ucr39_hist_beta ucr39_hist_simes, graphregion(color(white)) ycomm




