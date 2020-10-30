*Analysis code for "Population Mortality and Laws Encouraging Influenza Vaccinsation for Hospital Workers"
*by Mariana Carrera, Emily Lawler and Corey White (equal authorship)
*Code last edited 10/29/2020


clear all

global PATH "SET FILE PATH HERE"

set matsize 2000

set seed 12345 










************************************************************************************************
************************************************************************************************
****************************** TABLES
************************************************************************************************
************************************************************************************************











************************************************************************************************
***************Table 1:   Main Result; heterogeneity across years since law/seasons/age
***************Table S2:  Raw log coefficients, differential effects, and bootstrapped p-values
************************************************************************************************


****************Main Result (for Tables 1, S2)

*Import analysis file
cd "$PATH\Influenza National\"
clear 
use flu_hcw_national_analysis

*Drop flu-years with match rate<50%
drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 

*Get pre-treatment mean
sum pi_rate if max_offer==1 & offer_post==0 [aweight=meanpop]
global mean_pi = r(mean)

*Covariates: indicator for long term care laws, indicator for childcare laws, state-specific linear trends
global controls offer_ltc req_child  i.state#c.ym

*Set panel variable (for fixed effects) to be state
xtset state

*regression for main result
xtreg ln_pi_rate offer_post  i.ym $controls [aweight=meanpop] , fe cluster(state)

*report regression coefficient/p-value
lincom offer_post

*Calculate change in levels: percent change (exp(estimate in logs)-1) times baseline mean. (Table 1)
*Report monthly change in mortality
disp ((exp(r(estimate))-1)*$mean_pi )
*Report lower bound 95% CI
disp ((exp(r(estimate) - (invttail(r(df),.025))*r(se))-1)*$mean_pi )
*Report upper bound 95% CI
disp ((exp(r(estimate) + (invttail(r(df),.025))*r(se))-1)*$mean_pi )

*Get Boostrapped  P-values -- state clustered, treatment stratified -- (Table S2)
*Regress out controls so prior to estimation to allow bootstrap to run
qui reg ln_pi_rate $controls i.ym i.state [aw=meanpop]
predict pi_resid, resid
qui reg offer_post $controls i.ym i.state [aw=meanpop]
predict offer_resid, resid
program regression
reg pi_resid offer_resid [aw=meanpop]
end
bootstrap _b, reps(200) cluster(state) strata(max_offer): regression 





****************Heterogeneity across years since law (Tables 1, S2)
*See above for additional notes (duplicate comments are not repeated)

*Import analysis file
cd "$PATH\Influenza National\"
clear 
use flu_hcw_national_analysis

drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 

sum pi_rate if max_offer==1 & offer_post==0 [aweight=meanpop]
global mean_pi = r(mean)

*Create a series of indicators for the number of years since law implementation.
forvalues i = 0/4{
gen year`i' = (yearstooffer==`i')
}

*Indicator for +4 captures 4 or more years since to implementation (included in regression so it is not part of the omitted group, but not displayed)
replace year4 = 1 if yearstooffer>=4 & !missing(yearstooffer)

global controls offer_ltc req_child  i.state#c.ym
xtset state

xtreg ln_pi_rate  year0-year4  i.ym $controls [aweight=meanpop] , fe cluster(state)

*report regression result: year0 (Table 1)
lincom year0
disp ((exp(r(estimate))-1)*$mean_pi )
disp ((exp(r(estimate) - (invttail(r(df),.025))*r(se))-1)*$mean_pi )
disp ((exp(r(estimate) + (invttail(r(df),.025))*r(se))-1)*$mean_pi )

*report regression result: year1 (Table 1)
lincom year1
disp ((exp(r(estimate))-1)*$mean_pi )
disp ((exp(r(estimate) - (invttail(r(df),.025))*r(se))-1)*$mean_pi )
disp ((exp(r(estimate) + (invttail(r(df),.025))*r(se))-1)*$mean_pi )

*report regression result: year2 (Table 1)
lincom year2
disp ((exp(r(estimate))-1)*$mean_pi )
disp ((exp(r(estimate) - (invttail(r(df),.025))*r(se))-1)*$mean_pi )
disp ((exp(r(estimate) + (invttail(r(df),.025))*r(se))-1)*$mean_pi )

*report regression result: year3 (Table 1)
lincom year3
disp ((exp(r(estimate))-1)*$mean_pi )
disp ((exp(r(estimate) - (invttail(r(df),.025))*r(se))-1)*$mean_pi )
disp ((exp(r(estimate) + (invttail(r(df),.025))*r(se))-1)*$mean_pi )

*Bootstrapped p-values for each variable (Table S2)
qui reg ln_pi_rate $controls i.ym i.state [aw=meanpop]
predict pi_resid, resid
qui reg year0 $controls i.ym i.state [aw=meanpop]
predict year0_resid, resid
qui reg year1 $controls i.ym i.state [aw=meanpop]
predict year1_resid, resid
qui reg year2 $controls i.ym i.state [aw=meanpop]
predict year2_resid, resid
qui reg year3 $controls i.ym i.state [aw=meanpop]
predict year3_resid, resid
qui reg year4 $controls i.ym i.state [aw=meanpop]
predict year4_resid, resid
program regression_dynamics
reg pi_resid year0_resid year1_resid year2_resid year3_resid year4_resid [aw=meanpop]
end
bootstrap _b, reps(200) cluster(state) strata(max_offer): regression_dynamics 







****************Heterogeneity across seasons (Tables 1, S2)
*See above for additional notes (duplicate comments are not repeated)

cd "$PATH\Influenza National\"
clear 
use flu_hcw_national_analysis

*Generate interaction for law and peak season
gen offer_season = offer_post*season

drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 

*Get pre-treatment means
*peak months
sum pi_rate if max_offer==1 & offer_post==0 & season==1 [aweight=meanpop]
global mean_pi_season = r(mean)

*non-peak months
sum pi_rate if max_offer==1 & offer_post==0 & season==0 [aweight=meanpop]
global mean_pi_nonseason = r(mean)

*"sf_group" is an identifier for stateXseason groups. With the interaction model, we include stateXseason FE and stateXseason trends
global controls offer_ltc req_child  i.sf_group#c.ym
xtset sf_group

*Regression for season interaction model
xtreg ln_pi_rate offer_post offer_season i.ym $controls [aweight=meanpop] , fe cluster(state)

*Results: Effect of laws during non-peak season (Table 1)
lincom offer_post
disp ((exp(r(estimate))-1)*$mean_pi_nonseason )
disp ((exp(r(estimate) - (invttail(r(df),.025))*r(se))-1)*$mean_pi_nonseason )
disp ((exp(r(estimate) + (invttail(r(df),.025))*r(se))-1)*$mean_pi_nonseason )

*Results: Effect of laws during peak season (Table 1)
lincom offer_post + offer_season
disp ((exp(r(estimate))-1)*$mean_pi_season )
disp ((exp(r(estimate) - (invttail(r(df),.025))*r(se))-1)*$mean_pi_season )
disp ((exp(r(estimate) + (invttail(r(df),.025))*r(se))-1)*$mean_pi_season )

*Results: Differential effect for peak vs. non-peak (Table S2)
lincom offer_season

*Bootstrapped  P-values (Table S2)
qui reg ln_pi_rate $controls i.ym i.sf_group [aw=meanpop]
predict pi_resid, resid
qui reg offer_post $controls i.ym i.sf_group [aw=meanpop]
predict offer_resid, resid
qui reg offer_season $controls i.ym i.sf_group [aw=meanpop]
predict int_resid, resid
program regression_season
reg pi_resid offer_resid int_resid [aw=meanpop]
end
bootstrap (_b[offer_resid]) (_b[offer_resid]+_b[int_resid]) (_b[int_resid]) , reps(200) cluster(state) strata(max_offer): regression_season 





****************Heterogeneity across age (Tables 1, S2)
*See above for additional notes (duplicate comments are not repeated)

cd "$PATH\Influenza National\"
clear 
use flu_hcw_national_analysis

gen ln_pi_rate_age_065 = ln(pi_rate_age_065)
gen ln_pi_rate_age_65 = ln(pi_rate_age_65)

*To run an interacted model with age, we need to reshape the data to long format so that each observation is state-year-month-age; the next lines of code do this reshape.
keep ln_pi_rate_age_065 ln_pi_rate_age_65 pi_rate_age_065 pi_rate_age_65 max_offer yearstooffer season meanpop temp* hum* prcp* popshare* unemp pov offer_ltc req_child sf ym state flu_year offer_post
rename ln_pi_rate_age_065 ln_pi_rate1
rename ln_pi_rate_age_65 ln_pi_rate2
rename pi_rate_age_065 pi_rate1
rename pi_rate_age_65 pi_rate2
reshape long ln_pi_rate pi_rate, i(state ym) j(age)

drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 

*Generate interaction with 65+
gen o65 = (age==2)
gen offer_age = offer_post*o65

*Generate new variables for stateXage and timeXage FE
egen sa_group = group(state age)
egen yma_group = group(ym age)

*Note controls include the stateXage time trends
global controls offer_ltc req_child i.sa_group#c.ym
xtset sa_group

*Calculate means for each group
sum pi_rate if max_offer==1 & yearstooffer<0 & o65==1 [aweight=meanpop]
global mean_pi_o65 = r(mean)
sum pi_rate if max_offer==1 & yearstooffer<0 & o65==0 [aweight=meanpop]
global mean_pi_u65 = r(mean)

*Age heterogeneity regression
xtreg ln_pi_rate  offer_post offer_age i.yma $controls [aweight=meanpop] , fe cluster(state)

*Results: effect of laws for <65 (Table 1)
lincom offer_post
disp ((exp(r(estimate))-1)*$mean_pi_u65 )
disp ((exp(r(estimate) - (invttail(r(df),.025))*r(se))-1)*$mean_pi_u65 )
disp ((exp(r(estimate) + (invttail(r(df),.025))*r(se))-1)*$mean_pi_u65 )

*Results: effect of laws for >65 (Table 1)
lincom offer_post + offer_age
disp ((exp(r(estimate))-1)*$mean_pi_o65 )
disp ((exp(r(estimate) - (invttail(r(df),.025))*r(se))-1)*$mean_pi_o65 )
disp ((exp(r(estimate) + (invttail(r(df),.025))*r(se))-1)*$mean_pi_o65 )

*Results: Differential effect of laws for >65 vs. <65 (Table S2)
lincom offer_age

*Bootstrapped P-values (Table S2)
qui reg ln_pi_rate $controls i.yma i.sa_group [aw=meanpop]
predict pi_resid, resid
qui reg offer_post $controls i.yma i.sa_group [aw=meanpop]
predict offer_resid, resid
qui reg offer_age $controls i.yma i.sa_group [aw=meanpop]
predict int_resid, resid
program regression_age
reg pi_resid offer_resid int_resid [aw=meanpop]
end
bootstrap (_b[offer_resid]) (_b[offer_resid]+_b[int_resid]) (_b[int_resid]) , reps(200) cluster(state) strata(max_offer): regression_age











************************************************************************************************
***************Sensitivity Tests: Table 2 and Figure S4
************************************************************************************************

************************Varying Fixed Effects and Controls (Panel A)

cd "$PATH\Influenza National\"

clear
use flu_hcw_national_analysis

*Define covariates (not including trends, those are added below in some specifications)
global covariates offer_ltc req_child 

drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014

*excludes covariates, trends
xtset state
xtreg ln_pi_rate offer_post i.ym  [aweight=meanpop], fe cluster(state)
lincom offer_post

disp ((exp(r(estimate))-1)*$mean_pi )
disp ((exp(r(estimate) - (invttail(r(df),.025))*r(se))-1)*$mean_pi )
disp ((exp(r(estimate) + (invttail(r(df),.025))*r(se))-1)*$mean_pi )

*excludes trends
xtreg ln_pi_rate offer_post $covariates  i.ym  [aweight=meanpop], fe cluster(state)
lincom offer_post

disp ((exp(r(estimate))-1)*$mean_pi )
disp ((exp(r(estimate) - (invttail(r(df),.025))*r(se))-1)*$mean_pi )
disp ((exp(r(estimate) + (invttail(r(df),.025))*r(se))-1)*$mean_pi )

*excludes covariates
xtreg ln_pi_rate offer_post i.state#c.ym  i.ym  [aweight=meanpop], fe cluster(state)
lincom offer_post

disp ((exp(r(estimate))-1)*$mean_pi )
disp ((exp(r(estimate) - (invttail(r(df),.025))*r(se))-1)*$mean_pi )
disp ((exp(r(estimate) + (invttail(r(df),.025))*r(se))-1)*$mean_pi )


************************Varying Years (Panel B)

cd "$PATH\Influenza National\"

clear
use flu_hcw_national_analysis

global controls offer_ltc req_child i.state#c.ym

*Include bad match years
xtset state
xtreg ln_pi_rate offer_post  $controls i.ym  [aweight=meanpop], fe cluster(state)
lincom offer_post
disp ((exp(r(estimate))-1)*$mean_pi )
disp ((exp(r(estimate) - (invttail(r(df),.025))*r(se))-1)*$mean_pi )
disp ((exp(r(estimate) + (invttail(r(df),.025))*r(se))-1)*$mean_pi )

*Drop bad match years  for remaining regressions
drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014

*Drop H1N1 years
xtreg ln_pi_rate offer_post  $controls i.ym if flu_year!=2008 & flu_year != 2009   [aweight=meanpop], fe cluster(state)
lincom offer_post
disp ((exp(r(estimate))-1)*$mean_pi )
disp ((exp(r(estimate) - (invttail(r(df),.025))*r(se))-1)*$mean_pi )
disp ((exp(r(estimate) + (invttail(r(df),.025))*r(se))-1)*$mean_pi )


************************Poisson (Panel C)

cd "$PATH\Influenza National\"
clear 
use flu_hcw_national_analysis
drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 
sum pi_rate if max_offer==1 & offer_post==0 [aweight=meanpop]
global mean_pi = r(mean)
global controls offer_ltc req_child  i.state#c.ym
xtset state

*Poisson regression
xtpoisson pi offer_post  i.ym $controls  , exposure(pop_tot) fe vce(robust)

*Use model df
lincom offer_post
disp ((exp(r(estimate))-1)*$mean_pi )
disp ((exp(r(estimate) - (invttail(e(df_m),.025))*r(se))-1)*$mean_pi )
disp ((exp(r(estimate) + (invttail(e(df_m),.025))*r(se))-1)*$mean_pi )







****************Annual-Level Estimates + DIDm Estimator

cd "$PATH\Influenza National\"

clear
use flu_hcw_national_analysis

*Collapse to annual level
collapse (sum) pi pop_tot (mean) offer_ltc req_child offer_post max_offer offeryear meanpop, by(state flu_year)

*Re-generate P&I mortality rate
gen pi_rate = 100000*pi/pop_tot
gen ln_pi_rate=ln(pi_rate)

drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 

*Generate trends variables
tab state, gen(ST_)
forvalues i=1/51{
replace ST_`i' = ST_`i'*flu_year
} 

*Estimate standard 2fe regression (restricting to post-treatment years 0-4) to get degrees of freedom
xtset state
xtreg ln_pi_rate offer_post i.flu_year ST_* offer_ltc req_child  [aw=meanpop], fe cluster(state)
lincom offer_post
disp ((exp(r(estimate))-1)*$mean_pi )
disp ((exp(r(estimate) - (invttail(r(df),.025))*r(se))-1)*$mean_pi )
disp ((exp(r(estimate) + (invttail(r(df),.025))*r(se))-1)*$mean_pi )
global dfreedom = r(df)

xtset state
did_multiplegt ln_pi_rate state flu_year offer_post, controls(ST_* offer_ltc req_child) dynamic(3) breps(200) weight(meanpop) average_effect(simple) covariances 

*Display raw estimate, SE, t-stat, p-value
disp e(effect_average)
disp e(se_effect_average)
disp e(effect_average)/e(se_effect_average)
disp ttail($dfreedom ,abs(e(effect_average)/e(se_effect_average)))*2

*Report estimate/CI in levels
disp ((exp(e(effect_average))-1)*$mean_pi )
disp ((exp(e(effect_average) - (invttail($dfreedom ,.025))*e(se_effect_average))-1)*$mean_pi )
disp ((exp(e(effect_average) + (invttail($dfreedom ,.025))*e(se_effect_average))-1)*$mean_pi )

*Calculate sum of weights
twowayfeweights ln_pi_rate state flu_year offer_post, type(feTR) 

*Dynamic plot
did_multiplegt ln_pi_rate state flu_year offer_post, controls(ST_* offer_ltc req_child) placebo(7) dynamic(3) breps(200) weight(meanpop) average_effect(simple) covariances  

matrix DIDM = J(11,3,.)
matrix DIDM[1,1] = (exp(e(placebo_7))-1)*$mean_pi , (exp(e(placebo_7) - ((invttail($dfreedom ,.025))*e(se_placebo_7)))-1)*$mean_pi ,(exp(e(placebo_7) + ((invttail($dfreedom ,.025))*e(se_placebo_7)))-1)*$mean_pi
matrix DIDM[2,1] = (exp(e(placebo_6))-1)*$mean_pi , (exp(e(placebo_6) - ((invttail($dfreedom ,.025))*e(se_placebo_6)))-1)*$mean_pi ,(exp(e(placebo_6) + ((invttail($dfreedom ,.025))*e(se_placebo_6)))-1)*$mean_pi
matrix DIDM[3,1] = (exp(e(placebo_5))-1)*$mean_pi , (exp(e(placebo_5) - ((invttail($dfreedom ,.025))*e(se_placebo_5)))-1)*$mean_pi ,(exp(e(placebo_5) + ((invttail($dfreedom ,.025))*e(se_placebo_5)))-1)*$mean_pi
matrix DIDM[4,1] = (exp(e(placebo_4))-1)*$mean_pi , (exp(e(placebo_4) - ((invttail($dfreedom ,.025))*e(se_placebo_4)))-1)*$mean_pi ,(exp(e(placebo_4) + ((invttail($dfreedom ,.025))*e(se_placebo_4)))-1)*$mean_pi
matrix DIDM[5,1] = (exp(e(placebo_3))-1)*$mean_pi , (exp(e(placebo_3) - ((invttail($dfreedom ,.025))*e(se_placebo_3)))-1)*$mean_pi ,(exp(e(placebo_3) + ((invttail($dfreedom ,.025))*e(se_placebo_3)))-1)*$mean_pi
matrix DIDM[6,1] = (exp(e(placebo_2))-1)*$mean_pi , (exp(e(placebo_2) - ((invttail($dfreedom ,.025))*e(se_placebo_2)))-1)*$mean_pi ,(exp(e(placebo_2) + ((invttail($dfreedom ,.025))*e(se_placebo_2)))-1)*$mean_pi
matrix DIDM[7,1] = 0, 0, 0
matrix DIDM[8,1] = (exp(e(effect_0))-1)*$mean_pi , (exp(e(effect_0) - ((invttail($dfreedom ,.025))*e(se_effect_0)))-1)*$mean_pi , (exp(e(effect_0) +((invttail($dfreedom ,.025))*e(se_effect_0)))-1)*$mean_pi
matrix DIDM[9,1] = (exp(e(effect_1))-1)*$mean_pi , (exp(e(effect_1) - ((invttail($dfreedom ,.025))*e(se_effect_1)))-1)*$mean_pi , (exp(e(effect_1) + ((invttail($dfreedom ,.025))*e(se_effect_1)))-1)*$mean_pi
matrix DIDM[10,1] = (exp(e(effect_2))-1)*$mean_pi , (exp(e(effect_2) - ((invttail($dfreedom ,.025))*e(se_effect_2)))-1)*$mean_pi , (exp(e(effect_2) + ((invttail($dfreedom ,.025))*e(se_effect_2)))-1)*$mean_pi
matrix DIDM[11,1] = (exp(e(effect_3))-1)*$mean_pi , (exp(e(effect_3) - ((invttail($dfreedom ,.025))*e(se_effect_3)))-1)*$mean_pi , (exp(e(effect_3) + ((invttail($dfreedom ,.025))*e(se_effect_3)))-1)*$mean_pi

*Make All Month, and Peak/Non-Peak plots separately

 coefplot (matrix(DIDM[,1]),  ci((DIDM[,2] DIDM[,3])) lcolor(black) mcolor(black) mlcolor(black) msymbol(D) ylabel(, angle(0)) ///
 mlwidth(thin) connect(direct) ciopts(recast(rcap) lcolor(black))) ,  xlabel(1 "-7" 2 "-6" 3 "-5" 4 "-4" 5 "-3" 6 "-2" 7 "-1" 8 "0" 9 "1" 10 "2" 11 "3" ) ///
 leg(off)  xtitle("Influenza-Years Relative to Law") ytitle("Change in Monthly P&I Deaths per 100,000 Population") graphregion(color(white)) vertical  yline(0, lcolor(black) lwidth(medthin)) ///
 xline(8, lcolor(gs10) lwidth(thin)) name(didm_dynamic, replace)







************************************************************************************************
***************Table S3: Effects of Laws on Adolescent and Infant Vaccination Rates (NIS Data)
************************************************************************************************


*************** NIS Teen
*global root "DIRECTORY FOR NIS DATA"

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
global controls "offer_ltc req_child" //need to also add unemp and pov, maybe also MS & HS vaccine mandates?

*Excluding bad match years, with state and influenza year FEs, including time varying controls	ADDING LINEAR STATE TRENDS
	logit flu_pastyr  offer_post $controls i.flu_year i.state $linear if bad_match==0 [pw=provwt],  vce(cluster state)
margins , dydx(offer_post) atmeans





*************** NIS Child
*****************************************************
	use "$data\nis_analysisfile.dta", clear

*generate state-specific linear time trends at the year level
gen trend = flu_year- 2001
	
	xi i.state*trend 	
	global linear "_IstaXtren_*"

*define vector of time varying controls
global controls "offer_ltc req_child" 

*keep only the youngest age group
keep if agegrp==1

*************** regressions ***********************

*first do estimation for influenza vaccine, then for UTD PCV 
foreach var of varlist flu_any P_UTDPCV {

*Excluding bad match years, with state and influenza year FEs, including time varying controls	ADDING LINEAR STATE TRENDS
    logit `var' offer_post $controls   i.flu_year i.state $linear if bad_match==0  [pw=provwt],  vce(cluster state)	
	margins , dydx(offer_post) atmeans
	}
	

	*get baseline mean from ever adopters
	eststo: reg flu_any offer_post $controls   i.flu_year i.state $linear if bad_match==0  & agegrp==1 [aw=provwt],  vce(cluster state)
	
	sum flu_any if e(sample)==1 & yearstooffer<0 [aw=provwt] 
	
	eststo: reg P_UTDPCV offer_post $controls   i.flu_year i.state $linear if bad_match==0  & agegrp==1 [aw=provwt],  vce(cluster state)
	
	sum P_UTDPCV if e(sample)==1 & yearstooffer<0 [aw=provwt] 
	



	
	
	

************************************************************************************************
***************Table S4: Effects of Laws on Adult Vaccination Rates (BRFSS Data)
************************************************************************************************


cd "BRFSS FILE PATH"

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


global controls offer_ltc req_child 
xtset state


/*
foreach var in flu_vaxrate flu_u65 flu_o65 pvac pvac_u6 pvac_o65{
replace `var' = ln(`var')
}
*/

*Influenza vaccination rate regressions
foreach outcome in flu_vaxrate flu_u65 flu_o65 pvac pvac_u65 pvac_o65{
xtreg `outcome' offer_post  i.flu_year $controls [aweight=meanpop] , fe cluster(state)
sum `outcome' if max_offer==1 & offer_post==0 [aw=meanpop] 
local mean_outcome = r(mean)
lincom offer_post
*Report estimate
disp ((exp(r(estimate))-1)*`mean_outcome' )
*Report lower bound 95% CI
disp ((exp(r(estimate) - (invttail(r(df),.025))*r(se))-1)*`mean_outcome' )
*Report upper bound 95% CI
disp ((exp(r(estimate) + (invttail(r(df),.025))*r(se))-1)*`mean_outcome' )
}













************************************************************************************************
***************Table S5: E-Values
************************************************************************************************


cd "$PATH\Influenza National\"
clear 
use flu_hcw_national_analysis

*****Monthly Data (Row 1-2)
drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 
sum pi_rate if max_offer==1 & offer_post==0 [aweight=meanpop]
global mean_pi = r(mean)
global controls offer_ltc req_child  i.state#c.ym
xtset state
xtreg ln_pi_rate offer_post  i.ym $controls [aweight=meanpop] , fe cluster(state)

*Get treatment effect and SE
lincom offer_post
local effect = r(estimate)
local ste = r(se)

*get SD
sum ln_pi_rate [aw=meanpop]
local sd = r(sd)

*Get residualized SD
reg ln_pi_rate i.state i.flu_year $controls [aw=meanpop]
predict pi_res, resid
sum pi_res [aw=meanpop]
local sd_adj = r(sd)

*Calculate e-values and CIs (unadjusted)
local d = `effect'/`sd'
local RR = exp((0.91*`d'))
local RR_l = exp((0.91*`d'-1.78*`ste'))
local RR_u = exp((0.91*`d'+1.78*`ste'))

disp (1/`RR')+sqrt((1/`RR')*((1/`RR')-1))
disp (1/`RR_l')+sqrt((1/`RR_l')*((1/`RR_l')-1))
disp (1/`RR_u')+sqrt((1/`RR_u')*((1/`RR_u')-1))

*Calculate e-values and CIs (adjusted)
local d = `effect'/`sd_adj'
local RR = exp((0.91*`d'))
local RR_l = exp((0.91*`d'-1.78*`ste'))
local RR_u = exp((0.91*`d'+1.78*`ste'))

disp (1/`RR')+sqrt((1/`RR')*((1/`RR')-1))
disp (1/`RR_l')+sqrt((1/`RR_l')*((1/`RR_l')-1))
disp (1/`RR_u')+sqrt((1/`RR_u')*((1/`RR_u')-1))


******Annual Data (Row 3-4) -- same process as above; see above for comments
clear
use flu_hcw_national_analysis
collapse (sum) pi pop_tot (mean) offer_ltc req_child offer_post max_offer offeryear meanpop, by(state flu_year)
gen pi_rate = 100000*pi/pop_tot
gen ln_pi_rate=ln(pi_rate)
drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 
global controls offer_ltc req_child  i.state#c.flu_year
xtset state
xtreg ln_pi_rate offer_post i.flu_year $controls [aw=meanpop], fe cluster(state)

reg ln_pi_rate offer_post i.state i.flu_year $controls [aw=meanpop]

predict resid, resid

lincom offer_post
local effect = r(estimate)
local ste = r(se)

sum ln_pi_rate [aw=meanpop]
local sd = r(sd)

reg ln_pi_rate i.state i.flu_year $controls [aw=meanpop]
predict pi_res, resid
sum pi_res [aw=meanpop]
local sd_adj = r(sd)

local d = `effect'/`sd'
local RR = exp((0.91*`d'))
local RR_l = exp((0.91*`d'-1.78*`ste'))
local RR_u = exp((0.91*`d'+1.78*`ste'))

disp (1/`RR')+sqrt((1/`RR')*((1/`RR')-1))
disp (1/`RR_l')+sqrt((1/`RR_l')*((1/`RR_l')-1))
disp (1/`RR_u')+sqrt((1/`RR_u')*((1/`RR_u')-1))

local d = `effect'/`sd_adj'
local RR = exp((0.91*`d'))
local RR_l = exp((0.91*`d'-1.78*`ste'))
local RR_u = exp((0.91*`d'+1.78*`ste'))

disp (1/`RR')+sqrt((1/`RR')*((1/`RR')-1))
disp (1/`RR_l')+sqrt((1/`RR_l')*((1/`RR_l')-1))
disp (1/`RR_u')+sqrt((1/`RR_u')*((1/`RR_u')-1))








**************************************************************************
****************Table S6: Estimates by 39 causes of death. 
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
postfile mtpvalues me beta ste str20 estm using mtpvalues.dta, replace

*Loop over all 39 causes of death. Note that Cause 24 is P&I. Cause 2 is syphilis, and the model does not converge for that outcome.
local i =1
foreach var in 24 1 /*2*/ 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39{

*Note the UCR_PI_24 variable is the number of deaths with any primary/secondary diagnosis for P&I
*UCR_PI_* for other categories represents deaths with another primary category and no secondary diagosis for P&I

gen UCR_PI_rate_`var' = 100000*UCR_PI_`var'/pop_tot

sum UCR_PI_rate_`var' if offer_post==0 & max_offer==1 [aw=meanpop], meanonly
local UCR_mean = r(mean)

global controls i.state#c.ym

xtset state

*Run regression, for peak months only where largest effect for P&I is expected
 xtpoisson UCR_PI_`var' offer_post  i.ym $controls if season==1 , exposure(pop_tot) fe vce(robust)
lincom offer_post
local est1 = r(estimate)
local se1 = r(se)

post mtpvalues ((exp(`est1')-1)*`UCR_mean') (`est1') (`se1') ("est_`var'")

local i = `i'+1
}

postclose mtpvalues

*Open file with estimates, create p-values for each estimate. Create adjusted p-values correcting for multiple hypothesis testing
*Use a variety of adjustments. Note that bonferroni is the most conservative and simes is the most liberal among these.
clear 
use mtpvalues
gen z = beta/ste
gen p = 2*normal(-abs(z))
qqvalue p, method(bonferroni) qvalue(bonf)
qqvalue p, method(simes) qvalue(simes)

*Following Altman and Bland BMJ 2011 method for getting CI from p-value and estimate
gen z_bonf = -0.862 + sqrt(0.743 - 2.404*ln(bonf))
gen se_bonf = abs(me/z_bonf)
gen lb_bonf = me - 1.96*se_bonf 
gen ub_bonf = me + 1.96*se_bonf 

*Following Altman and Bland BMJ 2011 method for getting CI from p-value and estimate
gen z_simes = -0.862 + sqrt(0.743 - 2.404*ln(simes))
gen se_simes = abs(me/z_simes)
gen lb_simes = me - 1.96*se_simes 
gen ub_simes = me + 1.96*se_simes 

*Create table
sort simes
split estm, parse("_")
destring estm2, replace
foreach var in me lb_simes ub_simes bonf simes{
replace `var' = round(`var', 0.001)
}
local i =1
matrix C = J(38,6,.)
forvalues i=1/38{
matrix C[`i',1] = estm2[`i'], me[`i'], lb_simes[`i'], ub_simes[`i'], bonf[`i'], simes[`i']
}
mat colnames C = "Cause" "ME" "LB-95" "UB-95" "Bonferroni" "Simes"
esttab matrix(C)
esttab matrix(C) using "$PATH\Influenza HCW National\Draft\Export Files\cause_table.csv", replace





















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
 ytitle("Share of U.S. Population", axis(2)) legend(rows(3))  ylabel(0 "0" .2 "0.2" .4 "0.4" .6 "0.6" .8 "0.8" ,angle(0) axis(2)) ///
 yaxis(2) fcolor(ebblue) fintensity(25) lcolor(ebblue*.25) graphregion(color(white)) ///
|| line pi_rate flu_year if max_offer==1,  lcolor(black) yaxis(1) ytitle("Monthly P&I Deaths per 100,000 Population", axis(1)) ///
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
 ytitle("Share of U.S. Hospital Workers")  ylabel(0 .2 "0.2" .4 "0.4" .6 "0.6" .8 "0.8", angle(horizontal)) ///
 fcolor(ebblue) fintensity(25) lcolor(ebblue*.25) graphregion(color(white)) ///
|| line vacflush12m  flu_year if ind_cat==1,  lcolor(black)  ///
|| line vacflush12m flu_year if ind_cat==0 , lcolor(gs10) lpattern(shortdash)  ///
 legend(cols(1) order( 1 "Share of Hospital Workers Affected by Laws" 2 "Share Vaccinated: Hospital Workers" 3 "Share Vaccinated: Employed, Non-Healthcare")) ///
 xtitle("Influenza Year")


***************Combine Figure 1A and 1B
graph combine vax_trends mort_trends, xsize(4) ysize(2) graphregion(color(white))







  
  
**************************************************************************
*************** Figures 2, 3, and S2: Synthetic Control Trends by State, Synth Treatment Effect Graphs, Pre-treatmeant RMSPE
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

*Collapse to annual level (collapse population as a sum to get monthly death rates)
collapse (sum) pi pop_tot offer_post  (mean) meanpop max_offer, by(state flu_year)

*Re-generate offeryear variable
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
gen pi_rate = (100000*pi/pop_tot)
gen ln_pi_rate = ln(pi_rate)

*Set panel and time variables (for synth command)
tsset state flu_year

*Get pre-treatment mean to construct marginal effect
sum pi_rate if max_offer==1 & offer_post==0 [aw=meanpop]
global mean_pi = r(mean)

*Set up file to store regression estimates for each state
cap postclose synth_estimates
postfile synth_estimates me pval pval_rmspe rmspe_pre beta state year treated str20 stateabr using synth_estimates.dta, replace

*set up file to store placebo estimates (for p-values)
cap postclose synth_placebo
postfile synth_placebo beta rmspe_ratio state year treated using synth_placebo.dta, replace


*Run placebo estimates for all treatment years/donor states (for p-values)
foreach state in 2 4 5 9 10 12 15 16 18 19 20 21 22 26 27 28 29 30 32 34 35 37 38 39 41 42 45 46 48 49 50 51 53 54 55 56{
preserve
keep if (state==`state' | max_offer==0)
qui{
*2002
synth ln_pi_rate ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995) , trunit(`state') trperiod(2002) keep(synth_2002_`state', replace)
*2005
synth ln_pi_rate ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995) , trunit(`state') trperiod(2005) keep(synth_2005_`state', replace)
*2007
synth ln_pi_rate ln_pi_rate(2006) ln_pi_rate(2005) ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995), trunit(`state') trperiod(2007) keep(synth_2007_`state', replace)
*2008
synth ln_pi_rate ln_pi_rate(2007) ln_pi_rate(2006) ln_pi_rate(2005) ///
ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995)  , trunit(`state') trperiod(2008) keep(synth_2008_`state', replace)
*2009
synth ln_pi_rate ln_pi_rate(2008) ln_pi_rate(2007) ln_pi_rate(2006) ln_pi_rate(2005) ///
ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995) , trunit(`state') trperiod(2009) keep(synth_2009_`state', replace)
*2010
synth ln_pi_rate ln_pi_rate(2009) ln_pi_rate(2008) ln_pi_rate(2007) ln_pi_rate(2006) ln_pi_rate(2005) ///
ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995) , trunit(`state') trperiod(2010) keep(synth_2010_`state', replace)
*2011
synth ln_pi_rate  ln_pi_rate(2010) ln_pi_rate(2009) ln_pi_rate(2008) ln_pi_rate(2007) ln_pi_rate(2006) ln_pi_rate(2005) ///
ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995) , trunit(`state') trperiod(2011) keep(synth_2011_`state', replace)
*2012
synth ln_pi_rate ln_pi_rate(2011) ln_pi_rate(2010) ln_pi_rate(2009) ln_pi_rate(2008) ln_pi_rate(2007) ln_pi_rate(2006) ln_pi_rate(2005) ///
ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995) , trunit(`state') trperiod(2012) keep(synth_2012_`state', replace)
*2013
synth ln_pi_rate ln_pi_rate(2012) ln_pi_rate(2011) ln_pi_rate(2010) ln_pi_rate(2009) ln_pi_rate(2008) ln_pi_rate(2007) ln_pi_rate(2006) ln_pi_rate(2005) ///
ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995) , trunit(`state') trperiod(2013) keep(synth_2013_`state', replace)
}
restore
}


*Store results from placebo estimates (including the post/pre RMSPE ratio)
preserve
foreach trs in 2 4 5 9 10 12 15 16 18 19 20 21 22 26 27 28 29 30 32 34 35 37 38 39 41 42 45 46 48 49 50 51 53 54 55 56{
foreach trp in 2002 2005 2007 2008 2009 2010 2011 2012 2013{
qui{
*Calculate post/pre RMSPE ratio for Abadie (2020) inference method
clear
use synth_`trp'_`trs'
gen post_treat = 0 if _time<`trp'
replace post_treat = 1 if _time>=`trp' &  _time!=.
bysort post_treat: gen num = _N if _time!=.
gen diffsq = (_Y_treated- _Y_synthetic)^2
bysort post_treat: egen sumdiffsq = sum(diffsq) if _time!=.
gen rmspe_pre_temp = sqrt((1/num)*sumdiffsq) if post_treat==0
gen rmspe_post_temp = sqrt((1/num)*sumdiffsq) if post_treat==1
egen rmspe_pre = min(rmspe_pre_temp)
egen rmspe_post = min(rmspe_post_temp)
local rmspe_ratio = rmspe_post[1]/rmspe_pre[1]
*Get betas using monthly data
clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_`trp'_`trs'
drop _m
replace _W_Weight = 1 if state==`trs'
*assign placebo treatment
replace max_offer = 1 if state==`trs'
reg ln_pi_rate max_offer  if flu_year>=`trp'  [aweight=_W_Weight]
*store betas, RMSPE ratio, and state/time identifiers
post synth_placebo (_b[max_offer]) (`rmspe_ratio') (`trs') (`trp') (0)
}
}
}
postclose synth_placebo
restore



*Each state is treated separately, CA is first and is commented, the other states follow the same outline are are not commented
*CA
preserve
local trs=6
local trp=2007

*Keep only CA and non-adopting states
keep if (state==`trs' | max_offer==0)

*Run synth command, match on the outcome in each pre-treatment year (not average of the outcome over all pre-treatment years).
synth ln_pi_rate ln_pi_rate(2006) ln_pi_rate(2005) ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995), trunit(`trs') trperiod(`trp') keep(synth_`trs', replace)

*Pull in main data to run synthetic control regressions (i.e., using synth weights)
clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_`trs' 
drop _m
replace _W_Weight = 1 if state==`trs'
*Simple regression
reg ln_pi_rate max_offer  if flu_year>=`trp'  [aweight=_W_Weight]

*Collapse data and create figure
gen treatment = (state==`trs')
collapse pi_rate [aw=_W], by(flu_year treatment)
twoway line pi_rate flu_year if treatment==1, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white))  ylab(3 6 9, angle(0) nogrid) /// 
xtitle("") ytitle("") title("California (2007/08)", color(black)) name(CA, replace) ///
||  line pi_rate flu_year  if treatment==0, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(`trp',  lcolor(gs10) lwidth(medthin))

*Get RMSPE ratio for inference
clear
use synth_`trs'
gen post_treat = 0 if _time<`trp'
replace post_treat = 1 if _time>=`trp' &  _time!=.
bysort post_treat: gen num = _N if _time!=.
gen diffsq = (_Y_treated- _Y_synthetic)^2
bysort post_treat: egen sumdiffsq = sum(diffsq) if _time!=.
gen rmspe_pre_temp = sqrt((1/num)*sumdiffsq) if post_treat==0
gen rmspe_post_temp = sqrt((1/num)*sumdiffsq) if post_treat==1
egen rmspe_pre = min(rmspe_pre_temp)
egen rmspe_post = min(rmspe_post_temp)
local rmspe_ratio = rmspe_post[1]/rmspe_pre[1]
local rmspe_pre = rmspe_pre[1]

*Pull in placebo data and calculate p-values using the distribution of RMSPE ratios
clear
use synth_placebo
local abs_beta = abs(_b[max_offer])
disp `abs_beta'
count if abs(beta)>`abs_beta'
local pval = r(N)/_N
disp `pval'
count if rmspe_ratio>`rmspe_ratio'
local pval_rmspe = r(N)/_N
disp `pval_rmspe'

*Post output to file
post synth_estimates ((exp(_b[max_offer])-1)*$mean_pi ) (`pval') (`pval_rmspe') (`rmspe_pre') (_b[max_offer]) (`trs') (`trp') (1) ("CA")
*restore data and run for the remaining states
restore 

*CO
local trs=8
local trp=2012
preserve
keep if (state==`trs' | max_offer==0)
synth ln_pi_rate ln_pi_rate(2011) ln_pi_rate(2010) ln_pi_rate(2009) ln_pi_rate(2008) ln_pi_rate(2007) ln_pi_rate(2006) ln_pi_rate(2005) ///
ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995) , trunit(`trs') trperiod(`trp') keep(synth_`trs', replace)
clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_`trs' 
drop _m
replace _W_Weight = 1 if state==`trs'
reg ln_pi_rate max_offer  if flu_year>=`trp'  [aweight=_W_Weight]
gen treatment = (state==`trs')
collapse pi_rate [aw=_W], by(flu_year treatment)
twoway line pi_rate flu_year if treatment==1, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) ylab(3 6 9, angle(0) nogrid) /// 
xtitle("") ytitle("")  title("Colorado (2012/13)", color(black)) name(CO, replace) ///
||  line pi_rate flu_year if treatment==0, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(`trp',  lcolor(gs10) lwidth(medthin))
clear
use synth_`trs'
gen post_treat = 0 if _time<`trp'
replace post_treat = 1 if _time>=`trp' &  _time!=.
bysort post_treat: gen num = _N if _time!=.
gen diffsq = (_Y_treated- _Y_synthetic)^2
bysort post_treat: egen sumdiffsq = sum(diffsq) if _time!=.
gen rmspe_pre_temp = sqrt((1/num)*sumdiffsq) if post_treat==0
gen rmspe_post_temp = sqrt((1/num)*sumdiffsq) if post_treat==1
egen rmspe_pre = min(rmspe_pre_temp)
egen rmspe_post = min(rmspe_post_temp)
local rmspe_ratio = rmspe_post[1]/rmspe_pre[1]
local rmspe_pre = rmspe_pre[1]
clear
use synth_placebo
local abs_beta = abs(_b[max_offer])
disp `abs_beta'
count if abs(beta)>`abs_beta'
local pval = r(N)/_N
disp `pval'
count if rmspe_ratio>`rmspe_ratio'
local pval_rmspe = r(N)/_N
disp `pval_rmspe'
post synth_estimates ((exp(_b[max_offer])-1)*$mean_pi ) (`pval') (`pval_rmspe') (`rmspe_pre') (_b[max_offer]) (`trs') (`trp') (1) ("CO")
restore 

*DC
local trs=11
local trp=2008
preserve
keep if (state==`trs' | max_offer==0)
synth ln_pi_rate ln_pi_rate(2007) ln_pi_rate(2006) ln_pi_rate(2005) ///
ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995)  , trunit(`trs') trperiod(`trp') keep(synth_`trs', replace)
clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_`trs'
drop _m
replace _W_Weight = 1 if state==`trs'
reg ln_pi_rate max_offer  if flu_year>=`trp'  [aweight=_W_Weight]
gen treatment = (state==`trs')
collapse pi_rate [aw=_W], by(flu_year treatment)
twoway line pi_rate flu_year if treatment==1, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) ylab(3 6 9, angle(0) nogrid) /// 
xtitle("") ytitle("")  title("DC (2008/09)", color(black)) name(DC, replace) ///
||  line pi_rate flu_year if treatment==0, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(`trp',  lcolor(gs10) lwidth(medthin))
clear
use synth_`trs'
gen post_treat = 0 if _time<`trp'
replace post_treat = 1 if _time>=`trp' &  _time!=.
bysort post_treat: gen num = _N if _time!=.
gen diffsq = (_Y_treated- _Y_synthetic)^2
bysort post_treat: egen sumdiffsq = sum(diffsq) if _time!=.
gen rmspe_pre_temp = sqrt((1/num)*sumdiffsq) if post_treat==0
gen rmspe_post_temp = sqrt((1/num)*sumdiffsq) if post_treat==1
egen rmspe_pre = min(rmspe_pre_temp)
egen rmspe_post = min(rmspe_post_temp)
local rmspe_ratio = rmspe_post[1]/rmspe_pre[1]
local rmspe_pre = rmspe_pre[1]
clear
use synth_placebo
local abs_beta = abs(_b[max_offer])
disp `abs_beta'
count if abs(beta)>`abs_beta'
local pval = r(N)/_N
disp `pval'
count if rmspe_ratio>`rmspe_ratio'
local pval_rmspe = r(N)/_N
disp `pval_rmspe'
post synth_estimates ((exp(_b[max_offer])-1)*$mean_pi ) (`pval') (`pval_rmspe') (`rmspe_pre') (_b[max_offer]) (`trs') (`trp') (1) ("DC")
restore 

*GA
local trs=13
local trp=2010
preserve
keep if (state==`trs' | max_offer==0)
synth ln_pi_rate ln_pi_rate(2009) ln_pi_rate(2008) ln_pi_rate(2007) ln_pi_rate(2006) ln_pi_rate(2005) ///
ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995) , trunit(`trs') trperiod(`trp') keep(synth_`trs', replace)
local pval = e(pval_joint_post)
clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_`trs'
drop _m
replace _W_Weight = 1 if state==`trs'
reg ln_pi_rate max_offer  if flu_year>=`trp'  [aweight=_W_Weight]
gen treatment = (state==`trs')
collapse pi_rate [aw=_W], by(flu_year treatment)
twoway line pi_rate flu_year if treatment==1, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) ylab(3 6 9, angle(0) nogrid) /// 
xtitle("") ytitle("")  title("Georgia (2010/11)", color(black)) name(GA, replace) ///
||  line pi_rate flu_year if treatment==0, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(`trp',  lcolor(gs10) lwidth(medthin))
clear
use synth_`trs'
gen post_treat = 0 if _time<`trp'
replace post_treat = 1 if _time>=`trp' &  _time!=.
bysort post_treat: gen num = _N if _time!=.
gen diffsq = (_Y_treated- _Y_synthetic)^2
bysort post_treat: egen sumdiffsq = sum(diffsq) if _time!=.
gen rmspe_pre_temp = sqrt((1/num)*sumdiffsq) if post_treat==0
gen rmspe_post_temp = sqrt((1/num)*sumdiffsq) if post_treat==1
egen rmspe_pre = min(rmspe_pre_temp)
egen rmspe_post = min(rmspe_post_temp)
local rmspe_ratio = rmspe_post[1]/rmspe_pre[1]
local rmspe_pre = rmspe_pre[1]
clear
use synth_placebo
local abs_beta = abs(_b[max_offer])
disp `abs_beta'
count if abs(beta)>`abs_beta'
local pval = r(N)/_N
disp `pval'
count if rmspe_ratio>`rmspe_ratio'
local pval_rmspe = r(N)/_N
disp `pval_rmspe'
post synth_estimates ((exp(_b[max_offer])-1)*$mean_pi ) (`pval') (`pval_rmspe') (`rmspe_pre') (_b[max_offer]) (`trs') (`trp') (1) ("GA")
restore 

*IL
local trs=17
local trp=2010
preserve
keep if (state==`trs' | max_offer==0)
synth ln_pi_rate ln_pi_rate(2009) ln_pi_rate(2008) ln_pi_rate(2007) ln_pi_rate(2006) ln_pi_rate(2005) ///
ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995) , trunit(`trs') trperiod(`trp') keep(synth_`trs', replace)
clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_`trs'
drop _m
replace _W_Weight = 1 if state==`trs'
reg ln_pi_rate max_offer  if flu_year>=`trp'  [aweight=_W_Weight]
gen treatment = (state==`trs')
collapse pi_rate [aw=_W], by(flu_year treatment)
twoway line pi_rate flu_year if treatment==1, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) ylab(3 6 9, angle(0) nogrid) /// 
xtitle("") ytitle("")  title("Illinois (2010/11)", color(black)) name(IL, replace) ///
||  line pi_rate flu_year if treatment==0, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(`trp',  lcolor(gs10) lwidth(medthin))
clear
use synth_`trs'
gen post_treat = 0 if _time<`trp'
replace post_treat = 1 if _time>=`trp' &  _time!=.
bysort post_treat: gen num = _N if _time!=.
gen diffsq = (_Y_treated- _Y_synthetic)^2
bysort post_treat: egen sumdiffsq = sum(diffsq) if _time!=.
gen rmspe_pre_temp = sqrt((1/num)*sumdiffsq) if post_treat==0
gen rmspe_post_temp = sqrt((1/num)*sumdiffsq) if post_treat==1
egen rmspe_pre = min(rmspe_pre_temp)
egen rmspe_post = min(rmspe_post_temp)
local rmspe_ratio = rmspe_post[1]/rmspe_pre[1]
local rmspe_pre = rmspe_pre[1]
clear
use synth_placebo
local abs_beta = abs(_b[max_offer])
disp `abs_beta'
count if abs(beta)>`abs_beta'
local pval = r(N)/_N
disp `pval'
count if rmspe_ratio>`rmspe_ratio'
local pval_rmspe = r(N)/_N
disp `pval_rmspe'
post synth_estimates ((exp(_b[max_offer])-1)*$mean_pi ) (`pval') (`pval_rmspe') (`rmspe_pre') (_b[max_offer]) (`trs') (`trp') (1) ("IL")
restore 

*ME
local trs=23
local trp=2002
preserve
keep if (state==`trs' | max_offer==0)
synth ln_pi_rate ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995) , trunit(`trs') trperiod(`trp') keep(synth_`trs', replace)
clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_`trs'
drop _m
replace _W_Weight = 1 if state==`trs'
reg ln_pi_rate max_offer  if flu_year>=`trp'  [aweight=_W_Weight]
gen treatment = (state==`trs')
collapse pi_rate [aw=_W], by(flu_year treatment)
twoway line pi_rate flu_year if treatment==1, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) ylab(3 6 9, angle(0) nogrid) /// 
xtitle("") ytitle("")  title("Maine (2002/03)", color(black)) name(ME, replace) ///
||  line pi_rate flu_year if treatment==0, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(`trp',  lcolor(gs10) lwidth(medthin))
clear
use synth_`trs'
gen post_treat = 0 if _time<`trp'
replace post_treat = 1 if _time>=`trp' &  _time!=.
bysort post_treat: gen num = _N if _time!=.
gen diffsq = (_Y_treated- _Y_synthetic)^2
bysort post_treat: egen sumdiffsq = sum(diffsq) if _time!=.
gen rmspe_pre_temp = sqrt((1/num)*sumdiffsq) if post_treat==0
gen rmspe_post_temp = sqrt((1/num)*sumdiffsq) if post_treat==1
egen rmspe_pre = min(rmspe_pre_temp)
egen rmspe_post = min(rmspe_post_temp)
local rmspe_ratio = rmspe_post[1]/rmspe_pre[1]
local rmspe_pre = rmspe_pre[1]
clear
use synth_placebo
local abs_beta = abs(_b[max_offer])
disp `abs_beta'
count if abs(beta)>`abs_beta'
local pval = r(N)/_N
disp `pval'
count if rmspe_ratio>`rmspe_ratio'
local pval_rmspe = r(N)/_N
disp `pval_rmspe'
post synth_estimates ((exp(_b[max_offer])-1)*$mean_pi ) (`pval') (`pval_rmspe') (`rmspe_pre') (_b[max_offer]) (`trs') (`trp') (1) ("ME")
restore 

*MD
local trs=24
local trp=2008
preserve
keep if (state==`trs' | max_offer==0)
synth ln_pi_rate ln_pi_rate(2007) ln_pi_rate(2006) ln_pi_rate(2005) ///
ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995) , trunit(`trs') trperiod(`trp') keep(synth_`trs', replace)
clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_`trs'
drop _m
replace _W_Weight = 1 if state==`trs'
reg ln_pi_rate max_offer  if flu_year>=`trp'  [aweight=_W_Weight]
gen treatment = (state==`trs')
collapse pi_rate [aw=_W], by(flu_year treatment)
twoway line pi_rate flu_year if treatment==1, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) ylab(3 6 9, angle(0) nogrid) /// 
xtitle("") ytitle("")  title("Maryland (2008/09)", color(black)) name(MD, replace) ///
||  line pi_rate flu_year if treatment==0, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(`trp',  lcolor(gs10) lwidth(medthin))
clear
use synth_`trs'
gen post_treat = 0 if _time<`trp'
replace post_treat = 1 if _time>=`trp' &  _time!=.
bysort post_treat: gen num = _N if _time!=.
gen diffsq = (_Y_treated- _Y_synthetic)^2
bysort post_treat: egen sumdiffsq = sum(diffsq) if _time!=.
gen rmspe_pre_temp = sqrt((1/num)*sumdiffsq) if post_treat==0
gen rmspe_post_temp = sqrt((1/num)*sumdiffsq) if post_treat==1
egen rmspe_pre = min(rmspe_pre_temp)
egen rmspe_post = min(rmspe_post_temp)
local rmspe_ratio = rmspe_post[1]/rmspe_pre[1]
local rmspe_pre = rmspe_pre[1]
clear
use synth_placebo
local abs_beta = abs(_b[max_offer])
disp `abs_beta'
count if abs(beta)>`abs_beta'
local pval = r(N)/_N
disp `pval'
count if rmspe_ratio>`rmspe_ratio'
local pval_rmspe = r(N)/_N
disp `pval_rmspe'
post synth_estimates ((exp(_b[max_offer])-1)*$mean_pi ) (`pval') (`pval_rmspe') (`rmspe_pre') (_b[max_offer]) (`trs') (`trp') (1) ("MD")
restore 

*MA
local trs=25
local trp=2009
preserve
keep if (state==`trs' | max_offer==0)
synth ln_pi_rate ln_pi_rate(2008) ln_pi_rate(2007) ln_pi_rate(2006) ln_pi_rate(2005) ///
ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995) , trunit(`trs') trperiod(`trp') keep(synth_`trs', replace)
clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_`trs'
drop _m
replace _W_Weight = 1 if state==`trs'
reg ln_pi_rate max_offer  if flu_year>=`trp'  [aweight=_W_Weight]
gen treatment = (state==`trs')
collapse pi_rate [aw=_W], by(flu_year treatment)
twoway line pi_rate flu_year if treatment==1, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) ylab(3 6 9, angle(0) nogrid) /// 
xtitle("") ytitle("")  title("Massachusetts (2009/10)", color(black)) name(MA, replace) ///
||  line pi_rate flu_year if treatment==0, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(`trp',  lcolor(gs10) lwidth(medthin))
clear
use synth_`trs'
gen post_treat = 0 if _time<`trp'
replace post_treat = 1 if _time>=`trp' &  _time!=.
bysort post_treat: gen num = _N if _time!=.
gen diffsq = (_Y_treated- _Y_synthetic)^2
bysort post_treat: egen sumdiffsq = sum(diffsq) if _time!=.
gen rmspe_pre_temp = sqrt((1/num)*sumdiffsq) if post_treat==0
gen rmspe_post_temp = sqrt((1/num)*sumdiffsq) if post_treat==1
egen rmspe_pre = min(rmspe_pre_temp)
egen rmspe_post = min(rmspe_post_temp)
local rmspe_ratio = rmspe_post[1]/rmspe_pre[1]
local rmspe_pre = rmspe_pre[1]
clear
use synth_placebo
local abs_beta = abs(_b[max_offer])
disp `abs_beta'
count if abs(beta)>`abs_beta'
local pval = r(N)/_N
disp `pval'
count if rmspe_ratio>`rmspe_ratio'
local pval_rmspe = r(N)/_N
disp `pval_rmspe'
post synth_estimates ((exp(_b[max_offer])-1)*$mean_pi ) (`pval') (`pval_rmspe') (`rmspe_pre') (_b[max_offer]) (`trs') (`trp') (1) ("MA")
restore 

*NE
local trs=31
local trp=2011
preserve
keep if (state==`trs' | max_offer==0)
synth ln_pi_rate  ln_pi_rate(2010) ln_pi_rate(2009) ln_pi_rate(2008) ln_pi_rate(2007) ln_pi_rate(2006) ln_pi_rate(2005) ///
ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995) , trunit(`trs') trperiod(`trp') keep(synth_`trs', replace)
clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_`trs'
drop _m
replace _W_Weight = 1 if state==`trs'
reg ln_pi_rate max_offer  if flu_year>=`trp'  [aweight=_W_Weight]
gen treatment = (state==`trs')
collapse pi_rate [aw=_W], by(flu_year treatment)
twoway line pi_rate flu_year if treatment==1, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) ylab(3 6 9, angle(0) nogrid) /// 
xtitle("") ytitle("")  title("Nebraska (2011/12)", color(black)) name(NE, replace) ///
||  line pi_rate flu_year if treatment==0, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(`trp',  lcolor(gs10) lwidth(medthin))
clear
use synth_`trs'
gen post_treat = 0 if _time<`trp'
replace post_treat = 1 if _time>=`trp' &  _time!=.
bysort post_treat: gen num = _N if _time!=.
gen diffsq = (_Y_treated- _Y_synthetic)^2
bysort post_treat: egen sumdiffsq = sum(diffsq) if _time!=.
gen rmspe_pre_temp = sqrt((1/num)*sumdiffsq) if post_treat==0
gen rmspe_post_temp = sqrt((1/num)*sumdiffsq) if post_treat==1
egen rmspe_pre = min(rmspe_pre_temp)
egen rmspe_post = min(rmspe_post_temp)
local rmspe_ratio = rmspe_post[1]/rmspe_pre[1]
local rmspe_pre = rmspe_pre[1]
clear
use synth_placebo
local abs_beta = abs(_b[max_offer])
disp `abs_beta'
count if abs(beta)>`abs_beta'
local pval = r(N)/_N
disp `pval'
count if rmspe_ratio>`rmspe_ratio'
local pval_rmspe = r(N)/_N
disp `pval_rmspe'
post synth_estimates ((exp(_b[max_offer])-1)*$mean_pi ) (`pval') (`pval_rmspe') (`rmspe_pre') (_b[max_offer]) (`trs') (`trp') (1) ("NE")
restore 


*NH
local trs=33
local trp=2005
preserve
keep if (state==`trs' | max_offer==0)
synth ln_pi_rate ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995) , trunit(`trs') trperiod(`trp') keep(synth_`trs', replace)
clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_`trs'
drop _m
replace _W_Weight = 1 if state==`trs'
reg ln_pi_rate max_offer  if flu_year>=`trp'  [aweight=_W_Weight]
gen treatment = (state==`trs')
collapse pi_rate [aw=_W], by(flu_year treatment)
twoway line pi_rate flu_year if treatment==1, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) ylab(3 6 9, angle(0) nogrid) /// 
xtitle("") ytitle("")  title("New Hampshire (2005/06)", color(black)) name(NH, replace) ///
||  line pi_rate flu_year if treatment==0, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(`trp',  lcolor(gs10) lwidth(medthin))
clear
use synth_`trs'
gen post_treat = 0 if _time<`trp'
replace post_treat = 1 if _time>=`trp' &  _time!=.
bysort post_treat: gen num = _N if _time!=.
gen diffsq = (_Y_treated- _Y_synthetic)^2
bysort post_treat: egen sumdiffsq = sum(diffsq) if _time!=.
gen rmspe_pre_temp = sqrt((1/num)*sumdiffsq) if post_treat==0
gen rmspe_post_temp = sqrt((1/num)*sumdiffsq) if post_treat==1
egen rmspe_pre = min(rmspe_pre_temp)
egen rmspe_post = min(rmspe_post_temp)
local rmspe_ratio = rmspe_post[1]/rmspe_pre[1]
local rmspe_pre = rmspe_pre[1]
clear
use synth_placebo
local abs_beta = abs(_b[max_offer])
disp `abs_beta'
count if abs(beta)>`abs_beta'
local pval = r(N)/_N
disp `pval'
count if rmspe_ratio>`rmspe_ratio'
local pval_rmspe = r(N)/_N
disp `pval_rmspe'
post synth_estimates ((exp(_b[max_offer])-1)*$mean_pi ) (`pval') (`pval_rmspe') (`rmspe_pre') (_b[max_offer]) (`trs') (`trp') (1) ("NH")
restore 


*NY
local trs=36
local trp=2013
preserve
keep if (state==`trs' | max_offer==0)
synth ln_pi_rate ln_pi_rate(2012) ln_pi_rate(2011) ln_pi_rate(2010) ln_pi_rate(2009) ln_pi_rate(2008) ln_pi_rate(2007) ln_pi_rate(2006) ln_pi_rate(2005) ///
ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995) , trunit(`trs') trperiod(`trp') keep(synth_`trs', replace)
clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_`trs'
drop _m
replace _W_Weight = 1 if state==`trs'
reg ln_pi_rate max_offer  if flu_year>=`trp'  [aweight=_W_Weight]
gen treatment = (state==`trs')
collapse pi_rate [aw=_W], by(flu_year treatment)
twoway line pi_rate flu_year if treatment==1, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) ylab(3 6 9, angle(0) nogrid) /// 
 xtitle("") ytitle("") title("New York (2013/14)", color(black)) name(NY, replace) ///
||  line pi_rate flu_year if treatment==0, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(`trp',  lcolor(gs10) lwidth(medthin))
clear
use synth_`trs'
gen post_treat = 0 if _time<`trp'
replace post_treat = 1 if _time>=`trp' &  _time!=.
bysort post_treat: gen num = _N if _time!=.
gen diffsq = (_Y_treated- _Y_synthetic)^2
bysort post_treat: egen sumdiffsq = sum(diffsq) if _time!=.
gen rmspe_pre_temp = sqrt((1/num)*sumdiffsq) if post_treat==0
gen rmspe_post_temp = sqrt((1/num)*sumdiffsq) if post_treat==1
egen rmspe_pre = min(rmspe_pre_temp)
egen rmspe_post = min(rmspe_post_temp)
local rmspe_ratio = rmspe_post[1]/rmspe_pre[1]
local rmspe_pre = rmspe_pre[1]
clear
use synth_placebo
local abs_beta = abs(_b[max_offer])
disp `abs_beta'
count if abs(beta)>`abs_beta'
local pval = r(N)/_N
disp `pval'
count if rmspe_ratio>`rmspe_ratio'
local pval_rmspe = r(N)/_N
disp `pval_rmspe'
post synth_estimates ((exp(_b[max_offer])-1)*$mean_pi ) (`pval') (`pval_rmspe') (`rmspe_pre') (_b[max_offer]) (`trs') (`trp') (1) ("NY")
restore 

*OK
local trs=40
local trp=2009
preserve
keep if (state==`trs' | max_offer==0)
synth ln_pi_rate ln_pi_rate(2008) ln_pi_rate(2007) ln_pi_rate(2006) ln_pi_rate(2005) ///
ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995) , trunit(`trs') trperiod(`trp') keep(synth_`trs', replace)
clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_`trs'
drop _m
replace _W_Weight = 1 if state==`trs'
reg ln_pi_rate max_offer  if flu_year>=`trp'  [aweight=_W_Weight]
gen treatment = (state==`trs')
collapse pi_rate [aw=_W], by(flu_year treatment)
twoway line pi_rate flu_year if treatment==1, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) ylab(3 6 9, angle(0) nogrid) /// 
xtitle("") ytitle("")  title("Oklahoma (2009/10)", color(black)) name(OK, replace) ///
||  line pi_rate flu_year if treatment==0, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(`trp',  lcolor(gs10) lwidth(medthin))
clear
use synth_`trs'
gen post_treat = 0 if _time<`trp'
replace post_treat = 1 if _time>=`trp' &  _time!=.
bysort post_treat: gen num = _N if _time!=.
gen diffsq = (_Y_treated- _Y_synthetic)^2
bysort post_treat: egen sumdiffsq = sum(diffsq) if _time!=.
gen rmspe_pre_temp = sqrt((1/num)*sumdiffsq) if post_treat==0
gen rmspe_post_temp = sqrt((1/num)*sumdiffsq) if post_treat==1
egen rmspe_pre = min(rmspe_pre_temp)
egen rmspe_post = min(rmspe_post_temp)
local rmspe_ratio = rmspe_post[1]/rmspe_pre[1]
local rmspe_pre = rmspe_pre[1]
clear
use synth_placebo
local abs_beta = abs(_b[max_offer])
disp `abs_beta'
count if abs(beta)>`abs_beta'
local pval = r(N)/_N
disp `pval'
count if rmspe_ratio>`rmspe_ratio'
local pval_rmspe = r(N)/_N
disp `pval_rmspe'
post synth_estimates ((exp(_b[max_offer])-1)*$mean_pi ) (`pval') (`pval_rmspe') (`rmspe_pre') (_b[max_offer]) (`trs') (`trp') (1) ("OK")
restore 


*RI
local trs=44
local trp=2012
preserve
keep if (state==`trs' | max_offer==0)
synth ln_pi_rate ln_pi_rate(2011) ln_pi_rate(2010) ln_pi_rate(2009) ln_pi_rate(2008) ln_pi_rate(2007) ln_pi_rate(2006) ln_pi_rate(2005) ///
ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995) , trunit(`trs') trperiod(`trp') keep(synth_`trs', replace)
clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_`trs'
drop _m
replace _W_Weight = 1 if state==`trs'
reg ln_pi_rate max_offer  if flu_year>=`trp'  [aweight=_W_Weight]
gen treatment = (state==`trs')
collapse pi_rate [aw=_W], by(flu_year treatment)
twoway line pi_rate flu_year if treatment==1, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) ylab(3 6 9, angle(0) nogrid) /// 
xtitle("") ytitle("")  title("Rhode Island (2012/13)", color(black)) name(RI, replace) ///
||  line pi_rate flu_year if treatment==0, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(`trp',  lcolor(gs10) lwidth(medthin))
clear
use synth_`trs'
gen post_treat = 0 if _time<`trp'
replace post_treat = 1 if _time>=`trp' &  _time!=.
bysort post_treat: gen num = _N if _time!=.
gen diffsq = (_Y_treated- _Y_synthetic)^2
bysort post_treat: egen sumdiffsq = sum(diffsq) if _time!=.
gen rmspe_pre_temp = sqrt((1/num)*sumdiffsq) if post_treat==0
gen rmspe_post_temp = sqrt((1/num)*sumdiffsq) if post_treat==1
egen rmspe_pre = min(rmspe_pre_temp)
egen rmspe_post = min(rmspe_post_temp)
local rmspe_ratio = rmspe_post[1]/rmspe_pre[1]
local rmspe_pre = rmspe_pre[1]
clear
use synth_placebo
local abs_beta = abs(_b[max_offer])
disp `abs_beta'
count if abs(beta)>`abs_beta'
local pval = r(N)/_N
disp `pval'
count if rmspe_ratio>`rmspe_ratio'
local pval_rmspe = r(N)/_N
disp `pval_rmspe'
post synth_estimates ((exp(_b[max_offer])-1)*$mean_pi ) (`pval') (`pval_rmspe') (`rmspe_pre') (_b[max_offer]) (`trs') (`trp') (1) ("RI")
restore 

*TN
local trs=47
local trp=2007
preserve
keep if (state==`trs' | max_offer==0)
synth ln_pi_rate ln_pi_rate(2006) ln_pi_rate(2005) ///
ln_pi_rate(2004) ln_pi_rate(2003) ln_pi_rate(2002) ln_pi_rate(2001) ln_pi_rate(2000) ln_pi_rate(1999) ///
ln_pi_rate(1998) ln_pi_rate(1997) ln_pi_rate(1996) ln_pi_rate(1995) , trunit(`trs') trperiod(`trp') keep(synth_`trs', replace)
clear
use flu_hcw_national_analysis
gen _Co_Number = state
merge m:1 _Co_Number using synth_`trs'
drop _m
replace _W_Weight = 1 if state==`trs'
reg ln_pi_rate max_offer  if flu_year>=`trp'  [aweight=_W_Weight]
gen treatment = (state==`trs')
collapse pi_rate [aw=_W], by(flu_year treatment)
twoway line pi_rate flu_year if treatment==1, lcolor(black) lwidth(thick) ylab(#3) graphregion(color(white)) ylab(3 6 9, angle(0) nogrid) /// 
xtitle("") ytitle("") title("Tennessee (2007/08)", color(black)) name(TN, replace) ///
||  line pi_rate flu_year if treatment==0, lcolor(gs10) lwidth(medthick) lpattern(shortdash) xline(`trp',  lcolor(gs10) lwidth(medthin))
clear
use synth_`trs'
gen post_treat = 0 if _time<`trp'
replace post_treat = 1 if _time>=`trp' &  _time!=.
bysort post_treat: gen num = _N if _time!=.
gen diffsq = (_Y_treated- _Y_synthetic)^2
bysort post_treat: egen sumdiffsq = sum(diffsq) if _time!=.
gen rmspe_pre_temp = sqrt((1/num)*sumdiffsq) if post_treat==0
gen rmspe_post_temp = sqrt((1/num)*sumdiffsq) if post_treat==1
egen rmspe_pre = min(rmspe_pre_temp)
egen rmspe_post = min(rmspe_post_temp)
local rmspe_ratio = rmspe_post[1]/rmspe_pre[1]
local rmspe_pre = rmspe_pre[1]
clear
use synth_placebo
local abs_beta = abs(_b[max_offer])
disp `abs_beta'
count if abs(beta)>`abs_beta'
local pval = r(N)/_N
disp `pval'
count if rmspe_ratio>`rmspe_ratio'
local pval_rmspe = r(N)/_N
disp `pval_rmspe'
post synth_estimates ((exp(_b[max_offer])-1)*$mean_pi ) (`pval') (`pval_rmspe') (`rmspe_pre') (_b[max_offer]) (`trs') (`trp') (1) ("TN")
restore 

postclose synth_estimates


*Plot all synthetic control plots on same graph for Figure 2, locate plots and edit legend manually
grc1leg2 ME NH CA TN MD DC MA OK GA IL NE CO RI NY, graphregion(color(white)) 


*************** Figure 3/S2

*Open postfile of regression estimates
clear
use synth_estimates

*Merge in mean population data for weighted estimates
merge m:1 state using  "$PATH\Data\Population\1969-2017 SEER\mean_state_pop"
keep if _m==3

*Get average (unweighted and weighted) regression estimates
sum me 
sum me [aw=pop_tot]
sum beta [aw=pop_tot]

*Generate law categories
gen lawcat = 0
replace lawcat = 1 if stateab=="NH"  | stateab=="GA" | stateab=="NE" 
replace lawcat = 2 if stateab=="ME" | stateab=="CA" | stateab=="TN" | stateab=="MD" | stateab=="MA"  | stateab=="OK"  | stateab=="IL"  | stateab=="DC"
replace lawcat = 3 if stateab=="CO" | stateab=="RI" | stateab=="NY" 

*Plot effect sizes, delete y-axis and add text manually
gsort  - me
gen n= _n

*Plot effect sizes with different colors for different law types
gen me_label = me
replace me_label = me-.2 if me<0
gen plabel = me
replace plabel = me-.3 if me<0
replace plabel = me+.15 if me>0
replace pval_rmspe = round(pval_rmspe, 0.01)
tostring pval_rmspe, replace
gen pval2 = "p=0"
egen pval3 = concat(pval2 pval_rmspe)

*Figure 3: treatment effects by state
*Color code by law category, add mean treatment effect manually.
twoway bar me n if lawcat==1 , graphregion(color(white)) legend(row(4)) xtitle("") xlab("") xscale(lstyle(none))  ylab(-2.5(.5).5)  ///
 ytitle("Difference in Monthly P&I Deaths per 100,000") yline(-0.467, lcolor(gs6) lpattern(dash))  yline(0, lcolor(black)) ///
  color(ebblue) fintensity(25) lcolor(ebblue*.25) barw(.65) ylabel(,angle(0)) ///
  legend(cols(1) order( 1 "No Declination" 2 "Declination" 3 "Declination & Mask")) ///
 || bar me n if lawcat==2 , color(ebblue) fintensity(75) lcolor(ebblue*.75) barw(.65) ///
 || bar me n if lawcat==3 , color(ebblue) fintensity(150) lcolor(ebblue*1.5) barw(.65) ///
 || scat me_label n , mlab(stateab) msym(none) mlabpos(12) mlabcolor(black) mlabsize(medsmall) ///
 || scat plabel n , mlab(pval3) msym(none) mlabpos(12) mlabcolor(black) mlabsize(vsmall)

*Figure S2: pre-treatment RMSPE by state
 sort rmspe
 gen n2= _n
 sort rmspe_pre 
 twoway bar rmspe_pre n2, graphregion(color(white)) color(ebblue) fintensity(25) lcolor(ebblue*.25) barw(.65) ylabel(,angle(0)) ///
  ytitle("Pre-Treatment Root Mean Squared Prediction Error") xtitle("") xlab("") ///
 || scat rmspe_pre n2 , mlab(stateab) msym(none) mlabpos(12) mlabcolor(black) mlabsize(medsmall) legend(off)
 
 
 
  
  
  
  
  
  
  
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
*************** Figure S3: Predicted Trends and Annual Deaths Averted
**************************************************************************

*Import analysis file
cd "$PATH\Influenza National\"
clear
use flu_hcw_national_analysis

*Preserve data that includes the full sample including low match years(for graphs), but run regression only on the analysis sample
preserve
drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 
global controls offer_ltc req_child i.state#c.ym
xtset state
reg ln_pi_rate offer_post  i.ym i.state $controls [aweight=meanpop] ,  cluster(state)
*calculate estimates+CIs in percentage terms
global pct_change = exp(_b[offer_post]) - 1
global pct_change_lci = exp(_b[offer_post] - (invttail(e(df_r),.025))*_se[offer_post])-1
global pct_change_uci = exp(_b[offer_post] + (invttail(e(df_r),.025))*_se[offer_post])-1
restore

*Get predicted values in logs; don't make predictions for the years that are dropped from the sample 
predict predicted_pi if (flu_year!=1997 & flu_year!=2003 & flu_year!=2004 & flu_year !=2007 & flu_year!=2014 )

*Get counterfactual predicted values
gen offer_post2=offer_post
replace offer_post=0
predict predicted_pi_counterfactual if (flu_year!=1997 & flu_year!=2003 & flu_year!=2004 & flu_year !=2007 & flu_year!=2014 )
*confirm that we get the same difference in logs as our regressions
ttest predicted_pi=predicted_pi_counterfactual if offer_post2==1

*Convert predicted values in logs to predicted values in levels
replace predicted_pi = exp(predicted_pi)
replace predicted_pi_counterfactual = exp(predicted_pi_counterfactual)

*Get lower/upper bound predicted values in levels
gen predicted_pi_lci = predicted_pi_counterfactual
replace predicted_pi_lci = predicted_pi_counterfactual*($pct_change_lci + 1) if offer_post2==1

gen predicted_pi_uci = predicted_pi_counterfactual
replace predicted_pi_uci = predicted_pi_counterfactual*($pct_change_uci + 1) if offer_post2==1

*Generate a variable for the percent of the population affected
gen pop_offer = pop_tot*offer_post2

*Generate variable for the treatment-control difference for the right panel(counterfactual - treatment so numbers are positive)
gen tc_diff = pop_tot*(predicted_pi_counterfactual-predicted_pi)/100000
gen tc_diff_lci = pop_tot*(predicted_pi_counterfactual-predicted_pi_lci)/100000
gen tc_diff_uci = pop_tot*(predicted_pi_counterfactual-predicted_pi_uci)/100000

*Collapse data to treatment-group level (weighted by state population), note we use means of the predicted values so it is in terms of monthly P&I mortality
collapse (mean) predicted* (rawsum) pi pop_tot pop_offer tc_diff* [aw=meanpop], by(year month max_offer)

gen flu_year=year
replace flu_year=year-1 if month<7

*Collapse data to influenza-year level
collapse (mean) pi pop_tot pop_offer predicted* (rawsum) tc_diff*, by(flu_year max_offer)

*Generate variable for percent of population subject to laws in each year, note that this is defined the same for adopting and non-adopting states
*bysort flu_year max_offer: egen sumpop=sum(pop_tot)
bysort flu_year: egen sumpop=sum(pop_tot)
gen pct_pop = pop_offer/sumpop


*After collapsing, these were set to equal 0. Make sure they're missing.
replace tc_diff = . if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 
replace tc_diff_lci = . if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 
replace tc_diff_uci = . if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 

 *Make plots
 label var predicted_pi "P&I Mortality Rate"
 label var flu_year "Influenza Year"
 label var tc_diff "Annual Difference in Deaths (Nationwide)"
  
 twoway bar pct_pop flu_year if max_offer==1, name(mort_trends, replace) ///
 xlabel(1995 "95/96" 1997 "97/98" 1999 "99/00" 2001 "01/02" 2003 "03/04" 2005 "05/06" 2007 "07/08" 2009 "09/10" 2011 "11/12" 2013 "13/14" 2015 "15/16", angle(20)) ///
 ytitle("Share of U.S. Population", axis(2)) legend(rows(3))  ylabel(0 "0" .2 "0.2" .4 "0.4" .6 "0.6" .8 "0.8" ,angle(0) axis(2)) ///
 yaxis(2) fcolor(ebblue) fintensity(25) lcolor(ebblue*.25) graphregion(color(white)) ///
|| connected predicted_pi flu_year if max_offer==1, msymb(O) msize(small) mcolor(black) lcolor(black) yaxis(1) ytitle("Monthly P&I Mortality Rate per 100,000 Population", axis(1)) ///
|| connected predicted_pi_counterfactual flu_year  if max_offer==1,  msymb(O) msize(small) mcolor(black) lcolor(black) lpattern(shortdash) yaxis(1) ylabel(,angle(0) axis(1)) ///
 legend(cols(1) rows(4) order( 1 "Share of Population Affected by Laws" 2 "Predicted Outcome in Treated States" 3 "Counterfactual Predicted Outcome in Treated States")) ///
 xtitle("Influenza Year")

  twoway scat tc_diff flu_year if max_offer==1, mcolor(black) msymbol(D)  legend(cols(1) rows(4) order( 1 "Annual Deaths Averted Relative to No-Law Counterfactual" 2 "95% Confidence Interval")) ///
 xtitle("Influenza Year") name(mort_trends_2, replace) ytitle("Total Annual Deaths Averted (Nationwide)") ///
 xlabel(1995 "95/96" 1997 "97/98" 1999 "99/00" 2001 "01/02" 2003 "03/04" 2005 "05/06" 2007 "07/08" 2009 "09/10" 2011 "11/12" 2013 "13/14" 2015 "15/16", angle(20)) ///
  ylabel( ,angle(0))  graphregion(color(white))  fcolor(ebblue) fintensity(25) lcolor(ebblue*.25) ///
  || rcap tc_diff_lci tc_diff_uci flu_year if max_offer==1, lcolor(black) lwidth(medthin)

 
 
 graph combine mort_trends mort_trends_2, graphregion(color(white)) row(1) xsize(2) ysize(1)
 
 
 
 
 
 
 

 
 
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
 ytitle("Monthly P&I Mortality Rate per 100,000 Population", size(small)) ylabel(,angle(0))  ylabel(0(2)8) ///
 legend(cols(2) order(1 "Never" 2 "2002/03 - 2005/06" 3 "2006/07 - 2009/10" 4 "2010/11 - 2013/14" )) lcolor(black) graphregion(color(white)) lwidth(thick) ///
|| connected pi_rate flu_year  if yearbin==1, lcolor(gs10) mcolor(blue) msize(medsmall) msymb(Oh) ///
|| connected pi_rate flu_year  if yearbin==2, lcolor(gs10) mcolor(red) msize(small)  msymb(T) ///
|| connected pi_rate flu_year  if yearbin==3, lcolor(gs10) mcolor(midgreen) msize(small)  msymb(S) 





















************************************************************************************************
************************************************************************************************
****************************** OTHER ANALYSES NOT REPORTED IN TABLES/FIGURES
************************************************************************************************
************************************************************************************************



























************************************************************************************
***************************************** Bilinski and Hatfield Parallel trends test
************************************************************************************

cd "$PATH\Influenza National\"
clear 
use flu_hcw_national_analysis

drop if flu_year==1997 | flu_year==2003 | flu_year==2004 | flu_year ==2007 | flu_year==2014 

global controls offer_ltc req_child  i.state#c.ym

tab yearstooffer

sum pi_rate if max_offer==1 & offer_post==0 [aweight=meanpop]
global mean_pi = r(mean)

*Construct post-treatment indicators
forvalues i = 0/4{
gen year`i' = (yearstooffer==`i')
}
replace year4 = 1 if yearstooffer>=4 & !missing(yearstooffer)

*Construct treatment-group trend
gen treat_trend = max_offer*ym

*Get average value of ym in the pre-and post period
sum ym if offer_post==0 & max_offer==1 [aweight=meanpop]
local meanpre=r(mean)
sum ym if offer_post==1 & max_offer==1 [aweight=meanpop]
local meanpost=r(mean)

xtset state

*Following Bilinski and Hatfield, estimate model without and with treatment group trend, for model with dynamic treatment effects
xtreg ln_pi_rate year0-year4 i.ym [aweight=meanpop] , fe cluster(state)
xtreg ln_pi_rate treat_trend year0-year4  i.ym [aweight=meanpop] , fe cluster(state)

*Report monthly trend in mortality 
lincom _b[treat_trend]
disp ((exp(r(estimate))-1)*$mean_pi )
disp ((exp(r(estimate) - (invttail(r(df),.025))*r(se))-1)*$mean_pi )
disp ((exp(r(estimate) + (invttail(r(df),.025))*r(se))-1)*$mean_pi )







