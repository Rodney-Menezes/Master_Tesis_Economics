/* --------------------------------------

AUTHORS: Pablo Ottonello Thomas Winberry

PURPOSE: compute empirical moments using model-simulated data from steady state

DATE CREATED: 3/4/2017
DATE UPDATED: 3/16/2020

NOTES:

--------------------------------------- */

************************************************************
**   Housekeeping                                         **
************************************************************

clear all
estimates clear

set more off
pause on
graph set ps logo off


************************************************************
**   Import and clean data                                **
************************************************************

* Import data
insheet using mQuarterlyPanel.csv, comma
rename v1 firm_id
rename v2 quarter_id
rename v3 in_sample
rename v4 balanced_panel
rename v5 investment
rename v6 cash_flow
rename v7 market_value
rename v8 capital
rename v9 cash
rename v10 debt
rename v11 employment
xtset firm_id quarter_id

* Compute variables of interest
gen i_k 	= investment / capital
gen tobin_q 	= market_value / capital
gen cf_k	= L.cash_flow / capital
gen lev 	= debt / capital
replace lev 	= debt / (capital - debt) if debt < 0

gen pos_debt	= (debt > 0)
gen gross_lev 	= debt * pos_debt / capital

* Save data
save simulated_data.dta, replace


************************************************************
**  Misc. panel stats on Compustat sample           	  **
************************************************************

gen in_compustat = 1 if quarter_id >= 28

corr lev L.lev if in_sample == 1
corr lev L.lev if in_sample == 1 & in_compustat == 1

corr i_k lev if in_sample == 1
corr i_k lev if in_sample == 1 & in_compustat == 1


************************************************************
**   Compare compustat vs. non-compustat sample (table X) **
************************************************************

gen emp_dhs = (employment - L.employment) / ((employment + L.employment) / 2)
summ emp_dhs if in_sample == 1 & in_compustat == .
summ emp_dhs if in_sample == 1 & in_compustat == 1
