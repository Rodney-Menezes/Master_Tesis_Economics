
************************************************************
**   Housekeeping                                         **
************************************************************

clear all
estimates clear

set more off
pause on
graph set ps logo off


************************************************************
**   Run the files                                       **
************************************************************

foreach tPre in 28 36 44 52 60 68 76 84 92 100 108 116 {

	* Import data
	insheet using mTransitionPanel_`tPre'.csv, comma
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
	rename v12 interest_rate
	rename v13 unconstrained

	* Drop useless observations
	keep if in_sample == 1

	* Set as panel
	xtset firm_id quarter_id

	* Compute variables of interest
	gen i_k 	= investment / capital
	gen tobin_q = market_value / capital
	gen cf_k	= cash_flow / capital
	gen lev 	= F.debt / F.capital
	gen lev_gross = lev if lev > 0
	replace lev_gross = 0 if lev_gross == .
	replace lev = (F.debt / (F.capital - F.debt)) if lev < 0
	gen log_capital = log(capital)

	* Create LHS variables
	forvalues h = 1/8 {
		gen dlogk_`h' = F`h'.log_capital - L.log_capital
	}

	* Estimate firm effect for financial position
	egen mean_lev = mean(lev), by(firm_id)
	egen mean_lev_gross = mean(lev), by(firm_id)

	* Within-firm variation in financial position
	gen lev_within = lev - mean_lev
	gen lev_within_gross = lev_gross - mean_lev_gross

	* Drop if outside Compustat sample
	drop if quarter_id <= `tPre' - 4

	* Generate monetary shock
	gen epsilon_m = 0
	forvalues t = 1 / 10 {
		replace epsilon_m = -0.0025 * (.5^(`t' - 1)) if quarter_id == `tPre' + `t'
	}

	* Flip the sign and annualize as in regressions
	replace epsilon_m = -4 * epsilon_m

	* Standardize variables
	egen lev_within_std = std(lev_within)
	egen lev_within_std_gross = std(lev_within_gross)

	* Generate interaction with shock
	gen lev_shock = L.lev_within * epsilon_m
	gen lev_shock_std = L.lev_within_std * epsilon_m
	gen lev_shock_within_gross = L.lev_within_gross * epsilon_m
	gen lev_shock_within_gross_std = L.lev_within_std_gross * epsilon_m


	************************************************************
	**   Run regressions				                      **
	************************************************************

	***
	* Capital with standardized gross leverage
	***

	forvalues h = 1/8 {
		quietly areg dlogk_`h' lev_shock_within_gross_std L.lev_within_std_gross L.cf_k L.log_capital i.quarter_id, a(firm_id)
		quietly est sto m`h'
	}

	esttab m1 m2 m3 m4 m5 m6 m7 m8 ///
		using "../Results/panel_simulations/grosslev_std.csv", append plain ///
		keep(lev_shock_within_gross_std) se r2 compress collabels(, none)


	**
	* Capital without standardized leverage
	***

	forvalues h = 1/8 {
		quietly areg dlogk_`h' lev_shock_within_gross L.lev_within_gross L.cf_k L.log_capital i.quarter_id, a(firm_id)
		quietly est sto m`h'
	}

	esttab m1 m2 m3 m4 m5 m6 m7 m8 ///
		using "../Results/panel_simulations/grosslev.csv", append plain ///
		keep(lev_shock_within_gross) se r2 compress collabels(, none)


	***
	* Clear all
	***

	clear all

}
