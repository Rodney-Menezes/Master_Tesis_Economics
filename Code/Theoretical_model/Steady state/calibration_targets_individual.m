function vTargets = calibration_targets_individual(vParms)

%----------------------------------------------------------------
% Declare global variables
%----------------------------------------------------------------

global ttheta nnu ddelta rrhoProd ssigmaProd aalpha ppsi_0 ppsi_1 ppiExit kInitial b0 ...
    ggamma qSS bbeta nSS ssigmaM rrhoM pSS kRepSS yRepSS wRepSS inflationSS k0 A tthetaHat ...
	prodMin prodMax cashMin cashMax capitalMin capitalMax debtMin debtMax capitalMinDist cashPower capitalPower debtPower ...
	nProd nCash nCapital nDebt nCashDist nCapitalDist nDebtDist vProdGrid mProdTransition vProdErgodic ...
    prodMinMult prodMaxMult cashMinMult cashMaxMult capitalMinMult capitalMaxMult ...
    debtMinMult debtMaxMult capitalMinDistMult cashPower capitalPower debtPower nProdEnt ...
	mmuEnt ssigmaEnt cashPowerDist option_continuous_optimization_VFI option_continuous_optimization ...
	nShocks nOmega vOmegaGrid vOmegaWeights mShocksGrid mShocksTransition ssigmaOmega omegaMin omegaMax ...
	mProdGrid mCashGrid option_calibration

% Options we need to set
option_continuous_optimization_VFI	         = 0;
option_continuous_optimization 		         = 1;
option_calibration							 = 1;


%----------------------------------------------------------------
% Set parameters
%----------------------------------------------------------------

% Parameters of economic model
set_parameters_model;

% Parameters of numerical approximations
set_parameters_numerical;

%----------------------------------------------------------------
% Compute steady state
%----------------------------------------------------------------

%%%
% Solve for the market clearing wage using bisection
%%%

% Preliminaries
resid                           = 100;
iteration_market_clearing       = 1;
tolerance_market_clearing       = 5e-3;
max_iterations_market_clearing  = 20;

while abs(resid) > tolerance_market_clearing && iteration_market_clearing <= max_iterations_market_clearing

	% Compute residual
	wageSS  = .5 * (wageMin + wageMax);
	resid   = market_clearing_residual(wageSS);
	if isnan(resid)
		resid = -1e4;       % wage so high that no firms enter
	end

	% Update bounds
	if resid > 0        	% labor demand too high ==> increase wage
		wageMin = wageSS;
	else                	% labor demand too low ==> decrease wage
		wageMax = wageSS;
	end

	% Update iteration
	iteration_market_clearing       = iteration_market_clearing + 1;

end

% Set wage
wage		= wageSS;


%%%
% Compute steady state objects
%%%

% Compute objects
core_steady_state;
steady_state_aggregates;

% Calibrate disutility of labor supply
cchi                		= wageSS * (aggregateConsumption ^ (-ssigma));


%%%
% Calibration targets and other analysis
%%%

steady_state_lifecycle;
steady_state_calibration_stats;


%----------------------------------------------------------------
% Return residual
%----------------------------------------------------------------

vTargets    = [100*sdIKAnnual;100*averageGrossLeverage;400*meanDefault;100*fracPosDebt;...
              100*vEmploymentShareAnnual(1,1);100*sum(vEmploymentShareAnnual(2:10,1));...
              100*vExtensiveMarginAnnual(1,1);100*vExtensiveMarginAnnual(2,1)];
