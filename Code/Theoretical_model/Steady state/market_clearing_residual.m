function residual = market_clearing_residual(wage);

%----------------------------------------------------------------
% Declare global variables
%----------------------------------------------------------------

global ttheta nnu ddelta rrhoProd ssigmaProd aalpha ppsi_0 ppiExit kInitial b0 ...
    ggamma qSS bbeta nSS pSS kRepSS yRepSS wRepSS inflationSS k0 A tthetaHat prodMin prodMax ...
    capitalMinMult capitalMaxMult cashMinMult cashMaxMult debtMinMult debtMaxMult ...
	cashMin cashMax capitalMin capitalMax debtMin debtMax capitalMinDist cashPower capitalPower debtPower ...
	nProd nCash nCapital nDebt nCashDist nCapitalDist nDebtDist vProdGrid mProdTransition vProdErgodic ...
	cashPower capitalPower debtPower nProdEnt mmuEnt ssigmaEnt cashPowerDist ...
	option_continuous_optimization_VFI option_continuous_optimization ...
	nShocks nOmega vOmegaGrid vOmegaWeights mShocksGrid mShocksTransition ssigmaOmega omegaMin omegaMax

%----------------------------------------------------------------
% Perform core tasks
%----------------------------------------------------------------

core_steady_state;

%----------------------------------------------------------------
% Return market clearing residual
%----------------------------------------------------------------

% Compute aggregates in steady state
steady_state_aggregates;

% Compute residual
residual	= aggregateLabor - nSS;
