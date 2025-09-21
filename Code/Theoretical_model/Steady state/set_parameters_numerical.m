%%%
% Declare global variables
%%%

global nProd nCash nCapital nDebt nCashDist nCapitalDist nDebtDist vProdGrid mProdTransition vProdErgodic ...
		prodMinMult prodMaxMult cashMinMult cashMaxMult capitalMinMult capitalMaxMult ...
		debtMinMult debtMaxMult capitalMinDistMult cashPower capitalPower debtPower cashPowerDist ...
		prodMin prodMax capitalMinDist T wageMin wageMax vProdGrid mProdTransition vProdErgodic ...
		vOmegaGrid vOmegaWeights mShocksGrid omegaMin omegaMax nShocks nOmega mShocksTransition lambda0 alpha


%%%
% Size of grids
%%%

% Number of grid points for decision rules
nProd					 = 10;								% idiosyncratic shocks
nCash					 = 50;								% cash on hand
nCapital			     = 200;								% capital
nDebt					 = 200;								% debt

% Number of grid points for distribution
nCashDist				 = 150;								% cash on hand
nCapitalDist		     = 100;								% capital
nDebtDist				 = 100;								% debt

% Number of grid points for capital quality shocks
nOmega					 = 10;

%%%
% Parameters controlling the grid spacing
%%%

% Scale of the grids
prodMinMult					= 2.5;                              % lower bound = -prodMin * standard deviations
prodMaxMult 				= 2.5;                              % upper bound = prodMax * standard deviations
omegaMinMult				= 4;								% lower bound = -omegaMinMult * standard deviations
omegaMaxMult				= 0;								% upper bound = omegaMaxMult * standard deviations
cashMinMult					= 1.05;                             % lower bound = cashMin * min(default cutoff); will be computed later
cashMaxMult					= 1.05;                             % upper bound = cashMax * max(unconstrained cash); will be computed later
capitalMinMult				= .001;                           	% lower bound = capitalMin * min(unconstrained capital); will be computed later (has to be small so k^{\prime} = 0 is a choice)
capitalMaxMult				= 1.2;                            	% upper bound = capitalMax * max(unconstrained capital); will be computed later
debtMinMult					= 1.05;                             % lower bound = min(b0,debtMin * min(minimum savings policy)); will be computed later
debtMaxMult					= 1.05;                             % upper bound = debtMax * max(unconstrained capital); will be computed later
capitalMinDistMult			= 1;								% lower bound = capitalMinDist * k0; upper bound is the same as for choice grid

% Nonlinearity of the grids; 1 = linear; smaller number means more nonlinear
cashPower					= .3;								% cash on hand
cashPowerDist				= .55;								% (may want to be more linear, e.g. 0.75)
capitalPower				= .6;								% capital
debtPower					= .6;								% debt




%%%
% Compute bounds of grids which are already known
%%%

% Productivity shocks
prodMin         		= -prodMinMult * ssigmaProd / sqrt(1 - rrhoProd ^ 2);
prodMax         		= prodMaxMult * ssigmaProd / sqrt(1 - rrhoProd ^ 2);

% Capital quality shocks
omegaMin				= -omegaMinMult * ssigmaOmega;
omegaMax				= omegaMaxMult * ssigmaOmega;

% Capital min for distribution
capitalMinDist			= capitalMinDistMult * k0;






%%%
% Compute transition matrix idiosyncratic shocks
%%%

% TFP shocks
[vProdGrid,mProdTransition,vProdErgodic] = ...
	compute_tauchen_transition(rrhoProd,ssigmaProd,nProd);

% Capital quality shocks
[vOmegaGrid,vOmegaWeights]	= compute_quadrature_rule(ssigmaOmega,nOmega);

% Make a tensor product grid
nShocks								= nProd * nOmega;
[mProdGridShocks,mOmegaGridShocks]	= ndgrid(vProdGrid,vOmegaGrid);
mShocksGrid							= [mProdGridShocks(:) mOmegaGridShocks(:)];

% Make a combined transition matrix
mShocksTransition		= reshape(repmat(reshape(mProdTransition,nProd,nProd,1),[1 1 nOmega]) .* repmat(reshape(vOmegaWeights,...
							1,1,nOmega),[nProd nProd 1]),nProd,nShocks);





%%%
% Parameters for simulating firms in steady state
%%%

% Lifecycle dynamics
T_lifecycle				= 100;

% Panel simulation
tAnnual					= 20;				% number of years
N						= 5*5000;			% number of firms
T_panel					= tAnnual * 4;
tBurn					= tAnnual - 12;		% drop first 12 years, as in Cooper-Haltiwanger


%%%
% Bounds for bisection method for steady state wage
%%%

wageMin         = .8 * wRepSS;
wageMax         = 1.05 * wRepSS;


% ========== NUEVO: parámetros de la restricción financiera =========
%%%
lambda0         = 5;    % coeficiente base de apalancamiento (λ₀)
alpha           = 2;    % sensibilidad del apalancamiento al riesgo soberano (α)
%%% ====================================================================