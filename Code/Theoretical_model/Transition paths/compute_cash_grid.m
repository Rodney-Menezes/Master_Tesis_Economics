% Declare global variables
global cashMin cashMax nState mProdGrid mCashGrid mStateGrid

% Convenient object for computing the upper bound on cash
vProfitUnconstrained	= A * (exp(vProdGrid) .^ (1 / (1 - nnu))) .* ((exp(omegaMax) * vCapitalUnconstrained) .^ tthetaHat) * ...
							(wage ^ (-nnu / (1 - nnu)));
vCashUnconstrained		= vProfitUnconstrained + qSS * (1 - ddelta) * exp(omegaMax) * vCapitalUnconstrained - (vDebtUnconstrained / inflationSS) - ppsi_0;


% Grid for cash
cashMin                 = cashMinMult * min(vDefaultCutoff); % 1e-4 es un buffer pequeño pero razonable para que no haya problemas.
cashMax                 = cashMaxMult * max(vCashUnconstrained);
vCashRaw				= linspace(0,1,nCash)';
vCashGrid				= cashMin + (cashMax - cashMin) * (vCashRaw .^ (1 / cashPower));

% Compute tensor product grid of choices (useful for computations later on)
nState                  = nProd * nCash;
[mProdGrid,mCashGrid]   = ndgrid(vProdGrid,vCashGrid);
mStateGrid              = [mProdGrid(:) mCashGrid(:)];
