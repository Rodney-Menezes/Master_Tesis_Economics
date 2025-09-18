global vProdDistEnt vCashGridDist nStateDist mProdGridDist mCashGridDist mStateGridDist vDistEnt

%----------------------------------------------------------------
% Distribution of new entrants
%----------------------------------------------------------------

%%%
% Distribution of total factor productivity
%%%

% Upper and lower thresholds
vPlus 					= zeros(1,nProd);
    vPlus(1,nProd) 			= 1e2;
    vPlus(1,1:nProd-1) 		= vProdGrid(2:nProd)';
vMinus					= zeros(1,nProd);
    vMinus(1,1) 			= -1e2;
    vMinus(1,2:nProd) 		= vProdGrid(1:nProd-1)';

% Compute distribution
vProdDistEnt			= normcdf(vPlus,mmuEnt,ssigmaEnt) - normcdf(vMinus,mmuEnt,ssigmaEnt);
vProdDistEnt			= vProdDistEnt' ./ sum(vProdDistEnt);


%%%
% Joint distribution of TFP and capital quality
%%%

% Distribution
vDistEnt			= reshape(repmat(vProdDistEnt,[1 nOmega]),nShocks,1) .* ...
				        reshape(repmat(vOmegaWeights',[nProd 1]),nShocks,1);

% Implied cash on hand for new entrants
vProfitImpliedEnt   = A * (exp(mShocksGrid(:,1)) .^ (1 / (1 - nnu))) .* ((exp(mShocksGrid(:,2)) .* k0) .^ tthetaHat) .* ...
                        (wage ^ (-nnu / (1 - nnu)));
vCashImpliedEnt     = vProfitImpliedEnt + qSS * exp(mShocksGrid(:,2)) * (1 - ddelta) * k0 - b0 - ppsi_0;

% Default decision for new entrants
vContinueExitEnt     = (vCashImpliedEnt >= 0);
vContinueContEnt     = (vCashImpliedEnt >= reshape(repmat(vDefaultCutoff,[1 nOmega]),nShocks,1));


%----------------------------------------------------------------
% Cash on hand/Net worth
%----------------------------------------------------------------

% Cash grid
vCashRaw			= linspace(0,1,nCashDist)';
vCashGridDist		= cashMin+1e-14 + (cashMax-1e-14 - (cashMin+1e-14)) * (vCashRaw .^ (1 / cashPowerDist));


% Compute tensor product grid (useful for computations later on)
nStateDist					    = nProd * nCashDist;
[mProdGridDist,mCashGridDist]   = ndgrid(vProdGrid,vCashGridDist);
mStateGridDist    				= [mProdGridDist(:) mCashGridDist(:)];


%----------------------------------------------------------------
% Compute decision rules over the distribution grids grids
%----------------------------------------------------------------

%%%
% Interpolate policy functions
%%%

% Capital accumulation policy
vCapitalPrimeDist       = interpn(mProdGrid,mCashGrid,mCapitalPrime,...
                            mStateGridDist(:,1),mStateGridDist(:,2));
mCapitalPrimeDist       = reshape(vCapitalPrimeDist,nProd,nCashDist);

% Debt accumulation policy
vDebtPrimeDist          = interpn(mProdGrid,mCashGrid,mDebtPrime,...
                            mStateGridDist(:,1),mStateGridDist(:,2));
mDebtPrimeDist          = reshape(vDebtPrimeDist,nProd,nCashDist);

% Dividend policy
vDividendsDist       	= interpn(mProdGrid,mCashGrid,mDividends,...
                            mStateGridDist(:,1),mStateGridDist(:,2));
mDividendsDist       	= reshape(vDividendsDist,nProd,nCashDist);

%%%
% Evolution of cash on hand
%%%

% Evolution of state variables; CONVENTION: rows = next period's shocks; columns = today's decisions
mCapitalPrimePrimeDist  = repmat(vCapitalPrimeDist',[nShocks 1]);
mDebtPrimePrimeDist		= repmat(vDebtPrimeDist',[nShocks 1]);
mProdPrimeDist			= repmat(mShocksGrid(:,1),[1 nStateDist]);
mOmegaPrimeDist			= repmat(mShocksGrid(:,2),[1 nStateDist]);

% Next period's cash on hand
mProfitPrimeDist		= A * (exp(mProdPrimeDist) .^ (1 / (1 - nnu))) .* ((exp(mOmegaPrimeDist) .* mCapitalPrimePrimeDist) .^ tthetaHat) .* ...
						  (wage ^ (-nnu / (1 - nnu)));
mCashPrimeDist			= mProfitPrimeDist + qSS * (1 - ddelta) * exp(mOmegaPrimeDist) .* mCapitalPrimePrimeDist - (mDebtPrimePrimeDist / inflationSS) - ppsi_0;
vCashPrimeDist			= mCashPrimeDist(:);

% Next period's default decision (to condition on survival)
mContinueExitDist		= (mCashPrimeDist >= 0);
mContinueContDist		= (mCashPrimeDist >= reshape(repmat(vDefaultCutoff,[1 nOmega * nStateDist]),nShocks,nStateDist));
