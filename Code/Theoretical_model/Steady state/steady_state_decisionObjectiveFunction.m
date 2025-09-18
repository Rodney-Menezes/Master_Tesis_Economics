function objective = steady_state_decisionObjectiveFunction(x,iProd,iCash,vDefaultCutoff,mValue,wage,sParms)

A          = sParms.A;
vProdGrid  = sParms.vProdGrid;
nnu        = sParms.nnu;
tthetaHat  = sParms.tthetaHat;
qSS        = sParms.qSS;
ddelta     = sParms.ddelta;
inflationSS= sParms.inflationSS;
ppsi_0     = sParms.ppsi_0;
aalpha     = sParms.aalpha;
ppiExit    = sParms.ppiExit;
nProd      = sParms.nProd;
bbeta      = sParms.bbeta;
mProdGrid  = sParms.mProdGrid;
mCashGrid  = sParms.mCashGrid;
mProdTransition = sParms.mProdTransition;
cashMin    = sParms.cashMin;
cashMax    = sParms.cashMax;
nShocks    = sParms.nShocks;
mShocksGrid= sParms.mShocksGrid;
mShocksTransition = sParms.mShocksTransition;
nOmega     = sParms.nOmega;
capitalMin = sParms.capitalMin;

% Extract candidate choices
kPrime	= max(x(1),capitalMin + 1e-10);
bPrime	= x(2);

%---------------------------------------------------------------
% Compute next period's cash for different realizations of shock
%---------------------------------------------------------------

vProfitsPrime           = A * (exp(mShocksGrid(:,1)) .^ (1 / (1 - nnu))) .* ((exp(mShocksGrid(:,2)) .* kPrime) .^ tthetaHat) * ...
                            (wage ^ (-nnu / (1 - nnu)));
vCashPrime              = vProfitsPrime + qSS * (1 - ddelta) * exp(mShocksGrid(:,2)) * kPrime - (bPrime / inflationSS) - ppsi_0;
vCashPrime              = max(cashMin * ones(nShocks,1), min(cashMax * ones(nShocks,1),vCashPrime));

%---------------------------------------------------------------
% Compute bond price
%---------------------------------------------------------------

% Recovery value if default
vRecoveryValue      = max(0,min(1,(aalpha * exp(mShocksGrid(:,2)) * qSS * kPrime) / (bPrime / inflationSS)));

% Default probability
vDefaultIndicator   = ppiExit * (vCashPrime <= zeros(nShocks,1)) + (1 - ppiExit) * ...
                        (vCashPrime <= reshape(repmat(vDefaultCutoff,[1 nOmega]),nShocks,1));

% Debt price
debtPrice           = min(bbeta,max(0,mShocksTransition(iProd,:) * (bbeta * (1 - vDefaultIndicator .* (1 - vRecoveryValue)))));

%---------------------------------------------------------------
% Compute expected value function next period
%---------------------------------------------------------------

% Interpolate value function over possible states next period
vValuePrime         = reshape(interpn(mProdGrid,mCashGrid,mValue,mShocksGrid(:,1),vCashPrime),nShocks,1);

% Integrate over exit shock
vValuePrime         = ppiExit * (1 - (vCashPrime <= zeros(nShocks,1))) .* vCashPrime + ...
                        (1 - ppiExit) * (1 - (vCashPrime <= reshape(repmat(vDefaultCutoff,[1 nOmega]),nShocks,1))) .* vValuePrime;

% Integrate over idiosyncratic productivity shock
expectedValue       = mShocksTransition(iProd,:) * vValuePrime;

%---------------------------------------------------------------
% Compute objective function
%---------------------------------------------------------------

objective       = -(mCashGrid(iProd,iCash) - qSS * kPrime + debtPrice * bPrime + bbeta * expectedValue);
