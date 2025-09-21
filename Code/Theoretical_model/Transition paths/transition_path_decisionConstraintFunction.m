function [c,ceq] = transition_path_decisionConstraintFunction(x,iProd,iCash,vDefaultCutoff_prime,...
                        mValue_prime,wage_prime,A_prime,q_prime,q,inflation_prime,sdf,sParms)

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

vProfitsPrime           = A_prime * (exp(mShocksGrid(:,1)) .^ (1 / (1 - nnu))) .* ((exp(mShocksGrid(:,2)) .* kPrime) .^ tthetaHat) * ...
                            (wage_prime ^ (-nnu / (1 - nnu)));
vCashPrime              = vProfitsPrime + q_prime * (1 - ddelta) * exp(mShocksGrid(:,2)) * kPrime - (bPrime / inflation_prime) - ppsi_0;
vCashPrime              = max(cashMin * ones(nShocks,1), min(cashMax * ones(nShocks,1),vCashPrime));


%---------------------------------------------------------------
% Compute bond price
%---------------------------------------------------------------

% Recovery value if default
vRecoveryValue      = max(0,min(1,(aalpha * exp(mShocksGrid(:,2)) * q_prime * kPrime) / (bPrime / inflation_prime)));

% Default probability
vDefaultIndicator   = ppiExit * (vCashPrime <= zeros(nShocks,1)) + (1 - ppiExit) * ...
                        (vCashPrime <= reshape(repmat(vDefaultCutoff_prime,[1 nOmega]),nShocks,1));

% Debt price
debtPrice           = (sdf / inflation_prime) * min(1,max(0,mShocksTransition(iProd,:) * ((1 - vDefaultIndicator .* (1 - vRecoveryValue)))));


%---------------------------------------------------------------
% Compute dividends and return constraint
%---------------------------------------------------------------

c		= q * kPrime - mCashGrid(iProd,iCash) - debtPrice * bPrime;
ceq		= [];
