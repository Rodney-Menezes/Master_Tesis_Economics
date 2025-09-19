function [c,ceq] = steady_state_decisionConstraintFunction(x,iProd,iCash,vDefaultCutoff,mValue,wage,sParms)

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


% NUEVO: parámetros de restricción financiera
    lambda0 = 5;    % λ₀: coeficiente base de apalancamiento
    alpha   = 2;      % α: sensibilidad al riesgo soberano
    s_t     = 0;        % riesgo soberano actual
 % fin NUEVO

% Extract candidate choices
kPrime	= max(x(1),capitalMin + 1e-10);
bPrime	= x(2);

% NUEVO: restricción b' ≤ λ(s_t) · n_jt
    n_j      = mCashGrid(iProd,iCash);            % patrimonio neto actual
    lambda_s = lambda0 * exp(-alpha * s_t);       % λ(s_t)
    c        = bPrime - lambda_s * n_j;           % obliga c ≤ 0
% fin NUEVO

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
% Compute dividends and return constraint
%---------------------------------------------------------------

ceq		= qSS * kPrime - mCashGrid(iProd,iCash) - debtPrice * bPrime;
%c		= [];
