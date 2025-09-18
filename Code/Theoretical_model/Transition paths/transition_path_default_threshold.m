%----------------------------------------------------------------
% % Pre-compute various objects involved in the debt price schedule; CONVENTION: rows = next period shocks, columns = choices
%----------------------------------------------------------------

% For computing price schedule
mProdPrime               = repmat(mShocksGrid( :,1),[1 nChoices]); % rows are realizations of productivity next period
mOmegaPrime		         = repmat(mShocksGrid( :,2),[1 nChoices]);	% rows are realizations of capital quality next period
mCapitalPrime            = repmat(mChoicesGrid(:,1)',[nShocks 1]); % columns are possible choices of capital
mDebtPrime               = repmat(mChoicesGrid(:,2)',[nShocks 1]); % columns are possible choices of debt
mDefaultCutoffExit       = zeros(nShocks,nChoices); % default cutoff for firms who receive exit shock = 0

% For computing choices
mCapitalPrimeChoices     = repmat(mChoicesGrid(:,1)',[nProd 1]);
mDebtPrimeChoices		 = repmat(mChoicesGrid(:,2)',[nProd 1]);


%----------------------------------------------------------------
% Compute the path of thresholds
%----------------------------------------------------------------

for t = T:-1:1

    %%%
    % Compute recovery value of capital
    %%%

    mRecoveryValue      = max(zeros(nShocks,nChoices),min(ones(nShocks,nChoices),((aalpha * vQ(t+1,1) * exp(mOmegaPrime) .* ...
                            mCapitalPrime) ./ (mDebtPrime / vAggregateInflation(t+1,1)))));   % lender's payoff if default


    %%%
  	% Compute cash prime
  	%%%%

    mProfitPrime        = vA(t+1,1) * (exp(mProdPrime) .^ (1 / (1 - nnu))) .* ((exp(mOmegaPrime) .* mCapitalPrime) .^ tthetaHat) * ...
                          (vWage(t+1,1) ^ (-nnu / (1 - nnu)));                % revenue net of labor cost next period
    mCashPrime          = mProfitPrime + vQ(t+1,1) * (1 - ddelta) * exp(mOmegaPrime) .* mCapitalPrime - (mDebtPrime / vAggregateInflation(t+1,1)) - ...
						              ppsi_0;    % next period's cash on hand


    %%%
    % Compute new default threshold
    %%%

    % Compute default probability
    mDefaultCutoffCont  = reshape(repmat(mDefaultCutoffSeries(:,t+1),[1 nOmega * nChoices]),nShocks,nChoices);  % default cutoff for firms who do not get exit shock
    mDefaultIndicator   = ppiExit * (mCashPrime <= mDefaultCutoffExit) + ...
                            (1 - ppiExit) * (mCashPrime <= mDefaultCutoffCont);     % indicator for default next period

    % Compute debt price schedule
    mDebtPrice              = (vSDF(t,1) / vAggregateInflation(t+1,1)) * min(ones(nProd,nChoices), max(zeros(nProd,nChoices),mShocksTransition * ...
    							           ((1 - mDefaultIndicator .* (ones(nShocks,nChoices) - mRecoveryValue)))));
  	mDebtPriceSeries(:,t)   = mDebtPrice(:);

    % Compute all possible choices of the RHS
    mRHSOptions         = vQ(t,1) * mCapitalPrimeChoices - mDebtPrice .* mDebtPrimeChoices;

    % Maximize
    mDefaultCutoffSeries(:,t)   = min(mRHSOptions,[],2);

end
