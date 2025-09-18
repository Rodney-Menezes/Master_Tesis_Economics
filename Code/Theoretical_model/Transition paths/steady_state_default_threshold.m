%----------------------------------------------------------------
% % Pre-compute various objects involved in the debt price schedule; CONVENTION: rows = shocks next period, columns = choices
%----------------------------------------------------------------

% For computing debt price schedule
mProdPrime			 = repmat(mShocksGrid( :,1),[1 nChoices]);	% rows are realizations of productivity next period
mOmegaPrime			 = repmat(mShocksGrid( :,2),[1 nChoices]);	% rows are realizations of capital quality next period
mCapitalPrime        = repmat(mChoicesGrid(:,1)',[nShocks 1]); 	% columns are possible choices of capital
mDebtPrime           = repmat(mChoicesGrid(:,2)',[nShocks 1]); 	% columns are possible choices of debt
mProfitPrime         = A * (exp(mProdPrime) .^ (1 / (1 - nnu))) .* ((exp(mOmegaPrime) .* mCapitalPrime) .^ tthetaHat) * ...
                        (wage ^ (-nnu / (1 - nnu))); % revenue net of labor cost next period
mCashPrime           = mProfitPrime + qSS * (1 - ddelta) .* exp(mOmegaPrime) .* mCapitalPrime - (mDebtPrime / inflationSS) - ...
						ppsi_0; % next period's cash on hand
mDefaultCutoffExit   = zeros(nShocks,nChoices); % default cutoff for firms who receive exit shock = 0
mRecoveryValue       = max(zeros(nShocks,nChoices),min(ones(nShocks,nChoices),aalpha * qSS * ...
						exp(mOmegaPrime) .* mCapitalPrime ./ mDebtPrime)); % lender's payoff if default

% For computing choices
mCapitalPrimeChoices  = repmat(mChoicesGrid(:,1)',[nProd 1]);
mDebtPrimeChoices	  = repmat(mChoicesGrid(:,2)',[nProd 1]);


%----------------------------------------------------------------
% Iterate on the default threshold
%----------------------------------------------------------------

% Initialize various objects used in the iteration
vDefaultCutoff			= vUnconstrainedCutoff;		% initial guess
iteration               = 1;                        % current iteration
tolerance               = 1e-12;                    % acceptable error before terminating iteration
maxDifference           = 100;                      % maximum difference between iterations over the state space
dampening               = 0;                        % weight placed on old iteration in updating
maxIterations           = 2000;                     % maximum number of iterations before terminating
vErr                    = zeros(maxIterations,1);   % store error at each step of the iteration for debugging

% Do the iteration
while maxDifference > tolerance && iteration <= maxIterations

    %%%
    % Minimize the RHS of expression for cutoff
    %%%

    % Compute default probability
    mDefaultCutoffCont  = reshape(repmat(reshape(vDefaultCutoff,nProd,1,1),[1 nOmega nChoices]),nShocks,nChoices);  % default cutoff for firms who do not get exit shock
    mDefaultIndicator   = ppiExit * (mCashPrime <= mDefaultCutoffExit) + ...
                            (1 - ppiExit) * (mCashPrime <= mDefaultCutoffCont);     % indicator for default next period

    % Compute debt price schedule
    mDebtPrice          = min(bbeta * ones(nProd,nChoices), max(zeros(nProd,nChoices),mShocksTransition * ...
					       (bbeta * (1 - mDefaultIndicator .* (ones(nShocks,nChoices) - mRecoveryValue)))));

    % Compute all possible choices of the RHS
    mRHSOptions         = qSS * mCapitalPrimeChoices - mDebtPrice .* mDebtPrimeChoices;

    % Maximize
    vDefaultCutoffNew   = min(mRHSOptions,[],2);

    %%%
    % Update iteration
    %%%

    maxDifference       = max(abs(vDefaultCutoffNew - vDefaultCutoff));
    vErr(iteration,1)   = maxDifference;
    iteration           = iteration + 1;
    vDefaultCutoff      = (1 - dampening) * vDefaultCutoffNew + dampening * vDefaultCutoff;

end

%----------------------------------------------------------------
% Save debt price that comes out of the iteration
%----------------------------------------------------------------

mDebtPrice				= reshape(mDebtPrice,nProd,nCapital,nDebt);
