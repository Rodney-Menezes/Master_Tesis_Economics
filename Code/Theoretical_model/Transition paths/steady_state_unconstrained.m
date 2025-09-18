%----------------------------------------------------------------
% Compute unconstrained capital stock
%----------------------------------------------------------------

vEProductivityTerm      = mProdTransition * (exp(vProdGrid) .^ (1 / (1- nnu))); 	% productivity term in unconstrained capital stock
EOmegaTerm				= sum(vOmegaWeights .* exp(vOmegaGrid .* tthetaHat));
EOmegaTerm2				= sum(vOmegaWeights .* exp(vOmegaGrid));
vCapitalUnconstrained   = ((A * tthetaHat * (wage ^ (-nnu / (1 - nnu))) * (vEProductivityTerm) * ...
							EOmegaTerm) / (qSS * ((1 / bbeta) - (1 - ddelta) * EOmegaTerm2))) .^ (1 / (1 - tthetaHat));

%----------------------------------------------------------------
% Compute minimum savings policy
%----------------------------------------------------------------

% Pre-compute some useful objects for the iteration; CONVENTION: rows = this period, columns = next period
mCapitalPrime			= repmat(vCapitalUnconstrained,[1 nProd]);  % capital choices for this period (in the rows)
mCapitalPrimePrime		= repmat(vCapitalUnconstrained',[nProd 1]); % capital choices for next period (in the columns)
mProdPrime              = repmat(vProdGrid',[nProd 1]);             % productivity shocks for next period (in the columns)
mProfitsPrime           = A * (exp(mProdPrime) .^ (1 / (1 - nnu))) .* ((exp(omegaMin) * mCapitalPrime) .^ tthetaHat) .* ...
                            (wage ^ (-nnu / (1 - nnu)));            % revenue net of labor costs next period
mRHS1                   = mProfitsPrime + qSS * (1 - ddelta) * exp(omegaMin) * mCapitalPrime - ppsi_0;  % part of RHS must pay whether or not exit

% Initialize variables for the iteration
vDebtUnconstrained      = zeros(nProd,1);           % minimum savings policy
iteration               = 1;                        % current iteration
tolerance               = 1e-12;                    % acceptable error before terminating iteration
maxDifference           = 100;                      % maximum difference between iterations over the state space
dampening               = 0;                        % weight placed on old iteration in updating
maxIterations           = 10000;                    % maximum number of iterations before terminating
vErr                    = zeros(maxIterations,1);   % store error at each step of the iteration fore debugging

% Do the iteration
while maxDifference > tolerance && iteration <= maxIterations

    % Do the minimization
    mDebtPrimePrime 		  = repmat(vDebtUnconstrained',[nProd 1]); % debt choices for next period (in the columns)
    mRHS 		              = inflationSS * (mRHS1 + min(zeros(nProd,nProd),-qSS * mCapitalPrimePrime + bbeta * mDebtPrimePrime)); % the RHS
    [vDebtUnconstrainedNew,~] = min(mRHS,[],2); % minimize over the rows

    % Update the iteration
    maxDifference             = max(abs(vDebtUnconstrainedNew - vDebtUnconstrained));
    vErr(iteration,1)         = maxDifference;
    vDebtUnconstrained        = (1 - dampening) * vDebtUnconstrainedNew + dampening * vDebtUnconstrained;
    iteration                 = iteration + 1;

end


%----------------------------------------------------------------
% Compute unconstrained cutoff
%----------------------------------------------------------------

vUnconstrainedCutoff		= qSS * vCapitalUnconstrained - bbeta * vDebtUnconstrained;


%----------------------------------------------------------------
% Compute unconstrained value function v*(z)
%----------------------------------------------------------------

% Compute flow payoff in RHS; CONVENTION: rows = this period, columns = next period
mCapitalPrime			= repmat(vCapitalUnconstrained,[1 nShocks]);  % capital choices for this period (in the rows)
mProdPrime              = repmat(mShocksGrid(:,1)',[nProd 1]);        % productivity shocks for next period (in the columns)
mOmegaPrime             = repmat(mShocksGrid(:,2)',[nProd 1]);        % capital quality shocks for next period (in the columns)
mProfitsPrime           = A * (exp(mProdPrime) .^ (1 / (1 - nnu))) .* ((exp(mOmegaPrime) .* mCapitalPrime) .^ tthetaHat) .* ...
                            (wage ^ (-nnu / (1 - nnu)));            % revenue net of labor costs next period
mFlowPayoffPrime   		= -qSS * mCapitalPrime + bbeta * (mProfitsPrime + qSS * (1 - ddelta) * exp(mOmegaPrime) .* mCapitalPrime - ...
							ppsi_0);	% flow payoff given shock

% Compute expectation
vFlowPayoff				= diag(mShocksTransition * mFlowPayoffPrime');

% Compute value function (don't have to iterate because can directly solve matrix equation)
vValueUnconstrained 	= (eye(nProd) - bbeta * (1 - ppiExit) * mProdTransition) \ vFlowPayoff;
