%----------------------------------------------------------------
% Compute minimum savings policy, unconstrained cutoff, and unconstrained value function
%----------------------------------------------------------------

% Useful matrices for later on
mProdPrime               = repmat(vProdGrid',[nProd 1]); % productivity shocks for next period (in the columns)
mProdPrime2				 = repmat(mShocksGrid(:,1),[1 nProd]);
mOmegaPrime2			 = repmat(mShocksGrid(:,2),[1 nProd]);


for t = T:-1:1

	%%%
	% Unconstrained capital stock
	%%%

	mCapitalUnconstrainedSeries(:,t)	= ((vA(t+1,1) .* tthetaHat .* (vWage(t+1,1) .^ (-nnu / (1 - nnu))) .* (vEProductivityTerm) * ...
											EOmegaTerm) ./ ((vQ(t,1) ./ vSDF(t,1)) - vQ(t+1,1) * (1 - ddelta) * EOmegaTerm2)) .^ (1 / (1 - tthetaHat));

	%%%
	% Minimum savings policy
	%%%

	% Pre-compute some useful objects
	mCapitalPrime	                  = repmat(mCapitalUnconstrainedSeries(:,t),[1 nProd]); % capital choices for this period (in the rows)
	mCapitalPrimePrime                = repmat(mCapitalUnconstrainedSeries(:,t+1)',[nProd 1]); % capital choices for next period (in the columns)
	mProfitsPrime                     = vA(t+1,1) * (exp(mProdPrime) .^ (1 / (1 - nnu))) .* ((exp(omegaMin) * mCapitalPrime) .^ tthetaHat) .* ...
	                    				(vWage(t+1,1) ^ (-nnu / (1 - nnu))); % revenue net of labor costs next period
	mRHS1                             = mProfitsPrime + vQ(t+1,1) * (1 - ddelta) * exp(omegaMin) * mCapitalPrime - ppsi_0; % part of RHS must pay whether or not exit

	% Do the minimization
	mDebtPrimePrime                   = repmat(mDebtUnconstrainedSeries(:,t+1)',[nProd 1]); % debt choices for next period (in the columns)
	mRHS                              = vAggregateInflation(t+1,1) * (mRHS1 + min(zeros(nProd,nProd),...
	                              		-vQ(t+1,1) * mCapitalPrimePrime + vSDF(t+1,1) * mDebtPrimePrime)); % the RHS
	[mDebtUnconstrainedSeries(:,t),~] = min(mRHS,[],2); % minimize over the rows

	%%%
	% Unconstrained cutoff
	%%%

	mUnconstrainedCutoffSeries(:,t) = vQ(t,1) * mCapitalUnconstrainedSeries(:,t) - vSDF(t,1) * mDebtUnconstrainedSeries(:,t);

	%%%
	% Unconstrained value function v*(z)
	%%%

	% Compute flow payoff in RHS
	mCapitalPrime			             = repmat(mCapitalUnconstrainedSeries(:,t)',[nShocks 1]); % capital choices for this period (in the columns)
	mProfitsPrime                        = vA(t+1,1) * (exp(mProdPrime2) .^ (1 / (1 - nnu))) .* ((exp(mOmegaPrime2) .* mCapitalPrime) .^ tthetaHat) .* ...
	                    					(vWage(t+1,1) ^ (-nnu / (1 - nnu))); % revenue net of labor costs next period
	mFlowPayoffPrime                     = -vQ(t,1) * mCapitalPrime + vSDF(t,1) * (mProfitsPrime + vQ(t+1,1) * (1 - ddelta) * exp(mOmegaPrime2) .* mCapitalPrime - ...
											ppsi_0);	% flow payoff given shock

	% Take expectation
	vFlowPayoff 						 = diag(mShocksTransition * mFlowPayoffPrime); 	% take expectation
	mValueUnconstrainedSeries(:,t)       = vFlowPayoff + vSDF(t,1) * (1 - ppiExit) * ...
	                                  		mProdTransition * mValueUnconstrainedSeries(:,t+1);

end
