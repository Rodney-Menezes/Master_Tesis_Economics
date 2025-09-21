%---------------------------------------------------------------
% Compute evolution of net worth for Khan-Senga-Thomas' "Type 1" firms
% (this will only be used as an initial guess for the numerical optimizer)
%---------------------------------------------------------------

% Compute capital and debt decisions
vCapitalType1       = reshape(repmat(mCapitalUnconstrainedSeries(:,t),[1 nCash]),nState,1);
vDebtType1          = (vQ(t,1) * vCapitalType1 - mStateGrid(:,2)) / vSDF(t,1);

% Compute evolution of cash on hand
mProdPrimeType1         = repmat(mShocksGrid(:,1),[1 nState]);
mOmegaPrimeType1        = repmat(mShocksGrid(:,2),[1 nState]);
mCapitalPrimePrimeType1 = repmat(vCapitalType1',[nShocks 1]);
mDebtPrimePrimeType1    = repmat(vDebtType1',[nShocks 1]);
mProfitsPrimeType1      = vA(t+1,1) * (exp(mProdPrimeType1) .^ (1 / (1 - nnu))) .* ((exp(mOmegaPrimeType1) .* mCapitalPrimePrimeType1) .^ tthetaHat) * ...
							(vWage(t+1,1) ^ (-nnu / (1 - nnu)));
mCashPrimeType1         = mProfitsPrimeType1 + vQ(t+1,1) * (1 - ddelta) * exp(mOmegaPrimeType1) .* mCapitalPrimePrimeType1 - ...
							(mDebtPrimePrimeType1 / vAggregateInflation(t+1,1)) - ppsi_0;

% Compute indicator for default
mDefaultPrimeType1Exit  = mCashPrimeType1 <= 0;
mDefaultPrimeType1Cont  = (mCashPrimeType1 <= reshape(repmat(mDefaultCutoffSeries(:,t+1),[1 nOmega * nState]),nShocks,nState));

% Compute indicator for Type 1
mType1Indicator         = 1 - reshape(min(ones(nState,1),sum(mDefaultPrimeType1Exit + mDefaultPrimeType1Cont,1)'),...
                            nProd,nCash);

% Reshape decisions
mCapitalType1           = reshape(vCapitalType1,nProd,nCash);
mDebtType1              = reshape(vDebtType1,nProd,nCash);


%---------------------------------------------------------------
% Other preliminaries
%---------------------------------------------------------------

% Pre-allocate some useful matrices
mCapitalPrimeContinuous = zeros(nProd,nCash);
mDebtPrimeContinuous    = zeros(nProd,nCash);
mDividendsContinuous    = zeros(nProd,nCash);
mDebtPriceContinuous    = zeros(nProd,nCash);

% Stay within bounds of choice space (for numerical stability)
lb                      = [capitalMin;debtMin];
ub                      = [capitalMax;debtMax];

%---------------------------------------------------------------
% Do the optimization
%---------------------------------------------------------------

options		= optimoptions('fmincon','Display','off','maxFunctionEvaluations',1000,'MaxIterations',500,...
                'algorithm','interior-point');

for iProd   = 1 : nProd
	parfor iCash   = 1 : nCash

	    if mCashGrid(iProd,iCash) <= mDefaultCutoffSeries(iProd,t)

	        %%%
	        % Check if default
	        %%%

	        mCapitalPrimeContinuous(iProd,iCash)    = capitalMin;
	        mDebtPrimeContinuous(iProd,iCash)       = debtMin;
	        mDividendsContinuous(iProd,iCash)       = 0;
	        mDebtPriceContinuous(iProd,iCash)       = 0;

	    elseif mCashGrid(iProd,iCash) >= mUnconstrainedCutoffSeries(iProd,t)

	        %%%
	        % Check if unconstrained
	        %%%

	        mCapitalPrimeContinuous(iProd,iCash)	= mCapitalUnconstrainedSeries(iProd,t);
	        mDebtPrimeContinuous(iProd,iCash)       = mDebtUnconstrainedSeries(iProd,t);
	        mDividendsContinuous(iProd,iCash)       = mCashGrid(iProd,iCash) - vQ(t,1) * mCapitalUnconstrainedSeries(iProd,t) + ...
														vSDF(t,1) * mDebtUnconstrainedSeries(iProd,t);
	        mDebtPriceContinuous(iProd,iCash)       = vSDF(t,1);

	    elseif mType1Indicator(iProd,iCash) && (mCashGrid(iProd,iCash) < mUnconstrainedCutoffSeries(iProd,t))

	        %%%
	        % Check if type 1 (still do continuous optimization, just with type 1 decision rules as initial guess)
	        %%%

	        % Define the functions for the optimization
	        myfun           = @(x) transition_path_decisionObjectiveFunction(x,iProd,iCash,mDefaultCutoffSeries(:,t+1),...
	                        	reshape(mValueSeries(:,t+1),nProd,nCash),vWage(t+1,1),vA(t+1,1),vQ(t+1,1),vQ(t,1),...
	                        	vAggregateInflation(t+1,1),vSDF(t,1),sParms);
	        mycon           = @(x) transition_path_decisionConstraintFunction(x,iProd,iCash,mDefaultCutoffSeries(:,t+1),...
	                        	reshape(mValueSeries(:,t+1),nProd,nCash),vWage(t+1,1),vA(t+1,1),vQ(t+1,1),vQ(t,1),...
	                        	vAggregateInflation(t+1,1),vSDF(t,1),sParms);

	        % Compute optimum
	        [xopt,fval,exitflag] = fmincon(myfun,[mCapitalPrimeSS(iProd,iCash);mDebtPrimeSS(iProd,iCash)],...
						[],[],[],[],lb,ub,mycon,options);

	        % Store output
	        if exitflag == -2           % optimization routine did not converge

				%%%
				% Try with inequality constraint (seems to work better for low levels of cash on hand)
				%%%

				% Define the constraint function
				 mycon           = @(x) transition_path_decisionConstraintFunction_inequality(x,iProd,iCash,mDefaultCutoffSeries(:,t+1),...
										reshape(mValueSeries(:,t+1),nProd,nCash),vWage(t+1,1),vA(t+1,1),vQ(t+1,1),vQ(t,1),...
										vAggregateInflation(t+1,1),vSDF(t,1),sParms);

				% Compute optimum
				[xopt,fval,exitflag] = fmincon(myfun,[mCapitalPrimeSS(iProd,iCash);mDebtPrimeSS(iProd,iCash)],...
										[],[],[],[],lb,ub,mycon,options);

				% If still did not converge, use initial guess
				if exitflag == -2

					mCapitalPrimeContinuous(iProd,iCash)    = mCapitalPrime(iProd,iCash);
					mDebtPrimeContinuous(iProd,iCash)       = mDebtPrime(iProd,iCash);
					mDividendsContinuous(iProd,iCash)       = mDividends(iProd,iCash);
					mDebtPriceContinuous(iProd,iCash)       = mDebtPriceOptimal(iProd,iCash);

				else

					mCapitalPrimeContinuous(iProd,iCash)    = xopt(1);
					mDebtPrimeContinuous(iProd,iCash)       = xopt(2);
					[mDividendsContinuous(iProd,iCash),~]   = mycon(xopt);
					mDividendsContinuous(iProd,iCash)       = -mDividendsContinuous(iProd,iCash);
					mDebtPriceContinuous(iProd,iCash)       = (vQ(t,1) * xopt(1) - mCashGrid(iProd,iCash) - ...
																mDividendsContinuous(iProd,iCash)) / xopt(2);

				end

				aExitFlagSeries(iProd,iCash,t)				= exitflag;

	     	else                        % optimization routine converged

	            mCapitalPrimeContinuous(iProd,iCash)    = xopt(1);
	            mDebtPrimeContinuous(iProd,iCash)       = xopt(2);
	            [~,mDividendsContinuous(iProd,iCash)]   = mycon(xopt);
	            mDividendsContinuous(iProd,iCash)       = -mDividendsContinuous(iProd,iCash);
	            mDebtPriceContinuous(iProd,iCash)       = (vQ(t,1) * xopt(1) - mCashGrid(iProd,iCash) - ...
															mDividendsContinuous(iProd,iCash)) / xopt(2);

				aExitFlagSeries(iProd,iCash,t)			= exitflag;

	        end

	    else 					% not "type 1" unconstrained

	        % Define the functions for the optimization
	        myfun           = @(x) transition_path_decisionObjectiveFunction(x,iProd,iCash,mDefaultCutoffSeries(:,t+1),...
	                        	reshape(mValueSeries(:,t+1),nProd,nCash),vWage(t+1,1),vA(t+1,1),vQ(t+1,1),vQ(t,1),...
	                        	vAggregateInflation(t+1,1),vSDF(t,1),sParms);
	        mycon           = @(x) transition_path_decisionConstraintFunction(x,iProd,iCash,mDefaultCutoffSeries(:,t+1),...
	                        	reshape(mValueSeries(:,t+1),nProd,nCash),vWage(t+1,1),vA(t+1,1),vQ(t+1,1),vQ(t,1),...
	                        	vAggregateInflation(t+1,1),vSDF(t,1),sParms);

	        % Compute optimum
	        [xopt,fval,exitflag] = fmincon(myfun,[mCapitalPrimeSS(iProd,iCash);mDebtPrimeSS(iProd,iCash)],...
																[],[],[],[],lb,ub,mycon,options);

	        % Store output
	        if exitflag == -2           % optimization routine did not converge

				%%%
				% Try with inequality constraint (seems to work better for low levels of cash on hand)
				%%%

				% Define the constraint function
				 mycon           = @(x) transition_path_decisionConstraintFunction_inequality(x,iProd,iCash,mDefaultCutoffSeries(:,t+1),...
									reshape(mValueSeries(:,t+1),nProd,nCash),vWage(t+1,1),vA(t+1,1),vQ(t+1,1),vQ(t,1),...
									vAggregateInflation(t+1,1),vSDF(t,1),sParms);

				% Compute optimum
				[xopt,fval,exitflag] = fmincon(myfun,[mCapitalPrimeSS(iProd,iCash);mDebtPrimeSS(iProd,iCash)],...
										[],[],[],[],lb,ub,mycon,options);

				% If still did not converge, use initial guess
				if exitflag == -2

					mCapitalPrimeContinuous(iProd,iCash)    = mCapitalPrime(iProd,iCash);
					mDebtPrimeContinuous(iProd,iCash)       = mDebtPrime(iProd,iCash);
					mDividendsContinuous(iProd,iCash)       = mDividends(iProd,iCash);
					mDebtPriceContinuous(iProd,iCash)       = mDebtPriceOptimal(iProd,iCash);

				else

					mCapitalPrimeContinuous(iProd,iCash)    = xopt(1);
					mDebtPrimeContinuous(iProd,iCash)       = xopt(2);
					[mDividendsContinuous(iProd,iCash),~]   = mycon(xopt);
					mDividendsContinuous(iProd,iCash)       = -mDividendsContinuous(iProd,iCash);
					mDebtPriceContinuous(iProd,iCash)       = (vQ(t,1) * xopt(1) - mCashGrid(iProd,iCash) - ...
																mDividendsContinuous(iProd,iCash)) / xopt(2);

				end

				aExitFlagSeries(iProd,iCash,t)				= exitflag;

		    else                        % optimization routine converged

	            mCapitalPrimeContinuous(iProd,iCash)    = xopt(1);
	            mDebtPrimeContinuous(iProd,iCash)       = xopt(2);
	            [~,mDividendsContinuous(iProd,iCash)]   = mycon(xopt);
	            mDividendsContinuous(iProd,iCash)       = -mDividendsContinuous(iProd,iCash);
	            mDebtPriceContinuous(iProd,iCash)       = (vQ(t,1) * xopt(1) - mCashGrid(iProd,iCash) - ...
															mDividendsContinuous(iProd,iCash)) / xopt(2);

				aExitFlagSeries(iProd,iCash,t)			= exitflag;

	        end

	    end                     % end the if statements for (iProd,iCash)

    end
end
