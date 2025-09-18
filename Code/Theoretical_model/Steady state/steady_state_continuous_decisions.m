%---------------------------------------------------------------
% Compute evolution of cash hand hand for Khan-Senga-Thomas' "Type 1" firms
% (this will only be used as an initial guess for the numerical optimizer)
%---------------------------------------------------------------

% Compute capital and debt decisions
vCapitalType1           = reshape(repmat(vCapitalUnconstrained,[1 nCash]),nState,1);
vDebtType1              = (qSS * vCapitalType1 - mStateGrid(:,2)) / bbeta;

% Compute evolution of cash on hand
mProdPrimeType1         = repmat(mShocksGrid(:,1),[1 nState]);
mOmegaPrimeType1        = repmat(mShocksGrid(:,2),[1 nState]);
mCapitalPrimePrimeType1 = repmat(vCapitalType1',[nShocks 1]);
mDebtPrimePrimeType1    = repmat(vDebtType1',[nShocks 1]);
mProfitsPrimeType1      = A * (exp(mProdPrimeType1) .^ (1 / (1 - nnu))) .* ((exp(mOmegaPrimeType1) .* mCapitalPrimePrimeType1) .^ tthetaHat) * ...
                            (wage ^ (-nnu / (1 - nnu)));
mCashPrimeType1         = mProfitsPrimeType1 + qSS * (1 - ddelta) * exp(mOmegaPrimeType1) .* mCapitalPrimePrimeType1 - ...
                            (mDebtPrimePrimeType1 / inflationSS) - ppsi_0;

% Compute indicator for default
mDefaultPrimeType1Exit  = (mCashPrimeType1 <= 0);
mDefaultPrimeType1Cont  = (mCashPrimeType1 <= reshape(repmat(vDefaultCutoff,[1 nOmega * nState]),nShocks,nState));

% Compute indicator for Type 1
mType1Indicator         = 1 - reshape(min(ones(nState,1),sum(mDefaultPrimeType1Exit + mDefaultPrimeType1Cont,1)'),...
                            nProd,nCash);

% Reshape decisions
mCapitalType1           = reshape(vCapitalType1,nProd,nCash);
mDebtType1              = reshape(vDebtType1,nProd,nCash);

% Stay within bounds of choice space (for numerical stability)
lb                      = [capitalMin;debtMin];
ub                      = [capitalMax;debtMax];

%---------------------------------------------------------------
% Pre-allocate some useful matrices
%---------------------------------------------------------------

mCapitalPrimeContinuous = zeros(nProd,nCash);
mDebtPrimeContinuous    = zeros(nProd,nCash);
mDividendsContinuous    = zeros(nProd,nCash);
mDebtPriceContinuous    = zeros(nProd,nCash);

mExitFlag1				= zeros(nProd,nCash);
mExitFlag2				= zeros(nProd,nCash);

%---------------------------------------------------------------
% Store global variables in a structure (need to use this rather than globals for parallelization)
%---------------------------------------------------------------

sParms            = struct;
sParms.A          = A;
sParms.vProdGrid  = vProdGrid;
sParms.nnu        = nnu;
sParms.tthetaHat  = tthetaHat;
sParms.qSS        = qSS;
sParms.ddelta     = ddelta;
sParms.inflationSS= inflationSS;
sParms.ppsi_0     = ppsi_0;
sParms.aalpha     = aalpha;
sParms.ppiExit    = ppiExit;
sParms.nProd      = nProd;
sParms.bbeta      = bbeta;
sParms.mProdGrid  = mProdGrid;
sParms.mCashGrid  = mCashGrid;
sParms.mProdTransition = mProdTransition;
sParms.cashMin    = cashMin;
sParms.cashMax    = cashMax;
sParms.nShocks    = nShocks;
sParms.mShocksGrid= mShocksGrid;
sParms.mShocksTransition = mShocksTransition;
sParms.nOmega     = nOmega;
sParms.capitalMin = capitalMin;

%---------------------------------------------------------------
% Do the optimization
%---------------------------------------------------------------

options		= optimoptions('fmincon','Display','off','maxFunctionEvaluations',1000,'MaxIterations',500,...
                'algorithm','sqp');

for iProd   = 1 : nProd
    parfor iCash = 1 : nCash

    if mCashGrid(iProd,iCash) <= vDefaultCutoff(iProd)

        %%%
        % Check if default
        %%%

        mCapitalPrimeContinuous(iProd,iCash)    = capitalMin;
        mDebtPrimeContinuous(iProd,iCash)       = debtMin;
        mDividendsContinuous(iProd,iCash)       = 0;
        mDebtPriceContinuous(iProd,iCash)       = 0;

    elseif mCashGrid(iProd,iCash) >= vUnconstrainedCutoff(iProd)

        %%%
        % Check if unconstrained
        %%%

        mCapitalPrimeContinuous(iProd,iCash)    = vCapitalUnconstrained(iProd);
        mDebtPrimeContinuous(iProd,iCash)       = vDebtUnconstrained(iProd);
        mDividendsContinuous(iProd,iCash)       = mCashGrid(iProd,iCash) - qSS * vCapitalUnconstrained(iProd) + ...
                                                    bbeta * vDebtUnconstrained(iProd);
        mDebtPriceContinuous(iProd,iCash)       = bbeta;

    elseif mType1Indicator(iProd,iCash) && (mCashGrid(iProd,iCash) < vUnconstrainedCutoff(iProd))

        %%%
        % Check if type 1 (still do continuous optimization, just with type 1 decision rules as initial guess)
        %%%

		% Define the functions for the optimization
        myfun           = @(x) steady_state_decisionObjectiveFunction(x,iProd,iCash,vDefaultCutoff,mValue,wage,sParms);
        mycon           = @(x) steady_state_decisionConstraintFunction(x,iProd,iCash,vDefaultCutoff,mValue,wage,sParms);

        % Compute optimum
        [xopt,fval,exitflag] = fmincon(myfun,[mCapitalPrime(iProd,iCash);mDebtPrime(iProd,iCash)],...
		                          [],[],[],[],lb,ub,mycon,options);

	    % For debugging
	    mExitFlag1(iProd,iCash)	= exitflag;

        % Store output
        if exitflag == -2           % optimization routine did not converge

            %%%
            % Try with inequality constraint
            %%%

            % Define the constraint function
            mycon           = @(x) steady_state_decisionConstraintFunction_inequality(x,iProd,iCash,vDefaultCutoff,mValue,wage,sParms);

            % Compute optimum
            [xopt,fval,exitflag] = fmincon(myfun,[mCapitalPrime(iProd,iCash);mDebtPrime(iProd,iCash)],...
            [],[],[],[],lb,ub,mycon,options);

            % For debugging
            mExitFlag2(iProd,iCash)	= exitflag;

            % If still did not converge, use initial guess
            if exitflag == -2

                mCapitalPrimeContinuous(iProd,iCash)    = mCapitalPrime(iProd,iCash);
                mDebtPrimeContinuous(iProd,iCash)       = mDebtPrime(iProd,iCash);
                mDividendsContinuous(iProd,iCash)       = mDividends(iProd,iCash);
                mDebtPriceContinuous(iProd,iCash)       = mDebtPriceOptimal(iProd,iCash);

            else

                mCapitalPrimeContinuous(iProd,iCash)    = xopt(1);
                mDebtPrimeContinuous(iProd,iCash)       = xopt(2);
                [mDividendsContinuous(iProd,iCash),~]   = steady_state_decisionConstraintFunction_inequality(xopt,iProd,iCash,vDefaultCutoff,mValue,wage,sParms);
                mDividendsContinuous(iProd,iCash)       = -mDividendsContinuous(iProd,iCash);
                mDebtPriceContinuous(iProd,iCash)       = (qSS * xopt(1) - mCashGrid(iProd,iCash) - ...
                										mDividendsContinuous(iProd,iCash)) / xopt(2);

            end

        else                        % optimization routine converged

            mCapitalPrimeContinuous(iProd,iCash)    = xopt(1);
            mDebtPrimeContinuous(iProd,iCash)       = xopt(2);
			[~,mDividendsContinuous(iProd,iCash)]   = steady_state_decisionConstraintFunction(xopt,iProd,iCash,vDefaultCutoff,mValue,wage,sParms);
            mDividendsContinuous(iProd,iCash)       = -mDividendsContinuous(iProd,iCash);
            mDebtPriceContinuous(iProd,iCash)       = (qSS * xopt(1) - mCashGrid(iProd,iCash) - ...
                                                        mDividendsContinuous(iProd,iCash)) / xopt(2);

        end

    else                    % not "type 1" or unconstrained

      % Define the functions for the optimization
      myfun           = @(x) steady_state_decisionObjectiveFunction(x,iProd,iCash,vDefaultCutoff,mValue,wage,sParms);
      mycon           = @(x) steady_state_decisionConstraintFunction(x,iProd,iCash,vDefaultCutoff,mValue,wage,sParms);

      % Compute optimum
      [xopt,fval,exitflag] = fmincon(myfun,[mCapitalPrime(iProd,iCash);mDebtPrime(iProd,iCash)],...
                            [],[],[],[],lb,ub,mycon,options);

	  % For debugging
	  mExitFlag1(iProd,iCash)	= exitflag;

      % Store output
      if exitflag == -2           % optimization routine did not converge

  			%%%
  			% Try with inequality constraint
  			%%%

  			% Define the constraint function
  			mycon           = @(x) steady_state_decisionConstraintFunction_inequality(x,iProd,iCash,vDefaultCutoff,mValue,wage,sParms);

  			% Compute optimum
  			[xopt,fval,exitflag] = fmincon(myfun,[mCapitalPrime(iProd,iCash);mDebtPrime(iProd,iCash)],...
  			[],[],[],[],lb,ub,mycon,options);

  			% For debugging
  			mExitFlag2(iProd,iCash)	= exitflag;

  			% If still did not converge, use initial guess
  			if exitflag == -2

  				mCapitalPrimeContinuous(iProd,iCash)    = mCapitalPrime(iProd,iCash);
  				mDebtPrimeContinuous(iProd,iCash)       = mDebtPrime(iProd,iCash);
  				mDividendsContinuous(iProd,iCash)       = mDividends(iProd,iCash);
  				mDebtPriceContinuous(iProd,iCash)       = mDebtPriceOptimal(iProd,iCash);

  			else

  				mCapitalPrimeContinuous(iProd,iCash)    = xopt(1);
  				mDebtPrimeContinuous(iProd,iCash)       = xopt(2);
  				[mDividendsContinuous(iProd,iCash),~]   = steady_state_decisionConstraintFunction_inequality(xopt,iProd,iCash,vDefaultCutoff,mValue,wage,sParms);
  				mDividendsContinuous(iProd,iCash)       = -mDividendsContinuous(iProd,iCash);
  				mDebtPriceContinuous(iProd,iCash)       = (qSS * xopt(1) - mCashGrid(iProd,iCash) - ...
  															mDividendsContinuous(iProd,iCash)) / xopt(2);

  			end

      else                        % optimization routine converged

            mCapitalPrimeContinuous(iProd,iCash)    = xopt(1);
            mDebtPrimeContinuous(iProd,iCash)       = xopt(2);
			[~,mDividendsContinuous(iProd,iCash)]   = steady_state_decisionConstraintFunction(xopt,iProd,iCash,vDefaultCutoff,mValue,wage,sParms);
            mDividendsContinuous(iProd,iCash)       = -mDividendsContinuous(iProd,iCash);
            mDebtPriceContinuous(iProd,iCash)       = (qSS * xopt(1) - mCashGrid(iProd,iCash) - ...
                                                        mDividendsContinuous(iProd,iCash)) / xopt(2);

      end

    end                     % end the if statements for (iProd,iCash)

  end

end							% end the for loop
