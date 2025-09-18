%---------------------------------------------------------------
% Compute decisions by backward iteration
%---------------------------------------------------------------

for t = T:-1:1

    %%%
    % Set up choice grids over the state space
    %%%

    % Compute all possible choices of dividends; CONVENTION: rows = state variables, columns = choices
    mDebtPriceChoices   = reshape(repmat(reshape(mDebtPriceSeries(:,t),nProd,1,nChoices),[1 nCash 1]),...
                            nState,nChoices);			% define debt price over entire state space
    mDividendsChoices   = repmat(mStateGrid(:,2),[1 nChoices]) + mDebtPriceChoices .* repmat(mChoicesGrid(:,2)',[nState 1]) - ...
                            vQ(t,1) * repmat(mChoicesGrid(:,1)',[nState 1]);

    % Compute the unconstrained value function and choices over the whole state space (i.e., defined over productivity and cash on hand)
    vValueUnconstrainedState    = reshape(repmat(mValueUnconstrainedSeries(:,t),[1 nCash]),nState,1) + ...
                                    mStateGrid(:,2);    % continuation value of a firm who becomes unconstrained
    vCapitalUnconstrainedState  = reshape(repmat(mCapitalUnconstrainedSeries(:,t),[1 nCash]),nState,1);    % capital choice
    vDebtUnconstrainedState     = reshape(repmat(mDebtUnconstrainedSeries(:,t),[1 nCash]),nState,1);       % debt choice
    vDividendsUnconstrainedState= mStateGrid(:,2) + vSDF(t,1) * vDebtUnconstrainedState - ...
                                        vQ(t,1) * vCapitalUnconstrainedState;                              % dividend payment

    % Compute cutoff over the whole state space
    vUnconstrainedCutoffState   = reshape(repmat(mUnconstrainedCutoffSeries(:,t),[1 nCash]),nState,1);     % cutoff for being unconstrained
    vDefaultCutoffState         = reshape(repmat(mDefaultCutoffSeries(:,t),[1 nCash]),nState,1);           % cutoff for defaulting

    % Compute evolution of cash on hand for computing expectations; CONVENTION: rows = productivity next period; columns = choices
    % mCashPrime computed above for the debt price
  	mProdPrime          = repmat(mShocksGrid(:,1),[1 nChoices]);       % rows are realizations of productivity next period
  	mOmegaPrime         = repmat(mShocksGrid(:,2),[1 nChoices]);       % rows are realizations of capital quality next period
  	mCapitalPrime       = repmat(mChoicesGrid(:,1)',[nShocks 1]);      % columns are possible choices of capital
  	mDebtPrime          = repmat(mChoicesGrid(:,2)',[nShocks 1]);      % columns are possible choices of debt
    mProfitPrime        = vA(t+1,1) * (exp(mProdPrime) .^ (1 / (1 - nnu))) .* ((exp(mOmegaPrime) .* mCapitalPrime) .^ tthetaHat) * ...
							(vWage(t+1,1) ^ (-nnu / (1 - nnu)));                % revenue net of labor cost next period
    mCashPrime          = mProfitPrime + vQ(t+1,1) * (1 - ddelta) * exp(mOmegaPrime) .* mCapitalPrime - (mDebtPrime / vAggregateInflation(t+1,1)) -...
						    ppsi_0;    % next period's cash on hand
    mCashPrime          = min(max(cashMin * ones(nShocks,nChoices),mCashPrime),...
                            cashMax * ones(nShocks,nChoices)); 			% ensures it stays inside bounds of state space; other interpolation gives error
    vCashPrime          = mCashPrime(:);                     			% collapse into large vector for interpolation
    vProdPrime          = reshape(repmat(mShocksGrid(:,1),[1 nChoices]),nShocks * nChoices,1);   % corresponding vector of productivity shocks
    mDefaultExitPrime   = reshape((vCashPrime <= 0),nShocks,nChoices);   % indicator for default, conditional on exit shock
    mDefaultContPrime   = reshape((vCashPrime <= reshape(repmat(mDefaultCutoffSeries(:,t+1),[1 nOmega * nChoices]),...
                            nShocks * nChoices,1)),nShocks,nChoices);      % indicator for default, conditional on no exit shock


    %%%
    % Compute decisions by numerical maximization
    %%%

    % Interpolate next period's value function over next period's possible states
    mValuePrime     = reshape(interpn(mProdGrid,mCashGrid,reshape(mValueSeries(:,t+1),nProd,nCash),...
                        vProdPrime,vCashPrime),nShocks,nChoices); % rows = next period's productivity, columns = choices of k' and b'

    % Compute expected value function next period over exit shock
    mValuePrime     = ppiExit * (1 - mDefaultExitPrime) .* mCashPrime + (1 - ppiExit) * (1 - mDefaultContPrime) .* mValuePrime;

    % Compute expected value function over next period's productivity shock
    mEValuePrime    = mShocksTransition * mValuePrime;
    mEValuePrime    = reshape(repmat(reshape(mEValuePrime,nProd,1,nChoices),...
                        [1 nCash 1]),nState,nChoices);  				% expand so that rows = state variables, columns = choices

    % Compute objective function for the firm
    mObjective                          = mDividendsChoices + vSDF(t,1) * mEValuePrime;
    mObjective(mDividendsChoices < 0)   = -1e4;             			% enforce non-negativity constraint on dividends

    % Do the optimization
    [vValueNew,vChoiceIndices]      = max(mObjective,[],2);

    % Extract policy functions
    vCapitalPrime                   = mChoicesGrid(vChoiceIndices,1);   % capital policy function
    vDebtPrime                      = mChoicesGrid(vChoiceIndices,2);   % debt policy function
    vDividends                      = zeros(nState,1);                  % dividend policy function
    vDebtPriceOptimal               = zeros(nState,1);                  % implied debt price
    for iState  = 1 : nState
        vDividends(iState,1)        = mDividendsChoices(iState,vChoiceIndices(iState,1));
        vDebtPriceOptimal(iState,1) = mDebtPriceChoices(iState,vChoiceIndices(iState,1));
    end


    %%%
    % Correct at the boundaries of the state space
    %%%

    % Unconstrained firms
    vValueNew(mStateGrid(:,2) >= vUnconstrainedCutoffState)		    = ...
		vValueUnconstrainedState(mStateGrid(:,2) >= vUnconstrainedCutoffState);
    vCapitalPrime(mStateGrid(:,2) >= vUnconstrainedCutoffState)	    = ...
		vCapitalUnconstrainedState(mStateGrid(:,2) >= vUnconstrainedCutoffState);
	vDebtPrime(mStateGrid(:,2) >= vUnconstrainedCutoffState)		= ...
		vDebtUnconstrainedState(mStateGrid(:,2) >= vUnconstrainedCutoffState);
	vDividends(mStateGrid(:,2) >= vUnconstrainedCutoffState)		= ...
		vDividendsUnconstrainedState(mStateGrid(:,2) >= vUnconstrainedCutoffState);
    vDebtPriceOptimal(mStateGrid(:,2) >= vUnconstrainedCutoffState)  = vSDF(t,1);


    % Defaulting firms
  	vValueNew(mStateGrid(:,2) <= vDefaultCutoffState)			= 0;
  	vCapitalPrime(mStateGrid(:,2) <= vDefaultCutoffState)		= capitalMin;
  	vDebtPrime(mStateGrid(:,2) <= vDefaultCutoffState)			= debtMin;
  	vDividends(mStateGrid(:,2) <= vDefaultCutoffState)			= 0;


    %%%
    % Refine decisions by continuous optimization, if requested
    %%%

    if option_continuous_optimization == 1

        % Rename and reshape some stuff
        mCapitalPrime			= reshape(vCapitalPrime,nProd,nCash);
        mDebtPrime				= reshape(vDebtPrime,nProd,nCash);
        mDividends				= reshape(vDividends,nProd,nCash);
        mDebtPriceOptimal       = reshape(vDebtPriceOptimal,nProd,nCash);

        % Compute continuous decisions
        transition_path_continuous_decisions;

        % Store some stuff
        mCapitalPrimeContinuousSeries(:,t) 		 = mCapitalPrimeContinuous(:);
        mDebtPrimeContinuousSeries(:,t) 		 = mDebtPrimeContinuous(:);
        mDividendsContinuousSeries(:,t) 		 = mDividendsContinuous(:);
        mDebtPriceOptimalContinuousSeries(:,t) 	 = mDebtPriceContinuous(:);

        vCapitalPrime     = mCapitalPrimeContinuous(:);
        vDebtPrime        = mDebtPrimeContinuous(:);
        vDividends        = mDividendsContinuous(:);
        vDebtPriceOptimal = mDebtPriceContinuous(:);


        %%%
        % Update value function with continuous decision rules
        %%%

        if option_continuous_optimization_VFI == 1

            % Compute evolution of cash on hand under optimal decisions; CONVENTION: rows = next period productivity, columns = choices (i.e. states)
            mCapitalPrimeOptimal 	 = repmat(vCapitalPrime',[nShocks 1]);
            mDebtPrimeOptimal 	     = repmat(vDebtPrime',[nShocks 1]);
            mProdPrimeOptimal 	     = repmat(mShocksGrid(:,1),[1 nState]);
            mOmegaPrimeOptimal       = repmat(mShocksGrid(:,2),[1 nState]);
            mProfitPrimeOptimal 	 = vA(t+1,1) * (exp(mProdPrimeOptimal) .^ (1 / (1 - nnu))) .* ((exp(mOmegaPrimeOptimal) .* mCapitalPrimeOptimal) .^ tthetaHat) * ...
    								    (vWage(t+1,1) ^ (-nnu / (1 - nnu)));
            mCashPrimeOptimal 	     = mProfitPrimeOptimal + vQ(t+1,1) * (1 - ddelta) * exp(mOmegaPrimeOptimal) .* mCapitalPrimeOptimal - (mDebtPrimeOptimal / vAggregateInflation(t+1,1)) - ...
    									ppsi_0;
            mCashPrimeOptimal 	     = min(max(cashMin * ones(nShocks,nState),mCashPrimeOptimal),cashMax * ones(nShocks,nState));

            % Compute next period's default decisions
            mDefaultExitPrimeOptimal  = (mCashPrimeOptimal <= 0);
            mDefaultContPrimeOptimal  = (mCashPrimeOptimal <= reshape(repmat(mDefaultCutoffSeries(:,t+1),[1 nOmega * nState]),nShocks,nState));

            % Interpolate next period's value function over next period's state
            mValuePrime     = reshape(interpn(mProdGrid,mCashGrid,reshape(mValueSeries(:,t+1),nProd,nCash),...
                              mProdPrimeOptimal(:),mCashPrimeOptimal(:)),nShocks,nState);

            % Compute expected value function next period over exit shock
            mValuePrime     = ppiExit * (1 - mDefaultExitPrimeOptimal) .* mCashPrimeOptimal + ...
                              (1 - ppiExit) * (1 - mDefaultContPrimeOptimal) .* mValuePrime;

            % Computed expected value function over next period's productivity shock
            mEValuePrime1   = reshape(mShocksTransition * mValuePrime,nProd,nProd,nCash);        % has repeats due to z's in both transition matrix and state space
            mEValuePrime    = zeros(nProd,nCash);                   							% will extract repeats here
            for iProd   = 1 : nProd
                mEValuePrime(iProd,:)   = mEValuePrime1(iProd,iProd,:);
            end

            % Update value function
            vValueNew       = vDividends + vSDF(t,1) * mEValuePrime(:);

        end

        %%%
      	% Correct the new decision rules at the boundaries of the state space
      	%%%

        % Unconstrained firms
    	vCapitalPrime(mStateGrid(:,2) >= vUnconstrainedCutoffState)	    = ...
    		vCapitalUnconstrainedState(mStateGrid(:,2) >= vUnconstrainedCutoffState);
    	vDebtPrime(mStateGrid(:,2) >= vUnconstrainedCutoffState)		= ...
    		vDebtUnconstrainedState(mStateGrid(:,2) >= vUnconstrainedCutoffState);
    	vDividends(mStateGrid(:,2) >= vUnconstrainedCutoffState)		= ...
    		vDividendsUnconstrainedState(mStateGrid(:,2) >= vUnconstrainedCutoffState);
    	vDebtPriceOptimal(mStateGrid(:,2) >= vUnconstrainedCutoffState)  = vSDF(t,1);

		% Defaulting firms
		vValueNew(mStateGrid(:,2) <= vDefaultCutoffState)			 = 0;
		vCapitalPrime(mStateGrid(:,2) <= vDefaultCutoffState)		 = capitalMin;
		vDebtPrime(mStateGrid(:,2) <= vDefaultCutoffState)			 = debtMin;
		vDividends(mStateGrid(:,2) <= vDefaultCutoffState)			 = 0;

    end

    %%%
  	% Record results
  	%%%

    mValueSeries(:,t)               = vValueNew;
    mCapitalPrimeSeries(:,t)        = vCapitalPrime;
    mDebtPrimeSeries(:,t)           = vDebtPrime;
    mDividendsSeries(:,t)           = vDividends;
    mDebtPriceOptimalSeries(:,t)    = vDebtPriceOptimal;

end
