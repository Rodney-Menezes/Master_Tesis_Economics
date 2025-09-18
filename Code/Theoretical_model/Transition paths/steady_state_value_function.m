%---------------------------------------------------------------
% Pre-compute some useful objects
%---------------------------------------------------------------

% Compute all possible choices of dividends; CONVENTION: rows = state variables, columns = choices
mDebtPriceChoices   = reshape(repmat(reshape(mDebtPrice,nProd,1,nChoices),[1 nCash 1]),...
                        nState,nChoices);			% define debt price over entire state space
mDividendsChoices   = repmat(mStateGrid(:,2),[1 nChoices]) + mDebtPriceChoices .* repmat(mChoicesGrid(:,2)',[nState 1]) - ...
                        qSS * repmat(mChoicesGrid(:,1)',[nState 1]);

% Compute the unconstrained value function and choices over the whole state space (i.e., defined over productivity and cash on hand)
vValueUnconstrainedState    = reshape(repmat(vValueUnconstrained,[1 nCash]),nState,1) + ...
                                mStateGrid(:,2);    % continuation value of a firm who becomes unconstrained
vCapitalUnconstrainedState  = reshape(repmat(vCapitalUnconstrained,[1 nCash]),nState,1);    % capital choice
vDebtUnconstrainedState     = reshape(repmat(vDebtUnconstrained,[1 nCash]),nState,1);       % debt choice
vDividendsUnconstrainedState= mStateGrid(:,2) + bbeta * vDebtUnconstrainedState - ...
                                qSS * vCapitalUnconstrainedState;                           % dividend payment

% Compute thresholds over the whole state space
vUnconstrainedCutoffState   = reshape(repmat(vUnconstrainedCutoff,[1 nCash]),nState,1);     % cutoff for being unconstrained
vDefaultCutoffState         = reshape(repmat(vDefaultCutoff,[1 nCash]),nState,1);           % cutoff for defaulting

% Compute evolution of cash on hand for computing expectations; CONVENTION: rows = productivity next period; columns = choices
% mCashPrime computed above for the debt price
mCashPrime          = min(max(cashMin * ones(nShocks,nChoices),mCashPrime),...
                        cashMax * ones(nShocks,nChoices)); 		% ensures it stays inside bounds of state space; other interpolation gives error
vCashPrime          = mCashPrime(:);                     		% collapse into large vector for interpolation
vProdPrime          = reshape(repmat(mShocksGrid(:,1),[1 nChoices]),nShocks * nChoices,1);   % corresponding vector of productivity shocks
mDefaultExitPrime   = reshape((vCashPrime <= 0),nShocks,nChoices);   % indicator for default, conditional on exit shock
mDefaultContPrime   = reshape((vCashPrime <= reshape(repmat(vDefaultCutoff,[1 nOmega * nChoices]),...
                        nShocks * nChoices,1)),nShocks,nChoices);      % indicator for default, conditional on no exit shock


%---------------------------------------------------------------
% Discrete value function iteration
%---------------------------------------------------------------

% Initialize variables for the iteration
vValue                 = vValueUnconstrainedState; 	            % initialize with unconstrained value function
iteration              = 1;                        		        % current iteration
tolerance              = 1e-10;                    				% acceptable error before terminating iteration
maxDifference          = 100;                      				% maximum difference between iterations over the state space
dampening              = 0;                        				% weight placed on old iteration in updating
maxIterations          = 30;                      				% maximum number of iterations before terminating
vErr                   = zeros(maxIterations,1);   		        % store error at each step of the iteration before debugging
accStart               = 5;	                        			% overall iteration in which to begin policy function iteration
accIterations          = 100;                      				% number of times to update the value function for the same policy function
continuousStart        = 10;                       				% iteration on which to begin continuous optimization

% Do the iteration
while maxDifference > tolerance && iteration <= maxIterations

    if (option_continuous_optimization_VFI == 0) || (iteration < continuousStart && option_continuous_optimization_VFI)

        %%%
        % Compute current period's decisions by discrete value function iteration
        %%%

        % Interpolate next period's value function over next period's possible states
        mValuePrime     = reshape(interpn(mProdGrid,mCashGrid,reshape(vValue,nProd,nCash),...
        				    vProdPrime,vCashPrime),nShocks,nChoices); % rows = next period's productivity, columns = choices of k' and b'

        % Compute expected value function next period over exit shock
        mValuePrime     = ppiExit * (1 - mDefaultExitPrime) .* mCashPrime + (1 - ppiExit) * (1 - mDefaultContPrime) .* mValuePrime;

        % Compute expected value function over next period's productivity shock
        mEValuePrime    = mShocksTransition * mValuePrime;
        mEValuePrime    = reshape(repmat(reshape(mEValuePrime,nProd,1,nChoices),...
        			      [1 nCash 1]),nState,nChoices);  				% expand so that rows = state variables, columns = choices

        % Compute objective function for the firm
        mObjective                          = mDividendsChoices + bbeta * mEValuePrime;
        mObjective(mDividendsChoices < 0)   = -1e4;             			% enforce non-negativity constraint on dividends

        % Do the optimization
        [vValueNew,vChoiceIndices]          = max(mObjective,[],2);

        % Extract policy functions
        vCapitalPrime      = mChoicesGrid(vChoiceIndices,1);   % capital policy function
        vDebtPrime         = mChoicesGrid(vChoiceIndices,2);   % debt policy function
        vDividends         = zeros(nState,1);                  % dividend policy function
        vDebtPriceOptimal  = zeros(nState,1);                  % implied debt price
        for iState  = 1 : nState
        	vDividends(iState,1)        = mDividendsChoices(iState,vChoiceIndices(iState,1));
        	vDebtPriceOptimal(iState,1) = mDebtPriceChoices(iState,vChoiceIndices(iState,1));
        end

    else

        % Rename and reshape some stuff
		mValue					 = reshape(vValue,nProd,nCash);
		mCapitalPrime			 = reshape(vCapitalPrime,nProd,nCash);
		mDebtPrime				 = reshape(vDebtPrime,nProd,nCash);
		mDividends				 = reshape(vDividends,nProd,nCash);
        mDebtPriceOptimal        = reshape(vDebtPriceOptimal,nProd,nCash);

		% Compute the continuous decisions
		steady_state_continuous_decisions;

		% Rename decisions
		vCapitalPrime 		     = mCapitalPrimeContinuous(:);
		vDebtPrime 		         = mDebtPrimeContinuous(:);
		vDividends 		         = mDividendsContinuous(:);
		vDebtPriceOptimal 		 = mDebtPriceContinuous(:);

		%%%
		% Update value function with continuous decision rules
		%%%

		% Compute evolution of cash on hand under optimal decisions; CONVENTION: rows = next period shocks, columns = choices (i.e. states)
		mCapitalPrimeOptimal    = repmat(vCapitalPrime',[nShocks 1]);
		mDebtPrimeOptimal       = repmat(vDebtPrime',[nShocks 1]);
		mProdPrimeOptimal       = repmat(mShocksGrid(:,1),[1 nState]);
		mOmegaPrimeOptimal		= repmat(mShocksGrid(:,2),[1 nState]);
		mProfitPrimeOptimal     = A * (exp(mProdPrimeOptimal) .^ (1 / (1 - nnu))) .* ((exp(mOmegaPrimeOptimal) .* mCapitalPrimeOptimal) .^ tthetaHat) * ...
									(wage ^ (-nnu / (1 - nnu)));
		mCashPrimeOptimal       = mProfitPrimeOptimal + qSS * (1 - ddelta) * exp(mOmegaPrimeOptimal) .* mCapitalPrimeOptimal - (mDebtPrimeOptimal / inflationSS) - ...
									ppsi_0;
		mCashPrimeOptimal       = min(max(cashMin * ones(nShocks,nState),mCashPrimeOptimal),cashMax * ones(nShocks,nState));

		% Compute next period's default decisions
		mDefaultExitPrimeOptimal    = (mCashPrimeOptimal <= 0);
		mDefaultContPrimeOptimal    = (mCashPrimeOptimal <= reshape(repmat(vDefaultCutoff,[1 nOmega * nState]),nShocks,nState));

		% Interpolate next period's value function over next period's state
		mValuePrime     = reshape(interpn(mProdGrid,mCashGrid,mValue,...
				            mProdPrimeOptimal(:),mCashPrimeOptimal(:)),nShocks,nState);

		% Compute expected value function next period over exit shock
		mValuePrime     = ppiExit * (1 - mDefaultExitPrimeOptimal) .* mCashPrimeOptimal + ...
							(1 - ppiExit) * (1 - mDefaultContPrimeOptimal) .* mValuePrime;

		% Computed expected value function over next period's productivity shock
		mEValuePrime1   = reshape(mShocksTransition * mValuePrime,nProd,nProd,nCash);       % has repeats due to z's in both transition matrix and state space
		mEValuePrime    = zeros(nProd,nCash);                   															% will extract repeats here
		for iProd   = 1 : nProd
			mEValuePrime(iProd,:)   = mEValuePrime1(iProd,iProd,:);
		end

		% Update value function
		vValueNew       = vDividends + bbeta * mEValuePrime(:);

    end


    %%%
    % Correct the unconstrained firms
    %%%

    vValueNew(mStateGrid(:,2) >= vUnconstrainedCutoffState)		   = ...
    	vValueUnconstrainedState(mStateGrid(:,2) >= vUnconstrainedCutoffState);
    vCapitalPrime(mStateGrid(:,2) >= vUnconstrainedCutoffState)	   = ...
    	vCapitalUnconstrainedState(mStateGrid(:,2) >= vUnconstrainedCutoffState);
    vDebtPrime(mStateGrid(:,2) >= vUnconstrainedCutoffState)	   = ...
    	vDebtUnconstrainedState(mStateGrid(:,2) >= vUnconstrainedCutoffState);
    vDividends(mStateGrid(:,2) >= vUnconstrainedCutoffState)	   = ...
    	vDividendsUnconstrainedState(mStateGrid(:,2) >= vUnconstrainedCutoffState);
    vDebtPriceOptimal(mStateGrid(:,2) >= vUnconstrainedCutoffState)= bbeta;


    %%%
    % Correct defaulting firms
    %%%

    vValueNew(mStateGrid(:,2) <= vDefaultCutoffState)		= 0;
    vCapitalPrime(mStateGrid(:,2) <= vDefaultCutoffState)	= capitalMin;     % a convention
    vDebtPrime(mStateGrid(:,2) <= vDefaultCutoffState)		= debtMin;
    vDividends(mStateGrid(:,2) <= vDefaultCutoffState)		= 0;


    %%%
    % Update iteration
    %%%

    maxDifference		= max(abs(vValueNew - vValue));
    vErr(iteration)		= maxDifference;
    vValue			    = (1 - dampening) * vValueNew + dampening * vValue;
    iteration			= iteration + 1;


    %%%
    % Howard's improvement step (update value function many times with same policy function)
    %%%

    if iteration > accStart && accIterations >= 1

        % Compute evolution of cash on hand under optimal decisions; CONVENTION: rows = next period shocks, columns = choices (i.e. states)
        mCapitalPrimeOptimal    = repmat(vCapitalPrime',[nShocks 1]);
        mDebtPrimeOptimal       = repmat(vDebtPrime',[nShocks 1]);
        mProdPrimeOptimal       = repmat(mShocksGrid(:,1),[1 nState]);
        mOmegaPrimeOptimal		= repmat(mShocksGrid(:,2),[1 nState]);
        mProfitPrimeOptimal     = A * (exp(mProdPrimeOptimal) .^ (1 / (1 - nnu))) .* ((exp(mOmegaPrimeOptimal) .* mCapitalPrimeOptimal) .^ tthetaHat) * ...
        						     (wage ^ (-nnu / (1 - nnu)));
        mCashPrimeOptimal       = mProfitPrimeOptimal + qSS * (1 - ddelta) * exp(mOmegaPrimeOptimal) .* mCapitalPrimeOptimal - (mDebtPrimeOptimal / inflationSS) - ...
        							  ppsi_0;
        mCashPrimeOptimal       = min(max(cashMin * ones(nShocks,nState),mCashPrimeOptimal),cashMax * ones(nShocks,nState));

        % Compute next period's default decisions
        mDefaultExitPrimeOptimal    = (mCashPrimeOptimal <= 0);
        mDefaultContPrimeOptimal    = (mCashPrimeOptimal <= reshape(repmat(vDefaultCutoff,[1 nOmega * nState]),nShocks,nState));

        % Initialize the iteration
        iterationAcc            	= 1;
        maxDifferenceAcc        	= 100;

        % Do the iteration
        while iterationAcc <= accIterations & maxDifferenceAcc >= tolerance

            % Interpolate next period's value function over next period's state
            mValuePrime     = reshape(interpn(mProdGrid,mCashGrid,reshape(vValue,nProd,nCash),...
                              mProdPrimeOptimal(:),mCashPrimeOptimal(:)),nShocks,nState);

            % Compute expected value function next period over exit shock
            mValuePrime     = ppiExit * (1 - mDefaultExitPrimeOptimal) .* mCashPrimeOptimal + ...
                                (1 - ppiExit) * (1 - mDefaultContPrimeOptimal) .* mValuePrime;

            % Computed expected value function over next period's productivity shock
            mEValuePrime1   = reshape(mShocksTransition * mValuePrime,nProd,nProd,nCash);        % has repeats due to z's in both transition matrix and state space
            mEValuePrime    = zeros(nProd,nCash);                   							 % will extract repeats here
            for iProd   = 1 : nProd
                mEValuePrime(iProd,:)   = mEValuePrime1(iProd,iProd,:);
            end

            % Update value function
            vValueNew       = vDividends + bbeta * mEValuePrime(:);

            % Correct unconstrained firms
            vValueNew(mStateGrid(:,2) >= vUnconstrainedCutoffState)		= ...
                vValueUnconstrainedState(mStateGrid(:,2) >= vUnconstrainedCutoffState);

            % Correct defaulting firms
            vValueNew(mStateGrid(:,2) <= vDefaultCutoffState)			= 0;

            % Update Howard improvement iteration
            maxDifferenceAcc    = max(abs(vValueNew - vValue));
            iterationAcc        = iterationAcc + 1;
            vValue              = (1 - dampening) * vValueNew + dampening * vValue;

        end

    end             % end Howard's improvement step here

end             % end VFI here

% Compute evolution of cash on hand under optimal decisions if Howard's improvement step is turned off
if iteration <= accStart

	% Compute evolution of cash on hand under optimal decisions; CONVENTION: rows = next period shocks, columns = choices (i.e. states)
	mCapitalPrimeOptimal    = repmat(vCapitalPrime',[nShocks 1]);
	mDebtPrimeOptimal       = repmat(vDebtPrime',[nShocks 1]);
	mProdPrimeOptimal       = repmat(mShocksGrid(:,1),[1 nState]);
	mOmegaPrimeOptimal		= repmat(mShocksGrid(:,2),[1 nState]);
	mProfitPrimeOptimal     = A * (exp(mProdPrimeOptimal) .^ (1 / (1 - nnu))) .* ((exp(mOmegaPrimeOptimal) .* mCapitalPrimeOptimal) .^ tthetaHat) * ...
							     (wage ^ (-nnu / (1 - nnu)));
	mCashPrimeOptimal       = mProfitPrimeOptimal + qSS * (1 - ddelta) * exp(mOmegaPrimeOptimal) .* mCapitalPrimeOptimal - (mDebtPrimeOptimal / inflationSS) - ...
								 ppsi_0;
	mCashPrimeOptimal       = min(max(cashMin * ones(nShocks,nState),mCashPrimeOptimal),cashMax * ones(nShocks,nState));

	% Compute next period's default decisions
	mDefaultExitPrimeOptimal= (mCashPrimeOptimal <= 0);
	mDefaultContPrimeOptimal= (mCashPrimeOptimal <= reshape(repmat(vDefaultCutoff,[1 nOmega * nState]),nShocks,nState));

end

%----------------------------------------------------------------
% Save output of the iterations
%----------------------------------------------------------------

% Save objects
mValue					= reshape(vValue,nProd,nCash);
mCapitalPrime			= reshape(vCapitalPrime,nProd,nCash);
mDebtPrime				= reshape(vDebtPrime,nProd,nCash);
mDividends				= reshape(vDividends,nProd,nCash);
mCashPrime				= reshape(mCashPrimeOptimal,nShocks,nProd,nCash);
mDebtPriceOptimal	    = reshape(vDebtPriceOptimal,nProd,nCash);
