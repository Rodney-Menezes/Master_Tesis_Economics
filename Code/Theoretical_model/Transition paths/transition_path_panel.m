%----------------------------------------------------------------
% Set parameters for the simulation
%----------------------------------------------------------------

N						= 50000;
tSim				    = 10;
T_panel					= tPre + tSim;


%----------------------------------------------------------------
% Pre-allocate matrices for storing results
%----------------------------------------------------------------

% Panel for tracking state variables
mStatePanel				= zeros(N,3);		% column 1 = productivity; column 2 = capital; column 3 = debt

% Exogenous productivity shocks
mProductivityShocks		= randn(N,T_panel);
csvwrite('mProductivityShocks_transition.csv',mProductivityShocks);
%mProductivityShocks 	= csvread('mProductivityShocks_transition.csv');

% Exogenous exit shocks
mExitShocks				= rand(N,T_panel);
csvwrite('mExitShocks_transition.csv',mExitShocks);
%mExitShocks			= csvread('mExitShocks_transition.csv');

% Initial productivity draws
vInitialProductivity	= randn(N,1);
csvwrite('vInitialProductivity_transition.csv',vInitialProductivity);
%vInitialProductivity	= csvread('vInitialProductivity_transition.csv');

% Capital quality shock
mCapitalQuality			= ssigmaOmega * randn(N,T_panel);
mCapitalQuality         = min(zeros(size(mCapitalQuality)),mCapitalQuality);
csvwrite('mCapitalQuality_transition.csv',mCapitalQuality);
%mCapitalQuality		= csvread('mCapitalQuality_transition.csv');


% Matrices to store results
mEmploymentPanel		     = zeros(N,T_panel);
mInvestmentPanel		     = zeros(N,T_panel);
mCapitalPanel				 = zeros(N,T_panel);
mDefaultPanel				 = zeros(N,T_panel);
mExitPanel					 = zeros(N,T_panel);
mMarketValuePanel		     = zeros(N,T_panel);
mCashFlowPanel			     = zeros(N,T_panel);
mCashPanel					 = zeros(N,T_panel);
mDebtPanel                   = zeros(N,T_panel);
mCapitalAdjustedPanel	     = zeros(N,T_panel);
mInterestRatePanel           = zeros(N,T_panel);
mUnconstrainedPanel		     = zeros(N,T_panel);

% Initialize state variables
mStatePanel(:,1)		     = mmuEnt + ssigmaEnt * vInitialProductivity;
mStatePanel(:,2)		     = k0;
mStatePanel(:,3)		     = b0;


%----------------------------------------------------------------
% Simulate the panel in steady state
%----------------------------------------------------------------

for t = 1 : tPre

	%%%
	% Interpolate decisions over state space implied by panel
	%%%

	% Compute implied cash on hand
	vProfitsPanel		=  A * (exp(mStatePanel(:,1)) .^ (1 / (1 - nnu))) .* ((exp(mCapitalQuality(:,t)) .* mStatePanel(:,2)) .^ tthetaHat) .* ...
							(wage ^ (-nnu / (1 - nnu)));
	vCashPanel			= vProfitsPanel + qSS * (1 - ddelta) * exp(mCapitalQuality(:,t)) .* mStatePanel(:,2) - mStatePanel(:,3) - ...
							ppsi_0;

	% Ensure that implied productivity and cash are within bounds of state space
	vProdPanel			= min(max(prodMin * ones(N,1),mStatePanel(:,1)),prodMax*ones(N,1));
	vCashPanel			= min(max(cashMin * ones(N,1),vCashPanel),cashMax*ones(N,1));

	% Interpolate
	vDefaultCutoffPanel			= interpn(vProdGrid,vDefaultCutoffSS,vProdPanel);
    vUnconstrainedCutoffPanel   = interpn(vProdGrid,vUnconstrainedCutoffSS,vProdPanel);
	vCapitalPrimePanel			= interpn(mProdGrid,mCashGrid,mCapitalPrimeSS,vProdPanel,vCashPanel);
    vDebtPricePanel         	= interpn(mProdGrid,mCashGrid,mDebtPriceOptimalSS,vProdPanel,vCashPanel);
	vDebtPrimePanel				= interpn(mProdGrid,mCashGrid,mDebtPrimeSS,vProdPanel,vCashPanel);
	vMarketPricePanel			= interpn(mProdGrid,mCashGrid,mValueSS,vProdPanel,vCashPanel);
	vDividendsPanel				= interpn(mProdGrid,mCashGrid,mDividendsSS,vProdPanel,vCashPanel);

	%%%
	% Compute decisions
	%%%

	% Exit shock realization
	mExitPanel(:,t)			= (mExitShocks(:,t) <= ppiExit);

	% Default decision
	mDefaultPanel(:,t)		= mExitPanel(:,t) .* (vCashPanel <= 0) + (1 - mExitPanel(:,t)) .* (vCashPanel <= vDefaultCutoffPanel);

	% Employment
	mEmploymentPanel(:,t)	= ((nnu * pSS * exp(vProdPanel) .* ((exp(mCapitalQuality(:,t)) .* ...
								mStatePanel(:,2)) .^ ttheta)) / wage) .^ (1 / (1 - nnu));

	% Investment and capital
	mInvestmentPanel(:,t)	= vCapitalPrimePanel - (1 - ddelta) * exp(mCapitalQuality(:,t)) .* mStatePanel(:,2);
	mCapitalPanel(:,t)		= mStatePanel(:,2);
	mCapitalAdjustedPanel(:,t)	= exp(mCapitalQuality(:,t)) .* mStatePanel(:,2);

    % Interest rate
    mInterestRatePanel(:,t) = min((1 ./ vDebtPricePanel),100);

	% Market value
	mMarketValuePanel(:,t)	= vMarketPricePanel + vDividendsPanel;

	% Cash flow
	mCashFlowPanel(:,t) 	= vProfitsPanel;

	% Internal resources
	mCashPanel(:,t)			= vCashPanel;

    % Debt
    mDebtPanel(:,t)   		= mStatePanel(:,3);

    % Unconstrained?
    mUnconstrainedPanel(:,t) = (vCashPanel > vUnconstrainedCutoffPanel);


	%%%
	% Update the state variables
	%%%

	mStatePanel(:,1)		= rrhoProd * mStatePanel(:,1) + ssigmaProd * mProductivityShocks(:,t);
	mStatePanel(:,2)		= vCapitalPrimePanel;
	mStatePanel(:,3)		= vDebtPrimePanel / inflationSS;

end


%----------------------------------------------------------------
% Simulate the panel during the transition
%----------------------------------------------------------------

for t = tPre+1 : T_panel

	%%%
	% Interpolate decisions over state space implied by panel
	%%%

	% Compute implied cash on hand
	vProfitsPanel		=  vA(t-tPre,1) * (exp(mStatePanel(:,1)) .^ (1 / (1 - nnu))) .* ((exp(mCapitalQuality(:,t)) .* mStatePanel(:,2)) .^ tthetaHat) .* ...
							(vWage(t-tPre,1) ^ (-nnu / (1 - nnu)));
	vCashPanel			= vProfitsPanel + vQ(t-tPre,1) * (1 - ddelta) * exp(mCapitalQuality(:,t)) .* mStatePanel(:,2) - mStatePanel(:,3) - ...
							ppsi_0;

	% Ensure that implied productivity and cash are within bounds of state space
	vProdPanel			= min(max(prodMin * ones(N,1),mStatePanel(:,1)),prodMax*ones(N,1));
	vCashPanel			= min(max(cashMin * ones(N,1),vCashPanel),cashMax*ones(N,1));

	% Interpolate
	vDefaultCutoffPanel			= interpn(vProdGrid,mDefaultCutoffSeries(:,t-tPre),vProdPanel);
 	vUnconstrainedCutoffPanel	= interpn(vProdGrid,mUnconstrainedCutoffSeries(:,t-tPre),vProdPanel);
	vCapitalPrimePanel			= interpn(mProdGrid,mCashGrid,reshape(mCapitalPrimeSeries(:,t-tPre),nProd,nCash),vProdPanel,vCashPanel);
    vDebtPricePanel         	= interpn(mProdGrid,mCashGrid,reshape(mDebtPriceOptimalSeries(:,t-tPre),nProd,nCash),vProdPanel,vCashPanel);
	vDebtPrimePanel				= interpn(mProdGrid,mCashGrid,reshape(mDebtPrimeSeries(:,t-tPre),nProd,nCash),vProdPanel,vCashPanel);
	vMarketPricePanel			= interpn(mProdGrid,mCashGrid,reshape(mValueSeries(:,t-tPre),nProd,nCash),vProdPanel,vCashPanel);
	vDividendsPanel				= interpn(mProdGrid,mCashGrid,reshape(mDividendsSeries(:,t-tPre),nProd,nCash),vProdPanel,vCashPanel);

	%%%
	% Compute decisions
	%%%

	% Exit shock realization
	mExitPanel(:,t)				= (mExitShocks(:,t) <= ppiExit);

	% Default decision
	mDefaultPanel(:,t)			= mExitPanel(:,t) .* (vCashPanel <= 0) + (1 - mExitPanel(:,t)) .* (vCashPanel <= vDefaultCutoffPanel);

	% Employment
	mEmploymentPanel(:,t)		= ((nnu * vP(t-tPre,1) * exp(vProdPanel) .* ((exp(mCapitalQuality(:,t)) .* ...
									mStatePanel(:,2)) .^ ttheta)) / vWage(t-tPre,1)) .^ (1 / (1 - nnu));

	% Investment and capital
	mInvestmentPanel(:,t)		= vCapitalPrimePanel - (1 - ddelta) * exp(mCapitalQuality(:,t)) .* mStatePanel(:,2);
	mCapitalPanel(:,t)			= mStatePanel(:,2);
	mCapitalAdjustedPanel(:,t)	= exp(mCapitalQuality(:,t)) .* mStatePanel(:,2);

	% Market value
	mMarketValuePanel(:,t)		= vMarketPricePanel + vDividendsPanel;

    % Interest rate
    mInterestRatePanel(:,t) 	= min((1 ./ vDebtPricePanel),100);

	% Cash flow
	mCashFlowPanel(:,t)			= vProfitsPanel;

	% Internal resources
	mCashPanel(:,t)				= vCashPanel;

    % Debt
    mDebtPanel(:,t)         	= mStatePanel(:,3);

    % Unconstrained?
    mUnconstrainedPanel(:,t) 	= (vCashPanel > vUnconstrainedCutoffPanel);

	%%%
	% Update the state variables
	%%%

	mStatePanel(:,1)		= rrhoProd * mStatePanel(:,1) + ssigmaProd * mProductivityShocks(:,t);
	mStatePanel(:,2)		= vCapitalPrimePanel;
	mStatePanel(:,3)		= vDebtPrimePanel / vAggregateInflation(t-tPre+1,1);

end


%----------------------------------------------------------------
% Clean up panel for analysis
%----------------------------------------------------------------

% Indicator for whether in the sample
mInSample			= ones(N,T_panel);
mInSampleDefault	= ones(N,T_panel);

for t = 1 : T_panel

	% Not including default
	mInSample(mExitPanel(:,t) == 1,t+1:end)	= 0;
	mInSample(mDefaultPanel(:,t) == 1,t:end)= 0;

	% Including default
	mInSampleDefault(mExitPanel(:,t) == 1,t+1:end)	= 0;
	mInSampleDefault(mDefaultPanel(:,t) == 1,t+1:end)= 0;

end

% Save quarterly variables in .csv file for analysis in STATA
mFirmIDPanel				= linspace(1,N,N)' * ones(1,T_panel);
mQuarterIDPanel				= ones(N,1) * linspace(1,T_panel,T_panel);
mBalancedPanel				= logical((mInSample(:,T_panel) == 1) * ones(1,T_panel));
mTransitionPanel			= [mFirmIDPanel(:) mQuarterIDPanel(:) mInSample(:) mBalancedPanel(:) mInvestmentPanel(:) ...
								mCashFlowPanel(:) mMarketValuePanel(:) mCapitalPanel(:) mCashPanel(:) ...
								mDebtPanel(:) mEmploymentPanel(:) mInterestRatePanel(:) mUnconstrainedPanel(:)];
