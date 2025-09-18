%----------------------------------------------------------------
% Simulate panel of firms to compute annual statistics and correct for selection
%----------------------------------------------------------------

% Panel for tracking state variables
mStatePanel				             = zeros(N,3);		% column 1 = productivity; column 2 = capital; column 3 = cash

% Exogenous productivity shocks
mProductivityShocks		         = randn(N,T_panel);
csvwrite('mProductivityShocks.csv',mProductivityShocks);
mProductivityShocks 	             = csvread('mProductivityShocks.csv');

% Exogenous exit shocks
mExitShocks				         = rand(N,T_panel);
csvwrite('mExitShocks.csv',mExitShocks);
mExitShocks			                 = csvread('mExitShocks.csv');

% Initial productivity draws
vInitialProductivity	             = randn(N,1);
csvwrite('vInitialProductivity.csv',vInitialProductivity);
vInitialProductivity	             = csvread('vInitialProductivity.csv');

% Capital quality shock
mCapitalQuality			         = -((ssigmaOmega ^ 2) / 2) + ssigmaOmega * randn(N,T_panel);
mCapitalQuality			         = ssigmaOmega * randn(N,T_panel);
csvwrite('mCapitalQuality.csv',mCapitalQuality);
mCapitalQuality		                 = csvread('mCapitalQuality.csv');
mCapitalQuality(mCapitalQuality > 0) = 0;

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
mCapitalAdjustedPanel        = zeros(N,T_panel);

% Initialize state variables
mStatePanel(:,1)		     = mmuEnt + ssigmaEnt * vInitialProductivity;
mStatePanel(:,2)		     = k0;
mStatePanel(:,3)		     = b0;


%%%
% Simulate the panel
%%%

for t = 1 : T_panel

	% Ensure that state variables are within grid bounds
	mStatePanel(:,1)		= min(max(prodMin * ones(N,1),mStatePanel(:,1)),prodMax*ones(N,1));
	mStatePanel(:,2)		= min(max(capitalMin * ones(N,1),mStatePanel(:,2)),capitalMax*ones(N,1));
	mStatePanel(:,3)		= min(max(debtMin * ones(N,1),mStatePanel(:,3)),debtMax*ones(N,1));

	%%%
	% Interpolate decisions over state space implied by panel
	%%%

	% Compute implied cash on hand
	vProfitsPanel		     = A * (exp(mStatePanel(:,1)) .^ (1 / (1 - nnu))) .* ((exp(mCapitalQuality(:,t)) .* mStatePanel(:,2)) .^ tthetaHat) .* ...
								(wage ^ (-nnu / (1 - nnu)));
	vCashPanel			     = vProfitsPanel + qSS * (1 - ddelta) * exp(mCapitalQuality(:,t)) .* mStatePanel(:,2) - mStatePanel(:,3) - ...
								ppsi_0;

	% Ensure that implied productivity and cash are within bounds of state space
	vProdPanel			     = min(max(prodMin * ones(N,1),mStatePanel(:,1)),prodMax*ones(N,1));
	vCashPanel			     = min(max(cashMin * ones(N,1),vCashPanel),cashMax*ones(N,1));

	% Interpolate
	vDefaultCutoffPanel	     = interpn(vProdGrid,vDefaultCutoff,vProdPanel);
	vCapitalPrimePanel	     = interpn(mProdGrid,mCashGrid,mCapitalPrime,vProdPanel,vCashPanel);
	vDebtPrimePanel			 = interpn(mProdGrid,mCashGrid,mDebtPrime,vProdPanel,vCashPanel);
	vMarketPricePanel		 = interpn(mProdGrid,mCashGrid,mValue,vProdPanel,vCashPanel);
	vDividendsPanel			 = interpn(mProdGrid,mCashGrid,mDividends,vProdPanel,vCashPanel);

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
	mInvestmentPanel(:,t)		= vCapitalPrimePanel - (1 - ddelta) * exp(mCapitalQuality(:,t)) .* mStatePanel(:,2);
	mCapitalPanel(:,t)			= mStatePanel(:,2);
	mCapitalAdjustedPanel(:,t)	= exp(mCapitalQuality(:,t)) .* mStatePanel(:,2);

	% Market value
	mMarketValuePanel(:,t)		= vMarketPricePanel + vDividendsPanel;

	% Cash flow
	mCashFlowPanel(:,t)			= vProfitsPanel;

	% Internal resources
	mCashPanel(:,t)				= vCashPanel;

  	% Debt
  	mDebtPanel(:,t)         	= mStatePanel(:,3);


	%%%
	% Update the state variables
	%%%

	mStatePanel(:,1)		= rrhoProd * mStatePanel(:,1) + ssigmaProd * mProductivityShocks(:,t);
	mStatePanel(:,2)		= vCapitalPrimePanel;
	mStatePanel(:,3)		= vDebtPrimePanel / inflationSS;

end


%%%
% Construct variables of interest
%%%

% Indicator for whether in the sample
mInSample					 = ones(N,T_panel);
mInSampleDefault	         = ones(N,T_panel);

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
mQuarterlyPanel				= [mFirmIDPanel(:) mQuarterIDPanel(:) mInSample(:) mBalancedPanel(:) mInvestmentPanel(:) ...
								mCashFlowPanel(:) mMarketValuePanel(:) mCapitalPanel(:) mCashPanel(:) mDebtPanel(:) mEmploymentPanel(:)];
csvwrite('mQuarterlyPanel.csv',mQuarterlyPanel);

% Pre allocate matrices for time aggregation
mInvestmentRateAnnual		         = zeros(N,tAnnual);
mExitAnnual							 = zeros(N,tAnnual);
mEmploymentGrowthAnnual	             = zeros(N,tAnnual);
mInSampleAnnual					     = zeros(N,tAnnual);
mInSampleAnnualDefault	             = zeros(N,tAnnual);
mLeverageAnnual					     = zeros(N,tAnnual);

% Do the time aggregation
for t	= 1 : tAnnual

	% Investment rate
	mInvestmentRateAnnual(:,t)	= sum(mInvestmentPanel(:,4*(t-1)+1:4*t),2) ./ mCapitalPanel(:,4*t);

	% Leverage
	mLeverageAnnual(:,t)		= mDebtPanel(:,4*t) ./ mCapitalPanel(:,4*t);

	% Exit
	mExitAnnual(:,t)			= min(ones(N,1),sum(mDefaultPanel(:,4*(t-1)+1:4*t) + mExitPanel(:,4*(t-1)+1:4*t),2));

	% Employment growth
	mEmploymentGrowthAnnual(:,t)= (mEmploymentPanel(:,4*t) - mEmploymentPanel(:,4*(t-1)+1)) ./ ...
									(.5 * (mEmploymentPanel(:,4*t) + mEmploymentPanel(:,4*(t-1)+1)));

	% In sample
	mInSampleAnnual(:,t)		= ((min(4*ones(N,1),sum(mInSample(:,4*(t-1)+1:4*t),2))/4) >= .8); % just some number greater than .75

	% In sample (including default)
	mInSampleAnnualDefault(:,t)	= sum(mInSampleDefault(:,4*(t-1)+1:4*t),2) > 0;

end

%%%
% Investment rate moments
%%%

% Indicator for being a member of the balanced panel
mBalancedPanelIndicator				= logical((mInSampleAnnual(:,tAnnual) == 1) * ones(1,tAnnual));
mBalancedPanelIndicator(:,1:tBurn)	= 0;
nBalancedPanel						= sum(mBalancedPanelIndicator(:,end));

% Investment rate statistics
vInvestmentRateBalancedPanel	= mInvestmentRateAnnual(mBalancedPanelIndicator);
meanIKAnnual					= mean(vInvestmentRateBalancedPanel);
sdIKAnnual						= std(vInvestmentRateBalancedPanel);

%%%
% Lifecycle moments
%%%

% Employment growth rates
vFirmsInSample			= sum(mInSampleAnnual,1)';
vGrowthRates			= zeros(tAnnual,1);
for t = 1 : tAnnual
	vGrowthRates(t,1)	= sum(mEmploymentGrowthAnnual(logical(mInSampleAnnual(:,t)),t)) / vFirmsInSample(t,1);
end


%%%
% Distribution of growth rates, as in Davis et al.
%%%

% Compute employment growth
mEmploymentGrowth	= (mEmploymentPanel(:,2:end) - mEmploymentPanel(:,1:end-1)) ./ ...
						(.5 * (mEmploymentPanel(:,2:end) + mEmploymentPanel(:,1:end-1)));
mEmploymentGrowth 	= [mEmploymentGrowth  mEmploymentGrowth(:,end)];

% Categories of employment growth rates
vCategories			= [-2;-.2;-.05;-.02;0;.02;.05;.2;2];
nCat 				= 9;
N_firms 			= sum(logical(mInSample(:)));
vMeanGrowth 		= zeros(nCat-1,1);
for iCat = 1 : nCat-1
	X 				=	(mEmploymentGrowth(logical(mInSample)) > vCategories(iCat)) & (mEmploymentGrowth(logical(mInSample)) < vCategories(iCat+1));
	vMeanGrowth(iCat,1)	= sum(X(:)) / N_firms;
end
