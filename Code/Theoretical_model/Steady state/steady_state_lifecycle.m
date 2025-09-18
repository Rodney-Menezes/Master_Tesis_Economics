%----------------------------------------------------------------
% Do the simulation
%----------------------------------------------------------------

% Preallocate matrices to store results
vMass				 = zeros(T_lifecycle,1);
vMassContinue		 = zeros(T_lifecycle,1);
vAvgProductivity     = zeros(T_lifecycle,1);
vAvgCash			 = zeros(T_lifecycle,1);
vCapital             = zeros(T_lifecycle,1);
vDebt                = zeros(T_lifecycle,1);
vOutput              = zeros(T_lifecycle,1);
vLabor               = zeros(T_lifecycle,1);
vMass                = zeros(T_lifecycle,1);
vAvgDebtPrice        = zeros(T_lifecycle,1);
vAvgGrossLev         = zeros(T_lifecycle,1);
vAvgNetLev           = zeros(T_lifecycle,1);
vFracPosDebt         = zeros(T_lifecycle,1);

% Compute debt price over distribution grid
mDebtPriceDistribution  = reshape(interpn(mProdGrid,mCashGrid,reshape(vDebtPriceOptimal,nProd,nCash),...
                            mStateGridDist(:,1),mStateGridDist(:,2)),nProd,nCashDist);

% Initialize distributions used in the computation
mDistribution_lifecycle_previous	 = reshape(vDistributionEntrants,nProd,nCashDist);
mDistribution_lifecycle				 = reshape(vDistributionEntrants,nProd,nCashDist);

% Initial capital and debt
vCapital(1,1)			= k0;
vDebt(1,1)              = b0;
vAvgGrossLev(1,1)       = b0 / k0;
vAvgNetLev(1,1)         = b0 / k0;
vFracPosDebt(1,1)       = 0;

% Integrands for t = 1
vLaborIntegrand			= (((pSS * nnu * exp(mShocksGrid(:,1)) .* ((exp(mShocksGrid(:,2)) * k0) .^ ttheta)) / ...
							       wage) .^ (1 / (1 - nnu)) .* vContinueEntrant);
vOutputIntegrand		= (exp(mShocksGrid(:,1)) .* ((exp(mShocksGrid(:,2)) * k0) .^ ttheta) .* ...
							       (vLaborIntegrand .^ nnu)) .* vContinueEntrant;
vLaborIntegrand			= sum(reshape(vLaborIntegrand,nProd,nOmega) .* repmat(vOmegaWeights',[nProd 1]),2);
vOutputIntegrand		= sum(reshape(vOutputIntegrand,nProd,nOmega) .* repmat(vOmegaWeights',[nProd 1]),2);


% Integrands for t > 1
aLaborIntegrand     = reshape((((pSS * nnu * exp(mProdPrimeDist) .* ((exp(mOmegaPrimeDist) .* mCapitalPrimePrimeDist) .^ ttheta)) / ...
							         wage) .^ (1 / (1 - nnu))) .* mContinueIncumbent,nShocks,nProd,nCashDist);
aOutputIntegrand    = reshape((exp(mProdPrimeDist) .* ((exp(mOmegaPrimeDist) .* mCapitalPrimePrimeDist) .^ ttheta) .* ...
							       (reshape(aLaborIntegrand,nShocks,nStateDist) .^ nnu)) .* ...
                                    mContinueIncumbent,nShocks,nProd,nCashDist);
mLaborIntegrand     = zeros(nProd,nCashDist);
mOutputIntegrand    = zeros(nProd,nCashDist);
for iProd = 1 : nProd
    mLaborIntegrand(iProd,:)     = mShocksTransition(iProd,:) * squeeze(aLaborIntegrand(:,iProd,:));
    mOutputIntegrand(iProd,:)    = mShocksTransition(iProd,:) * squeeze(aOutputIntegrand(:,iProd,:));
end

%%%
% Compute lifecycle simulation
%%%

for t = 1 : T_lifecycle

	%%%
	% Compute distribution of firms in production
	%%%

	mDistributionProduction_lifecycle			= reshape(vContinue .* mDistribution_lifecycle(:),nProd,nCashDist);

	%%%
	% Compute objects from this period's distribution
	%%%

	% Mass of firms in production
    vMass(t,1)				= sum(mDistributionProduction_lifecycle(:));

	% Average cash
	vAvgCash(t,1)			= sum(mCashGridDist(:) .* mDistributionProduction_lifecycle(:));

	% Average productivity
	vAvgProductivity(t,1)	= sum(exp(mProdGridDist(:)) .* mDistributionProduction_lifecycle(:));

	% Average debt price
	vAvgDebtPrice(t,1)		= (1 - ppiExit) * sum(vDebtPriceProduction .* vDistContContinue .* mDistribution_lifecycle(:)) ./ ...
							     ((1 - ppiExit) * sum(vDistContContinue .* mDistribution_lifecycle(:)));	% note it is already normalized

	% Next period's capital
    vCapital(t+1,1)         = (1 - ppiExit) * sum(vCapitalPrimeDist .* vDistContContinue .* ...
							     mDistribution_lifecycle(:));

	% Next period's debt
	vDebt(t+1,1)		 	= (1 - ppiExit) * sum(vDebtPrimeDist .* vDistContContinue .* ...
								 mDistribution_lifecycle(:));

    % Fraction with positive debt
    vFracPosDebt(t+1,1)     = (1 - ppiExit) * sum((vDebtPrimeDist > 0) .* vDistContContinue .* ...
								 mDistribution_lifecycle(:));

    % Moments of gross leverage
    vAvgGrossLev(t+1,1)     = (1 - ppiExit) * sum(vGrossLeverageDist .* vDistContContinue .* ...
								 mDistribution_lifecycle(:));

    % Moments of net leverage
    vAvgNetLev(t+1,1)       = (1 - ppiExit) * sum(vNetLeverageDist .* vDistContContinue .* ...
								 mDistribution_lifecycle(:));

	%%%
	% Compute objects from last period's distribution
	%%%

	% In first period, distribution comes from new entrants
	if t == 1

		vLabor(t,1)			= sum(vLaborIntegrand .* vProdDistEnt);
		vOutput(t,1)		= sum(vOutputIntegrand .* vProdDistEnt);

	% In later periods, distribution comes from last period's incumbents
	elseif t > 1

        vLabor(t,1)   = (1 - ppiExit) * sum(mLaborIntegrand(:) .* mDistribution_lifecycle_previous(:) .* ...
                                mDistContContinue(:));
        vOutput(t,1)  = (1 - ppiExit) * sum(mOutputIntegrand(:) .* mDistribution_lifecycle_previous(:) .* ...
                                mDistContContinue(:));

	end

	%%%
	% Update distribution
	%%%

	mDistribution_lifecycle_previous		= mDistribution_lifecycle;
	mDistribution_lifecycle					= reshape((1 - ppiExit) * mTransitionIncumbents * ...
                                                (vDistContContinue .* mDistribution_lifecycle(:)),nProd,nCashDist);

end


%----------------------------------------------------------------
% Normalize series for plotting
%----------------------------------------------------------------

%%%
% Aggregate series
%%%

vCapitalAggregate       	= vCapital(1:T_lifecycle,1) ./ vMass;
vDebtAggregate          	= vDebt(1:T_lifecycle,1) ./ vMass;
vLeverageAggregate      	= vDebtAggregate(1:T_lifecycle,1) ./ vCapitalAggregate;
vLaborAggregate         	= vLabor(1:T_lifecycle,1) ./ vMass;
vOutputAggregate        	= vOutput(1:T_lifecycle,1) ./ vMass;
vAvgProductivityAggregate	= vAvgProductivity(1:T_lifecycle,1) ./ vMass;
vAvgCashAggregate			= vAvgCash(1:T_lifecycle,1) ./ vMass;

vFracPosDebtAggregate       = vFracPosDebt(1:T_lifecycle,1) ./ vMass;
vAvgGrossLevAggregate       = vAvgGrossLev(1:T_lifecycle,1) ./ vMass;
vAvgNetLevAggregate         = vAvgNetLev(1:T_lifecycle,1) ./ vMass;

%%%
% Compute calibration targets
%%%

% For computing share of firms by age
vExtensiveMargin		= massEntrants * vMass;		% name comes from Clementi, Palazzo, and Wu (2017)

% Time aggregate to annual frequency
T_lifecycle_annual		= T_lifecycle / 4;
vExtensiveMarginAnnual	= zeros(T_lifecycle_annual,1);
for t = 1 : T_lifecycle_annual
	vExtensiveMarginAnnual(t,1)		= sum(vExtensiveMargin(4*(t-1)+1:4*t,1));
end

% Employment share by age
vEmploymentShare       = (vLaborAggregate .* vMass * massEntrants) / aggregateLabor;
vEmploymentShareAnnual = zeros(T_lifecycle_annual,1);
for t = 1 : T_lifecycle_annual
    vEmploymentShareAnnual(t,1) = sum(vLaborAggregate(4*(t-1)+1:4*t,1) .* vMass(4*(t-1)+1:4*t,1) * massEntrants) / aggregateLabor;
end


%%%
% Compute statistics of Compustat sample
%%%

cutoff          = 28;
vAge            = linspace(1,T_lifecycle,T_lifecycle)';

fracPosDebtCompustat        = sum(vFracPosDebt(cutoff:end)) / sum(vMass(cutoff:end));
avgGrossLeverageCompustat   = sum(vAvgGrossLev(cutoff:end)) / sum(vMass(cutoff:end));
avgNetLeverageCompustat     = sum(vAvgNetLev(cutoff:end)) / sum(vMass(cutoff:end));
employmentRatioCompustat    = (sum(vLabor(cutoff:end)) / sum(vMass(cutoff:end))) / (sum(vLabor(1:cutoff-1)) / sum(vMass(1:cutoff-1)));
ageRatioCompustat           = (sum(vAge(cutoff:end) .* vMass(cutoff:end)) / sum(vMass(cutoff:end))) / (sum(vAge(1:cutoff-1) .* vMass(1:cutoff-1)) / sum(vMass(1:cutoff-1)));
