%----------------------------------------------------------------
% Initializations
%----------------------------------------------------------------

%  Initialize distribution of incumbents at steady state
mDistribution_previous	     = mDistributionSS;
mDistributionSeries          = zeros(nStateDist,T);

% Preallocate some matrices
vAggregateInvestmentDemand	 = zeros(T,1);
vAggregateInvestmentSupply	 = zeros(T,1);
vAggregateConsumptionDemand	 = zeros(T,1);
vAggregateConsumptionSupply	 = zeros(T,1);
vNewCapitalProduction        = zeros(T,1);

vAggregateCapital            = aggregateCapitalSS * ones(T+1,1);
vAggregateOutput             = aggregateOutputSS * ones(T+1,1);
vAggregateLabor              = aggregateLaborSS * ones(T+1,1);
vAggregateTFP                = aggregateTFPSS * ones(T+1,1);
vAggregateMass               = ones(T+1,1);


%----------------------------------------------------------------
% Pre-compute heterogeneity groupings and steady state benchmarks
%----------------------------------------------------------------

% Evaluate steady state policy functions on the distribution grid
vCapitalPrimeDistSS   = interpn(mProdGrid,mCashGrid,mCapitalPrimeSS,...
                        mStateGridDist(:,1),mStateGridDist(:,2));
vDebtPrimeDistSS      = interpn(mProdGrid,mCashGrid,mDebtPrimeSS,...
                        mStateGridDist(:,1),mStateGridDist(:,2));

% Steady state leverage and distance to default
vLeverageSS           = vDebtPrimeDistSS ./ max(k0,vCapitalPrimeDistSS);
mDefaultCutoffDistSS  = repmat(vDefaultCutoffSS,[1 nCashDist]);
vDistanceDefaultSS    = mStateGridDist(:,2) - mDefaultCutoffDistSS(:);

% Weights used to compute steady state moments (continuing incumbents)
vWeightsSS            = vDistContContinueSS .* mDistributionSS(:);
if sum(vWeightsSS) == 0
    vWeightsSS        = ones(size(vWeightsSS));
end
vWeightsSS            = vWeightsSS ./ sum(vWeightsSS);

% Helper to compute weighted medians
[~,idxSortLev]        = sort(vLeverageSS);
cumWeightsLev         = cumsum(vWeightsSS(idxSortLev));
idxMedianLev          = find(cumWeightsLev >= 0.5,1,'first');
if isempty(idxMedianLev)
    idxMedianLev      = length(idxSortLev);
end
levThreshold          = vLeverageSS(idxSortLev(idxMedianLev));

[~,idxSortDist]       = sort(vDistanceDefaultSS);
cumWeightsDist        = cumsum(vWeightsSS(idxSortDist));
idxMedianDist         = find(cumWeightsDist >= 0.5,1,'first');
if isempty(idxMedianDist)
    idxMedianDist     = length(idxSortDist);
end
distThreshold         = vDistanceDefaultSS(idxSortDist(idxMedianDist));

% Indicator vectors for each group
vGroupLeverageLow     = (vLeverageSS <= levThreshold);
vGroupLeverageHigh    = (vLeverageSS >  levThreshold);
if sum(vWeightsSS .* vGroupLeverageHigh) == 0
    vGroupLeverageHigh= (vLeverageSS >= levThreshold);
end
if sum(vWeightsSS .* vGroupLeverageLow) == 0
    vGroupLeverageLow = (vLeverageSS < levThreshold);
end
if sum(vWeightsSS .* vGroupLeverageHigh) == 0
    vGroupLeverageHigh(idxSortLev(end)) = true;
end
if sum(vWeightsSS .* vGroupLeverageLow) == 0
    vGroupLeverageLow(idxSortLev(1)) = true;
end

vGroupDistanceClose   = (vDistanceDefaultSS <= distThreshold);
vGroupDistanceFar     = (vDistanceDefaultSS >  distThreshold);
if sum(vWeightsSS .* vGroupDistanceFar) == 0
    vGroupDistanceFar = (vDistanceDefaultSS >= distThreshold);
end
if sum(vWeightsSS .* vGroupDistanceClose) == 0
    vGroupDistanceClose = (vDistanceDefaultSS < distThreshold);
end
if sum(vWeightsSS .* vGroupDistanceFar) == 0
    vGroupDistanceFar(idxSortDist(end)) = true;
end
if sum(vWeightsSS .* vGroupDistanceClose) == 0
    vGroupDistanceClose(idxSortDist(1)) = true;
end

vGroupLeverageLow     = double(vGroupLeverageLow);
vGroupLeverageHigh    = double(vGroupLeverageHigh);
vGroupDistanceClose   = double(vGroupDistanceClose);
vGroupDistanceFar     = double(vGroupDistanceFar);

% Steady state averages for each group
massLeverageLowSS     = sum(vWeightsSS .* vGroupLeverageLow);
massLeverageHighSS    = sum(vWeightsSS .* vGroupLeverageHigh);
massDistanceCloseSS   = sum(vWeightsSS .* vGroupDistanceClose);
massDistanceFarSS     = sum(vWeightsSS .* vGroupDistanceFar);

vCapitalByLeverageSS              = zeros(1,2);
vCapitalByDefaultDistanceSS       = zeros(1,2);
if massLeverageLowSS > 0
    vCapitalByLeverageSS(1,1)     = sum(vCapitalPrimeDistSS .* vWeightsSS .* vGroupLeverageLow) / massLeverageLowSS;
else
    vCapitalByLeverageSS(1,1)     = NaN;
end
if massLeverageHighSS > 0
    vCapitalByLeverageSS(1,2)     = sum(vCapitalPrimeDistSS .* vWeightsSS .* vGroupLeverageHigh) / massLeverageHighSS;
else
    vCapitalByLeverageSS(1,2)     = NaN;
end
if massDistanceCloseSS > 0
    vCapitalByDefaultDistanceSS(1,1) = sum(vCapitalPrimeDistSS .* vWeightsSS .* vGroupDistanceClose) / massDistanceCloseSS;
else
    vCapitalByDefaultDistanceSS(1,1) = NaN;
end
if massDistanceFarSS > 0
    vCapitalByDefaultDistanceSS(1,2) = sum(vCapitalPrimeDistSS .* vWeightsSS .* vGroupDistanceFar) / massDistanceFarSS;
else
    vCapitalByDefaultDistanceSS(1,2) = NaN;
end

% Preallocate time series for heterogeneity analysis
mCapitalByLeverage                = NaN(T,2);
mMassByLeverage                   = NaN(T,2);
mCapitalByDefaultDistance         = NaN(T,2);
mMassByDefaultDistance            = NaN(T,2);


%----------------------------------------------------------------
% Compute the path of objects
%----------------------------------------------------------------

for t = 1 : T


    %----------------------------------------------------------------
    % Compute objects that require knowing last period's distribution
    %---------------------------------------------------------------

    %%%
    % Interpolate policy functions
    %%%

    % Decide which policy functions to use
    if t == 1           % use steady state decision rules

    	mCapitalPrime       = reshape(mCapitalPrimeSeries(:,T+1),nProd,nCash);
    	mDebtPrime          = reshape(mDebtPrimeSeries(:,T+1),nProd,nCash);
        vDistContContinue   = vDistContContinueSS;

    elseif t > 1        % use previous period's decision rules

        mCapitalPrime     = reshape(mCapitalPrimeSeries(:,t-1),nProd,nCash);
        mDebtPrime        = reshape(mDebtPrimeSeries(:,t-1),nProd,nCash);
        vDistContContinue = (mStateGridDist(:,2) >= reshape(repmat(mDefaultCutoffSeries(:,t-1),...
                            [1 nCashDist]),nStateDist,1));

    end

    vContinue               = ppiExit * (mStateGridDist(:,2) >= 0) + (1 - ppiExit) * vDistContContinue;  % NB: why do we need this?

    % Capital accumulation policy
    vCapitalPrimeDist       = interpn(mProdGrid,mCashGrid,mCapitalPrime,...
                                mStateGridDist(:,1),mStateGridDist(:,2));
    mCapitalPrimeDist       = reshape(vCapitalPrimeDist,nProd,nCashDist);

    % Debt accumulation policy
    vDebtPrimeDist          = interpn(mProdGrid,mCashGrid,mDebtPrime,...
                                mStateGridDist(:,1),mStateGridDist(:,2));
    mDebtPrimeDist          = reshape(vDebtPrimeDist,nProd,nCashDist);


    %%%
    % Compute evolution of cash on hand for incumbents
    %%%

    % Evolution of state variables; CONVENTION: rows = this period's shocks; columns = last period's state
    mCapitalPrimePrimeDist		 = repmat(vCapitalPrimeDist',[nShocks 1]);
    mDebtPrimePrimeDist			 = repmat(vDebtPrimeDist',[nShocks 1]);
    mProdPrimeDist				 = repmat(mShocksGrid(:,1),[1 nStateDist]);
    mOmegaPrimeDist				 = repmat(mShocksGrid(:,2),[1 nStateDist]);

    % This period's cash on hand
    mProfitPrimeDist		 = vA(t,1) * (exp(mProdPrimeDist) .^ (1 / (1 - nnu))) .* ((exp(mOmegaPrimeDist) .* mCapitalPrimePrimeDist) .^ tthetaHat) .* ...
    							(vWage(t,1) ^ (-nnu / (1 - nnu)));
    mCashPrimeDist			 = mProfitPrimeDist + vQ(t,1) * (1 - ddelta) * exp(mOmegaPrimeDist) .* mCapitalPrimePrimeDist - ...
                                (mDebtPrimePrimeDist / vAggregateInflation(t,1)) - ppsi_0;
    vCashPrimeDist			 = mCashPrimeDist(:);

    % This period's default decision (to condition on survival)
    mContinueExitDist		 = (mCashPrimeDist >= 0);
    mContinueContDist		 = (mCashPrimeDist >= reshape(repmat(mDefaultCutoffSeries(:,t),[1 nOmega * nStateDist]),nShocks,nStateDist));
    mContinueIncumbent       = ppiExit * mContinueExitDist + (1 - ppiExit) * mContinueContDist;


    %%%
    % Cash on hand for new entrants
    %%%

    % Compute implied cash on hand for new entrants
    vProfitImpliedEnt   = vA(t,1) * (exp(mShocksGrid(:,1)) .^ (1 / (1 - nnu))) .* ((exp(mShocksGrid(:,2)) * k0) .^ tthetaHat) .* ...
                            (vWage(t,1) ^ (-nnu / (1 - nnu)));
    vCashImpliedEnt     = vProfitImpliedEnt + vQ(t,1) * (1 - ddelta) * exp(mShocksGrid(:,2)) * k0 - b0 - ppsi_0;

    % Default decision for new entrants
    vContinueExitEnt    = (vCashImpliedEnt >= 0);
    vContinueContEnt    = (vCashImpliedEnt >= reshape(repmat(mDefaultCutoffSeries(:,t),[1 nOmega]),nShocks,1));
    vContinueEntrant	= ppiExit * vContinueExitEnt + (1 - ppiExit) * vContinueContEnt;


    %%%
    % Compute transition matrices and mass of new entrants from t-1 to t
    %%%

    transition_path_update_distribution;


    %%%
    % Employment
    %%%

    % Incumbents
    aLaborIntegrandIncumbent	= reshape((((vP(t,1) * nnu * exp(mProdPrimeDist) .* ((exp(mOmegaPrimeDist) .* mCapitalPrimePrimeDist) .^ ttheta)) / ...
                                    vWage(t,1)) .^ (1 / (1 - nnu))) .* mContinueIncumbent,nShocks,nProd,nCashDist);
    mLaborIntegrandIncumbent  = zeros(nProd,nCashDist);
    for iProd = 1 : nProd
        mLaborIntegrandIncumbent(iProd,:)   = mShocksTransition(iProd,:) * squeeze(aLaborIntegrandIncumbent(:,iProd,:));
    end

    % New entrants
    vLaborIntegrandEntrant		 = (((vP(t,1) * nnu * exp(mShocksGrid(:,1)) .* ((exp(mShocksGrid(  :,2)) * k0) .^ ttheta)) / ...
                                    vWage(t,1)) .^ (1 / (1 - nnu))) .* vContinueEntrant;

    % Aggregate
    vAggregateLabor(t,1)         = (1 - ppiExit) * sum(mLaborIntegrandIncumbent(:) .* mDistribution_previous(:) .* ...
                                    vDistContContinue) + massEntrants * sum(vLaborIntegrandEntrant .* vDistEnt);


    %%%
    % Output
    %%%

    % Incumbents
    aOutputIntegrandIncumbent	= reshape(exp(mProdPrimeDist) .* ((exp(mOmegaPrimeDist) .* mCapitalPrimePrimeDist) .^ ttheta) .* ...
                                    (reshape(aLaborIntegrandIncumbent,nShocks,nStateDist) .^ nnu) .* ...
                                    mContinueIncumbent,nShocks,nProd,nCashDist);
    mOutputIntegrandIncumbent   = zeros(nProd,nCashDist);
    for iProd = 1 : nProd
        mOutputIntegrandIncumbent(iProd,:)   = mShocksTransition(iProd,:) * squeeze(aOutputIntegrandIncumbent(:,iProd,:));
    end

    % New entrants
    vOutputIntegrandEntrant		= exp(mShocksGrid(:,1)) .* ((exp(mShocksGrid(:,2)) * k0) .^ ttheta) .* (vLaborIntegrandEntrant .^ nnu) .* vContinueEntrant;


    % Aggregate
    vAggregateOutput(t,1)       = (1 - ppiExit) * sum(mOutputIntegrandIncumbent(:) .* mDistribution_previous(:) .* vDistContContinue) + ...
    							     massEntrants * sum(vOutputIntegrandEntrant .* vDistEnt);


    %----------------------------------------------------------------
    % Compute objects involving current distribution
    %----------------------------------------------------------------

    %%%
    % Compute current distribution
    %%%

    mDistribution           = reshape((1 - ppiExit) * mTransitionIncumbents * (vDistContContinue .* mDistribution_previous(:)) + ...
    							massEntrants * vDistributionEntrants,nProd,nCashDist);
    mDistribution_previous  = mDistribution;
    mDistributionSeries(:,t)= mDistribution(:);

    %%%
    % Interpolate policy functions from this period
    %%%

    mCapitalPrime       = reshape(mCapitalPrimeSeries(:,t),nProd,nCash);
    mDebtPrime          = reshape(mDebtPrimeSeries(:,t),nProd,nCash);

    vDistContContinue   = (mStateGridDist(:,2) >= reshape(repmat(mDefaultCutoffSeries(:,t),...
                          [1 nCashDist]),nStateDist,1));
    vContinue           = ppiExit * (mStateGridDist(:,2) >= 0) + (1 - ppiExit) * vDistContContinue;

    % Capital accumulation policy
    vCapitalPrimeDist   = interpn(mProdGrid,mCashGrid,mCapitalPrime,...
                          mStateGridDist(:,1),mStateGridDist(:,2));
    mCapitalPrimeDist   = reshape(vCapitalPrimeDist,nProd,nCashDist);

    % Debt accumulation policy
    vDebtPrimeDist      = interpn(mProdGrid,mCashGrid,mDebtPrime,...
                          mStateGridDist(:,1),mStateGridDist(:,2));
    mDebtPrimeDist      = reshape(vDebtPrimeDist,nProd,nCashDist);


    %%%
    % Heterogeneity: average capital by leverage and default distance
    %%%

    vWeightsCurrent     = mDistribution(:) .* vDistContContinue;

    massLowLev          = sum(vWeightsCurrent .* vGroupLeverageLow);
    massHighLev         = sum(vWeightsCurrent .* vGroupLeverageHigh);
    massCloseDefault    = sum(vWeightsCurrent .* vGroupDistanceClose);
    massFarDefault      = sum(vWeightsCurrent .* vGroupDistanceFar);

    mMassByLeverage(t,1)            = massLowLev;
    mMassByLeverage(t,2)            = massHighLev;
    mMassByDefaultDistance(t,1)     = massCloseDefault;
    mMassByDefaultDistance(t,2)     = massFarDefault;

    if massLowLev > 0
        mCapitalByLeverage(t,1)     = sum(vCapitalPrimeDist .* vWeightsCurrent .* vGroupLeverageLow) / massLowLev;
    end
    if massHighLev > 0
        mCapitalByLeverage(t,2)     = sum(vCapitalPrimeDist .* vWeightsCurrent .* vGroupLeverageHigh) / massHighLev;
    end
    if massCloseDefault > 0
        mCapitalByDefaultDistance(t,1) = sum(vCapitalPrimeDist .* vWeightsCurrent .* vGroupDistanceClose) / massCloseDefault;
    end
    if massFarDefault > 0
        mCapitalByDefaultDistance(t,2) = sum(vCapitalPrimeDist .* vWeightsCurrent .* vGroupDistanceFar) / massFarDefault;
    end

   
    %%%
    % Aggregate Capital
    %%%

    vAggregateCapital(t+1,1)  = (1 - ppiExit) * sum(vCapitalPrimeDist .* vDistContContinue .* mDistribution(:)) + ...
    							  massEntrants * sum(vCapitalPrimeDist .* vDistContContinue .* vDistributionEntrants);


end


%%%
% Investment market clearing
%%%

% Demand
vNewCapitalProduction(1:T,1)        = vAggregateCapital(2:T+1,1) - (1 - ddelta) * EOmegaTerm2 * vAggregateCapital(1:T,1) + ...
                                        (1 - (1 - ddelta) * EOmegaTerm2) * k0 * massEntrants;
vAggregateInvestmentDemand(1:T,1)   = ((((vNewCapitalProduction(1:T) ./ vAggregateCapital(1:T,1)) + (ddeltaHat / (pphiCapital - 1))) .* ...
                                        ((1 - (1 / pphiCapital)) / (ddeltaHat ^ (1 / pphiCapital)))) .^ (pphiCapital / (pphiCapital - 1))) .* vAggregateCapital(1:T,1);

% Supply
vAggregateInvestmentSupply(1:T,1)	= (vAggregatePriceCapital(1:T,1) .^ pphiCapital) .* ddeltaHat .* vAggregateCapital(1:T,1);


%%%
% Consumption market clearing
%%%

% Demand
vAggregateConsumptionDemand(1:T,1)	= vAggregateMarginalUtility(1:T,1) .^ (-1 / ssigma);

% Supply
vAggregateConsumptionSupply(1:T,1)	= vAggregateOutput(1:T,1) - vAggregateInvestmentDemand(1:T,1) - ppsi_0;


%%%
% Other aggregates
%%%

vAggregateConsumption	= vAggregateConsumptionDemand;
vAggregateInvestment	= vAggregateInvestmentSupply;
vAggregateR_real		= vAggregateR_nom(1:T,1) ./ vAggregateInflation(2:T+1,1);
