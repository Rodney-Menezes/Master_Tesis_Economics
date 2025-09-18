%----------------------------------------------------------------
% Compute error in the iteration
%----------------------------------------------------------------

err_consumption		= max(abs(vAggregateConsumptionDemand(1:T,1) - vAggregateConsumptionSupply(1:T,1)));
err_investment		= max(abs(vAggregateInvestmentDemand(1:T,1) - vAggregateInvestmentSupply(1:T,1)));
errMIT				= max([err_consumption err_investment]);


%----------------------------------------------------------------
% Update series
%----------------------------------------------------------------

% Preallocate to be the myAD type (don't think this is necessary)
vAggregateMarginalUtilityNew		= vAggregateMarginalUtility;
vAggregatePriceCapitalNew			= vAggregatePriceCapital;
vAggregateOutputNew                 = vAggregateOutput;

% Compute new series
vAggregateMarginalUtilityNew(1:T,1)	= vAggregateConsumptionSupply(1:T,1) .^ (-1 / ssigma);
vAggregateOutputNew(1:T,1)			= vAggregateOutputSupply(1:T,1);
vAggregatePriceCapitalNew(1:T,1)	= (vAggregateInvestmentDemand(1:T,1) ./ (ddelta * vAggregateCapital(1:T,1))) .^ (1 / pphiCapital);

% Update
vAggregatePriceCapital(1:T,1) 	= dampeningQ * vAggregatePriceCapital(1:T,1) + (1 - dampeningQ) * vAggregatePriceCapitalNew(1:T,1);
vAggregateMarginalUtility(1:T,1)= dampeningMu * vAggregateMarginalUtility(1:T,1) + (1 - dampeningMu) * vAggregateMarginalUtilityNew(1:T,1);
vAggregateOutput(1:T,1)			= dampeningY * vAggregateOutput(1:T,1) + (1 - dampeningY) * vAggregateOutputNew(1:T,1);

% Update counter
iterationMIT			= iterationMIT + 1;

% Update other series
vAggregateConsumption	= vAggregateConsumptionDemand;
vAggregateInvestment	= vAggregateInvestmentSupply;
vAggregateR_real		= vAggregateR_nom(1:T,1) ./ vAggregateInflation(2:T+1,1);
