%%%
% Initialize aggregate series
%%%

% Steady state output and marginal utility
aggregateOutputSteadyState		= aggregateOutputSS;
aggregateMarginalSteadyState	= aggregateMarginalUtilitySS;
wageSteadyState					= wageSS;

% Initialize aggregate series
vAggregateOutput                                 = aggregateOutputSS * ones(T+1,1);
vAggregateConsumption 			                 = aggregateConsumptionSS * ones(T+1,1);
vAggregateInvestment 	                         = aggregateInvestmentSS * ones(T+1,1);
vAggregateLabor                                  = aggregateLaborSS * ones(T+1,1);
vAggregateCapital                                = aggregateCapitalSS * ones(T+1,1);
vAggregateTFP                                    = aggregateTFPSS * ones(T+1,1);
vAggregateDebt									 = aggregateDebtSS * ones(T+1,1);
vAggregateMass                                   = ones(T+1,1);
vAggregateMarginalUtility                        = (aggregateConsumptionSS ^ (-ssigma)) * ones(T+1,1);
vAggregateWage                                   = wageSS * ones(T+1,1);
%vAggregateCreditSpread					         %= %aggregateCreditSpreadSS * ones(T+1,1);%


% Serie AR(1) para el spread soberano
s_t = zeros(T+1,1);
s_t(1) = aggregateCreditSpreadSS;    % inicia en estado estacionario
for tt = 2 : T+1
    s_t(tt) = rrhoS * s_t(tt-1) + ssigmaS * randn;
end
vAggregateCreditSpread = s_t;



vAggregateInflation                              = ones(T+1,1);
vAggregatePriceIntermediate                      = pSS * ones(T+1,1);
vAggregatePriceCapital                           = ones(T+1,1);

vAggregateR_nom                                  = (1 / bbeta) * ones(T+1,1);
vAggregateR_real                                 = (1 / bbeta) * ones(T+1,1);

% Initialize individual decisions
mDefaultCutoffSeries                             = repmat(vDefaultCutoffSS,[1 T+1]);
mDebtPriceSeries                                 = repmat(mDebtPriceSS(:),[1 T+1]);

mUnconstrainedCutoffSeries                       = repmat(vUnconstrainedCutoffSS,[1 T+1]);
mCapitalUnconstrainedSeries                      = repmat(vCapitalUnconstrainedSS,[1 T+1]);
mDebtUnconstrainedSeries                         = repmat(vDebtUnconstrainedSS,[1 T+1]);
mValueUnconstrainedSeries                        = repmat(vValueUnconstrainedSS,[1 T+1]);

mValueSeries                                     = repmat(mValueSS(:),[1 T+1]);
mCapitalPrimeSeries                              = repmat(mCapitalPrimeSS(:),[1 T+1]);
mDebtPrimeSeries                                 = repmat(mDebtPrimeSS(:),[1 T+1]);
mDividendsSeries                                 = repmat(mDividendsSS(:),[1 T+1]);
mDebtPriceOptimalSeries                          = repmat(mDebtPriceOptimalSS(:),[1 T+1]);

% Exit flags for the continuous optimizer
aExitFlagSeries							         = zeros(nProd,nCash,T);

%%%
% Compute transition path
%%%

% Get prices and reset some prices to steady state
load prices_het_agent.mat;
vAggregatePriceCapital(1:T,1)        	= exp(vPrices(1:T,1));
vAggregateMarginalUtility(1:T,1)     	= exp(vPrices(T+1:2*T,1)) * aggregateMarginalUtilitySS;
transition_path_price_series;

% RESET SOME PRICES TO STEADY STATE
vQ(1:end,1)								= qSS;
vAggregateR_real(1:end,1)			    = 1/bbeta;
vSDF(1:end,1) 							= bbeta;

%%%
% Compute decisions and aggregates
%%%

% Decision rules
transition_path_decisions;

% Aggregate series
transition_path_update_aggregates;


%%%
% Compute elasticities across the state space
%%%

% Interpolate capital policy functions over distribution grid
vCapitalPrimeDist       = interpn(mProdGrid,mCashGrid,reshape(mCapitalPrimeSeries(:,1),nProd,nCash),...
                            mStateGridDist(:,1),mStateGridDist(:,2));
vCapitalPrimeDistSS     = interpn(mProdGrid,mCashGrid,reshape(mCapitalPrimeSS,nProd,nCash),...
                            mStateGridDist(:,1),mStateGridDist(:,2));

% Interpolate capital policy functions over distribution grid
vDebtPrimeDist          = interpn(mProdGrid,mCashGrid,reshape(mDebtPrimeSeries(:,1),nProd,nCash),...
                            mStateGridDist(:,1),mStateGridDist(:,2));
vDebtPrimeDistSS        = interpn(mProdGrid,mCashGrid,reshape(mDebtPrimeSS,nProd,nCash),...
                            mStateGridDist(:,1),mStateGridDist(:,2));

% Compute elasticities
mCapitalElasticity      = reshape((log(vCapitalPrimeDist) - log(vCapitalPrimeDistSS)) / (-4*vEpsilon_m(1,1)),nProd,nCashDist);
mDebtElasticity         = reshape((vDebtPrimeDist - vDebtPrimeDistSS) / (-4 * vEpsilon_m(1,1)),nProd,nCashDist);

% Distribution
mDistribution           = reshape(mDistributionSeries(:,1),nProd,nCashDist);

% Compute average elasticities w.r.t. productivity
vCapitalElasticityAverage   = sum(repmat(vProdErgodic,[1 nCashDist]) .* mCapitalElasticity,1)';
vDebtElasticityAverage      = sum(repmat(vProdErgodic,[1 nCashDist]) .* mDebtElasticity,1)';
vDistributionAverage        = sum(mDistribution,1)';
