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
%vAggregateCreditSpread					         = aggregateCreditSpreadSS * ones(T+1,1);



vAggregateInflation                              = ones(T+1,1);
vAggregatePriceIntermediate                      = pSS * ones(T+1,1);
vAggregatePriceCapital                           = ones(T+1,1);



% (Inicialización “semilla” para Taylor extendido)
vAggregateR_nom   = NaN(T+1,1);
vAggregateR_nom(1)= 1/bbeta;    % punto de partida igual al steady state

% vAggregateR_real lo calcularemos tras la regla extendida





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



%%% Regla de Taylor extendida con spread
for t = 2 : T+1
    vAggregateR_nom(t) = exp( ...
        rho * log(vAggregateR_nom(t-1)) + ...
        (1-rho) * ( pphiInflation * log(vAggregateInflation(t)) ...
                  + pphiOutput    * log(vAggregateOutput(t) / aggregateOutputSteadyState) ) ...
        + pphiS * vAggregateCreditSpread(t) ...
        + vEpsilon_m(t) ...
    );
end
% Ahora calculamos la tasa real
vAggregateR_real = vAggregateR_nom ./ vAggregateInflation;





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
