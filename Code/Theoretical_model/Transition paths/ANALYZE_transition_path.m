% Computes and analyzes response to MIT monetary shock
%
% Pablo Ottonello and Thomas Winberry
% This draft: June 27th, 2020

% Housekeeping
clear all
close all

%----------------------------------------------------------------
% Set options for what is to be done
%----------------------------------------------------------------

global option_continuous_optimization option_continuous_optimization_VFI option_updating_rule

% How to solve for decision rules
option_continuous_optimization			= 1;			% use continuous optimization after discretized VFI has converged
option_continuous_optimization_VFI		= 1;			% use continuous optimization in value function iteration

% Options for the MIT shock iteration
option_updating_rule					= 1;			% = 0 for naive updating
														% = 1 for quasi-Newton updating (using derivatives of rep agent model) [dampening = .75 works and dampening = .5 seems to work for T = 200]

option_continuous_optimization_VFI_init	= option_continuous_optimization_VFI;	% use continuous optimization at each step of the transition path


%----------------------------------------------------------------
% Set parameters
%----------------------------------------------------------------

% Parameters of economic model
set_parameters_model;

% Parameters of numerical approximations
set_parameters_numerical;

% Path of shocks
global vEpsilon_m
vEpsilon_m			        = zeros(T+1,1);
vEpsilon_m(1,1)	            = -1*ssigmaM;
for t = 1 : ceil(T/5)
	vEpsilon_m(t+1,1) = rrhoM * vEpsilon_m(t,1);
end


%----------------------------------------------------------------
% Compute steady state
%----------------------------------------------------------------

option_continuous_optimization_VFI		= 0;

%%%
% Solve for the market clearing wage using bisection
%%%

% Preliminaries
resid                           = 100;
iteration_market_clearing       = 1;
tolerance_market_clearing       = 5e-3;
max_iterations_market_clearing  = 20;

while abs(resid) > tolerance_market_clearing && iteration_market_clearing <= max_iterations_market_clearing

	% Compute residual
	wageSS  = .5 * (wageMin + wageMax);
	resid   = market_clearing_residual(wageSS);
	if isnan(resid)
		resid = -1e4;       % wage so high that no firms enter
	end

	% Update bounds on wage
	if resid > 0        	% labor demand too high ==> increase wage
		wageMin = wageSS;
	else                	% labor demand too low ==> decrease wage
		wageMax = wageSS;
	end

	% Update iteration
	iteration_market_clearing       = iteration_market_clearing + 1;

end

% Set wage
wage		= wageSS;

% Compute the steady state using continuous optimization
option_continuous_optimization_VFI 			= option_continuous_optimization_VFI_init;

% Compute objects
core_steady_state;
steady_state_aggregates;

% Calibrate disutility of labor supply
cchi                		= wageSS * (aggregateConsumption ^ (-ssigma));

% Calibrate effective depreciation
global ddeltaHat
ddeltaHat					= (1 - (1 - ddelta) * EOmegaTerm2) * (1 + (massEntrants * k0 / aggregateCapital));


%----------------------------------------------------------------
% Some computations involving the steady state
%----------------------------------------------------------------

%%%
% Rename steady state objects for analysis below
%%%

% Aggregates
aggregateCapitalSS			                         = aggregateCapital;
aggregateLaborSS				                     = aggregateLabor;
aggregateOutputSS				                     = aggregateOutput;
aggregateInvestmentSS		                         = aggregateInvestment;
aggregateConsumptionSS	                             = aggregateConsumption;
aggregateMarginalUtilitySS	                         = aggregateConsumption ^ (-ssigma);
aggregateDebtSS					                     = aggregateDebt;
aggregateTFPSS					                     = aggregateTFP;
wageSS									             = wage;
aggregateCreditSpreadSS								 = meanSpread;

% Default decisions
vDefaultCutoffSS			                         = vDefaultCutoff;
mDebtPriceSS					                     = mDebtPrice;

% Unconstrained decisions
vUnconstrainedCutoffSS		                         = vUnconstrainedCutoff;
vCapitalUnconstrainedSS		                         = vCapitalUnconstrained;
vDebtUnconstrainedSS			                     = vDebtUnconstrained;
vValueUnconstrainedSS			                     = vValueUnconstrained;

% Constrained decisions
mValueSS						                     = mValue;
mCapitalPrimeSS			                             = mCapitalPrime;
mDebtPrimeSS				                         = mDebtPrime;
mDividendsSS				                         = mDividends;
mDebtPriceOptimalSS	                                 = reshape(vDebtPriceOptimal,nProd,nCash);

% Distribution
mDistributionSS                                      = mDistribution;
vDistContContinueSS                                  = vDistContContinue;


%%%
% Recompute representative agent steady state given new value of cchi
%%%

% Solve for labor supply
global nRepSS
nRepSS		         = nSS;

% Depreciation rate
ddeltaHet 			 = ddelta;
ddelta 				 = ddeltaHat;

% Compute other aggregates in steady state
pSS					 = (ggamma - 1) / ggamma; 				% relative price of production firms' output
kRepSS			     = ((pSS * ttheta * (nRepSS ^ nnu)) / (qSS * ((1 / bbeta) - (1 - ddelta)))) ^ ...
                    	(1 / (1 - ttheta)); 					% capital stock
yRepSS			     = (kRepSS ^ ttheta) * (nRepSS ^ nnu); 	% output

% Additional quantities from representative agent steady state
cRepSS               = yRepSS - ddelta * kRepSS - ppsi_0;	% consumption
muRepSS			     = cRepSS ^ (-ssigma);					% marginal utility
wRepSS			     = nnu * pSS * yRepSS / nRepSS; 		% real wage
inflationSS		     = 1; 									% steady state inflation

% Scale of the economy
cchiRep		     	 = wRepSS * muRepSS;
cchiHet		     	 = cchi;


%----------------------------------------------------------------
% Compute the heterogeneous agent transition path
%----------------------------------------------------------------

%%%
% Set some parameters and initialize relevant objects
%%%

% Select appropriate scale of the economy
cchi		= cchiHet;
ddelta 		= ddeltaHet;

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
vAggregateDebt									 = aggregateDebtSS * ones(T+1,1);
vAggregateMarginalUtility                        = (aggregateConsumptionSS ^ (-ssigma)) * ones(T+1,1);
vAggregateWage                                   = wageSS * ones(T+1,1);
vAggregateCreditSpread					         = aggregateCreditSpreadSS * ones(T+1,1);

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

% Compute initial price vector
load prices_het_agent.mat;

%%%
% Compute aggregate prices from price path
%%%

vAggregatePriceCapital(1:T,1)        = exp(vPrices(1:T,1));
vAggregateMarginalUtility(1:T,1)     = exp(vPrices(T+1:2*T,1)) * aggregateMarginalUtilitySS;

transition_path_price_series;


%%%
% Compute aggregate series along transition path
%%%

% Decision rules
transition_path_decisions;

% Aggregate series
transition_path_update_aggregates;

%----------------------------------------------------------------
% Impulse response plots from paper
%----------------------------------------------------------------

%%%
% Aggregate transmission mechanism
%%%

figure

h               = gcf;
h.PaperUnits    = 'inches';
h.PaperPosition = [0 0 15 4];

vTime = linspace(1,T,T)';

subplot(1,3,1)
hold on
plot(vTime,400 * (vAggregateR_real(1:T,1) - (1 / bbeta)),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,400 * log(vAggregateInflation(1:T,1)),'linewidth',1.5,'linestyle','-.','color',[.5,.24,.5])
plot(vTime,400 * (vAggregateR_nom(1:T,1) - (1 / bbeta)),'linewidth',1.5,'linestyle','--','color',[178/255,34/255,34/255])
plot(vTime,zeros(T,1),'linewidth',1.5,'linestyle','--','color','k')
xlim([1 12])
h	 = legend('Real rate','Inflation','Nominal rate');
set(h,'interpreter','latex','location','northeast','fontsize',14)
set(gcf,'color','w')
xlabel('Quarters','interpreter','latex')
ylabel('Annualized p.p. deviation','interpreter','latex')
grid on
title('Interest Rates and Inflation','interpreter','latex','fontsize',14)
hold off

subplot(1,3,2)
hold on
plot(vTime,100 * log(vAggregateOutput(1:T,1) / aggregateOutputSS),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,100 * log(vAggregateConsumption(1:T,1) / aggregateConsumptionSS),'linewidth',1.5,'linestyle','-.','color',[.5,.24,.5])
plot(vTime,100 * log(vAggregateInvestment(1:T,1) / (aggregateInvestmentSS)),'linewidth',1.5,'linestyle','--','color',[178/255,34/255,34/255])
plot(vTime,zeros(T,1),'linewidth',1.5,'linestyle','--','color','k')
xlim([1 12])
h	 = legend('Output','Consumption','Investment');
set(h,'interpreter','latex','location','northeast','fontsize',14)
set(gcf,'color','w')
xlabel('Quarters','interpreter','latex')
ylabel('$\%$ deviation','interpreter','latex')
grid on
title('Aggregate Quantities','interpreter','latex','fontsize',14)
hold off

subplot(1,3,3)
hold on
plot(vTime,100 * log(vAggregatePriceIntermediate(1:T,1) / pSS),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,100 * log(vAggregatePriceCapital(1:T,1)),'linewidth',1.5,'linestyle','-.','color',[.5,.24,.5])
plot(vTime,100 * log(vAggregateWage(1:T,1) / (wageSS)),'linewidth',1.5,'linestyle','--','color',[178/255,34/255,34/255])
plot(vTime,zeros(T,1),'linewidth',1.5,'linestyle','--','color','k')
xlim([1 12])
h	 = legend('Intermediate Good Price','Capital Price','Real Wage');
set(h,'interpreter','latex','location','northeast','fontsize',14)
set(gcf,'color','w')
xlabel('Quarters','interpreter','latex')
ylabel('$\%$ deviation','interpreter','latex')
grid on
title('Prices','interpreter','latex','fontsize',14)
hold off

print('../Results/agg_transmission_mechanism.eps','-depsc')


%%%
% Comparison to representative agent model
%%%

load rep_agent_path.mat

figure

h               = gcf;
h.PaperUnits    = 'inches';
h.PaperPosition = [0 0 13 4];

subplot(1,2,1)
hold on
plot(vTime,400 * (vAggregateR_real(1:T,1) - (1 / bbeta)),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,400 * (vAggregateR_real_rep_agent(1:T,1) - (1 / bbeta)),'linewidth',1.5,'linestyle','--','color',[178/255,34/255,34/255])
plot(vTime,zeros(T,1),'linewidth',1.5,'linestyle','--','color','k')
xlim([1 12])
xlabel('Quarters','interpreter','latex')
ylabel('Annualized p.p. deviation','interpreter','latex')
h	 = legend('Het agent','Rep agent');
set(h,'interpreter','latex','location','southwest','fontsize',14)
set(gcf,'color','w')
grid on
title('Real Interest Rate','interpreter','latex','fontsize',14)
hold off

subplot(1,2,2)
hold on
plot(vTime,100 * log(vAggregateInvestment(1:T,1) / (aggregateInvestmentSS)),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,100 * log(vAggregateInvestment_rep_agent(1:T,1) / (ddeltaHat * kRepSS)),'linewidth',1.5,'linestyle','--','color',[178/255,34/255,34/255])
plot(vTime,zeros(T,1),'linewidth',1.5,'linestyle','--','color','k')
xlim([1 12])
xlabel('Quarters','interpreter','latex')
ylabel('$\%$ deviation','interpreter','latex')
set(gcf,'color','w')
grid on
title('Aggregate Investment','interpreter','latex','fontsize',14)
hold off

print('../Results/rep_firm_comparison.eps','-depsc')


%%%
% Heterogeneity across leverage and default distance
%%%

capitalSeriesLeverage  = [vCapitalByLeverageSS; mCapitalByLeverage(1:T,:)];
capitalSeriesDefault   = [vCapitalByDefaultDistanceSS; mCapitalByDefaultDistance(1:T,:)];
massSeriesLeverage     = [massLeverageLowSS massLeverageHighSS; mMassByLeverage(1:T,:)];
massSeriesDefault      = [massDistanceCloseSS massDistanceFarSS; mMassByDefaultDistance(1:T,:)];

totalCapitalLeverage   = capitalSeriesLeverage .* massSeriesLeverage;
totalCapitalDefault    = capitalSeriesDefault .* massSeriesDefault;

totalCapitalLeverage(~isfinite(totalCapitalLeverage))   = 0;
totalCapitalDefault(~isfinite(totalCapitalDefault))     = 0;

currentCapitalLeverage = totalCapitalLeverage(1:T,:);
nextCapitalLeverage    = totalCapitalLeverage(2:end,:);
currentCapitalDefault  = totalCapitalDefault(1:T,:);
nextCapitalDefault     = totalCapitalDefault(2:end,:);

denomLeverage = max(currentCapitalLeverage,eps);
denomDefault  = max(currentCapitalDefault,eps);

mIKByLeverage         = (nextCapitalLeverage - (1 - ddelta) * EOmegaTerm2 * currentCapitalLeverage) ./ denomLeverage;
mIKByDefaultDistance  = (nextCapitalDefault - (1 - ddelta) * EOmegaTerm2 * currentCapitalDefault) ./ denomDefault;

mIKByLeverage(~isfinite(mIKByLeverage))                       = NaN;
mIKByDefaultDistance(~isfinite(mIKByDefaultDistance))         = NaN;

zeroCapitalLeverage = (denomLeverage <= eps);
zeroCapitalDefault  = (denomDefault  <= eps);

mIKByLeverage(zeroCapitalLeverage)           = NaN;
mIKByDefaultDistance(zeroCapitalDefault)     = NaN;

ikBaseline = 1 - (1 - ddelta) * EOmegaTerm2;

irfIKLowLev    = 100 * (mIKByLeverage(:,1) / ikBaseline - 1);
irfIKHighLev   = 100 * (mIKByLeverage(:,2) / ikBaseline - 1);
irfIKCloseDef  = 100 * (mIKByDefaultDistance(:,1) / ikBaseline - 1);
irfIKFarDef    = 100 * (mIKByDefaultDistance(:,2) / ikBaseline - 1);

figure

h               = gcf;
h.PaperUnits    = 'inches';
h.PaperPosition = [0 0 13 4];

subplot(1,2,1)
hold on
plot(vTime,irfIKLowLev,'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,irfIKHighLev,'linewidth',1.5,'linestyle','--','color',[178/255,34/255,34/255])
plot(vTime,zeros(T,1),'linewidth',1.5,'linestyle','--','color','k')
xlim([1 12])
h        = legend('Low leverage','High leverage');
set(h,'interpreter','latex','location','northeast','fontsize',14)
set(gcf,'color','w')
xlabel('Quarters','interpreter','latex')
ylabel('$\%$ deviation','interpreter','latex')
grid on
title('Investment rate by leverage (I/K)','interpreter','latex','fontsize',14)
hold off

subplot(1,2,2)
hold on
plot(vTime,irfIKCloseDef,'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,irfIKFarDef,'linewidth',1.5,'linestyle','--','color',[178/255,34/255,34/255])
plot(vTime,zeros(T,1),'linewidth',1.5,'linestyle','--','color','k')
xlim([1 12])
h        = legend('Near default','Far from default');
set(h,'interpreter','latex','location','northeast','fontsize',14)
set(gcf,'color','w')
xlabel('Quarters','interpreter','latex')
ylabel('$\%$ deviation','interpreter','latex')
grid on
title('Investment rate by distance to default (I/K)','interpreter','latex','fontsize',14)
hold off

print('../Results/heterogeneity_channels.eps','-depsc')

validLeverage      = ~isnan(mIKByLeverage);
netInvestmentLev   = nextCapitalLeverage - (1 - ddelta) * EOmegaTerm2 * currentCapitalLeverage;
totalCapitalLevRaw = sum(currentCapitalLeverage .* validLeverage,2);
totalCapitalLev    = max(totalCapitalLevRaw,eps);
avgIKSeries        = sum(netInvestmentLev .* validLeverage,2) ./ totalCapitalLev;
avgIKSeries(totalCapitalLevRaw == 0) = NaN;

irfIKAverage       = 100 * (avgIKSeries / ikBaseline - 1);

figure

h               = gcf;
h.PaperUnits    = 'inches';
h.PaperPosition = [0 0 6.5 4];

hold on
plot(vTime,irfIKAverage,'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,zeros(T,1),'linewidth',1.5,'linestyle','--','color','k')
xlim([1 12])
set(gcf,'color','w')
xlabel('Quarters','interpreter','latex')
ylabel('$\%$ deviation','interpreter','latex')
grid on
title('Average investment rate (I/K)','interpreter','latex','fontsize',14)
hold off

print('../Results/heterogeneity_channels_average.eps','-depsc')


%----------------------------------------------------------------
% Simulate panel of firms along transition path
% (run the State file transition_path_regs.do to get regression results)
%----------------------------------------------------------------

vCutoffRange		= [7;9;11;13;15;17;19;21;23;25;27;29];
for iCutoff = 1:12
	tPre    = vCutoffRange(iCutoff) * 4; % select firms who have survived at least tPre years
        transition_path_panel;                           % simulate the transition path
        eval(sprintf('mTransitionPanel_%d = mTransitionPanel;',tPre));
end
% Uncomment to produce new panel simulations with new draws of shocks
%
csvwrite('mTransitionPanel_28.csv',mTransitionPanel_28);
csvwrite('mTransitionPanel_36.csv',mTransitionPanel_36);
csvwrite('mTransitionPanel_44.csv',mTransitionPanel_44);
csvwrite('mTransitionPanel_52.csv',mTransitionPanel_52);
csvwrite('mTransitionPanel_60.csv',mTransitionPanel_60);
csvwrite('mTransitionPanel_68.csv',mTransitionPanel_68);
csvwrite('mTransitionPanel_76.csv',mTransitionPanel_76);
csvwrite('mTransitionPanel_84.csv',mTransitionPanel_84);
csvwrite('mTransitionPanel_92.csv',mTransitionPanel_92);
csvwrite('mTransitionPanel_100.csv',mTransitionPanel_100);
csvwrite('mTransitionPanel_108.csv',mTransitionPanel_108);
csvwrite('mTransitionPanel_116.csv',mTransitionPanel_116);


%----------------------------------------------------------------
% State dependence counterfactual
%----------------------------------------------------------------

%%%
% Compute elasticities across the state space
%%%

% Interpolate capital policy functions over distribution grid
vCapitalPrimeDist       = interpn(mProdGrid,mCashGrid,reshape(mCapitalPrimeSeries(:,1),nProd,nCash),...
                            mStateGridDist(:,1),mStateGridDist(:,2));
vCapitalPrimeDistSS     = interpn(mProdGrid,mCashGrid,reshape(mCapitalPrimeSeries(:,40),nProd,nCash),...
                            mStateGridDist(:,1),mStateGridDist(:,2));

% Compute elasticities
mCapitalElasticity      = reshape(((vCapitalPrimeDist) - vCapitalPrimeDistSS) ./ vCapitalPrimeDistSS,nProd,nCashDist);

% Distribution
mDistribution           = reshape(mDistributionSeries(:,1),nProd,nCashDist);


%%%
% Compute state dependence counterfactual
%%%

% Reference distributions
mDistributionGood				= mDistributionProduction;
vProductivityMarginal			= sum(mDistributionGood,2);

mDistributionBad				= repmat(mDistributionProduction(2,:),[nProd 1]);
vWeight							= sum(mDistributionBad,2) ./ vProductivityMarginal;
mDistributionBad				= mDistributionBad ./ repmat(vWeight,[1 nCashDist]);

% Compute aggregate capital accumulation in each scenario
nRuns						 	= 20;
vWeight					     	= linspace(0,1,nRuns)';
vKResponse				     	= zeros(nRuns,1);
vRiskyConstrained		     	= zeros(nRuns,1);
vMeanCash				     	= zeros(nRuns,1);
vKResponse2			         	= zeros(nRuns,1);

for iRun	= 1 : nRuns
	mDistribution					= vWeight(iRun,1) * mDistributionBad + (1 - vWeight(iRun,1)) * mDistributionGood;
	vKResponse(iRun,1)				= sum((vCapitalPrimeDist - vCapitalPrimeDistSS) .* mDistribution(:));
	vRiskyConstrained(iRun,1)		= sum((vDebtPriceProduction < bbeta) .* mDistribution(:));
	vMeanCash(iRun,1)				= sum(mStateGridDist(:,2) .* mDistribution(:));
end

% Save statistics cited in the paper
vKResponseNormalized				= vKResponse ./ vKResponse(1,1);
vRiskyConstrainedNormalized			= vRiskyConstrained ./ vRiskyConstrained(1,1);
vMeanCashNormalized					= vMeanCash ./ vMeanCash(1,1);

rows 			  = {'Avg net investment response','Average net worth','Fraction risky constrained'};
columns			  = {'Bad','Medium'};
column1 		  = {vKResponseNormalized(nRuns,1);vMeanCashNormalized(nRuns,1);vRiskyConstrainedNormalized(nRuns,1)};
column2 		  = {vKResponseNormalized(nRuns/2,1);vMeanCashNormalized(nRuns/2,1);vRiskyConstrainedNormalized(nRuns/2,1)};
state_dependence_table = table(column1,column2,'RowNames',rows,'VariableNames',columns)
writetable(state_dependence_table,'../Results/state_dependence_table_het_agent_prices.xls','WriteRowNames',true)


%----------------------------------------------------------------
% Compute decomposition of transition paths
%----------------------------------------------------------------

% Real rate
transition_path_real_rate;
mCapitalElasticity_real_rate		= mCapitalElasticity;

% Capital price
transition_path_capital_price;
mCapitalElasticity_capital_price	= mCapitalElasticity;

% Other prices
transition_path_others;
mCapitalElasticity_others			= mCapitalElasticity;

% Plot decomposition of elasticities
fig = figure;

h               = gcf;
h.PaperUnits    = 'inches';
h.PaperPosition = [0 0 13 4];

left_color = [8/255,62/255,118/255];
right_color = [.5 .24 .5];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

subplot(1,2,1)
hold on
yyaxis left
plot(vCashGridDist,mCapitalElasticity_real_rate(ceil(nProd/3),:),'linewidth',1.5,'linestyle','-.','color',[8/255,62/255,118/255])
plot(vCashGridDist,mCapitalElasticity_capital_price(ceil(nProd/3),:),'linewidth',1.5,'linestyle','--','marker','x','color',[0 .4 0])
plot(vCashGridDist,mCapitalElasticity_others(ceil(nProd/3),:),'linewidth',1.5,'linestyle',':','color',[178/255,34/255,34/255])
xlabel('Net worth','interpreter','latex')
ylabel('Semi-Elasticity','interpreter','latex')
h	 = legend('Real rate only','Capital price only','Others');
set(h,'interpreter','latex','location','northeast','fontsize',14)
xlim([0.08 4])
yyaxis right
plot(vCashGridDist,mDistributionProduction(ceil(nProd/3),:) / sum(mDistributionProduction(ceil(nProd/3),:)),'linewidth',1.5,'linestyle','--','color',[.5,.24,.5])
ylabel('Mass of firms')
title('Low Productivity','interpreter','latex','fontsize',14)
set(gcf,'color','w')
grid on
hold off

subplot(1,2,2)
hold on
yyaxis left
plot(vCashGridDist,mCapitalElasticity_real_rate(ceil(2*nProd/3),:),'linewidth',1.5,'linestyle','-.','color',[8/255,62/255,118/255])
plot(vCashGridDist,mCapitalElasticity_capital_price(ceil(2*nProd/3),:),'linewidth',1.5,'linestyle','--','marker','x','color',[0 .4 0])
plot(vCashGridDist,mCapitalElasticity_others(ceil(2*nProd/3),:),'linewidth',1.5,'linestyle',':','color',[178/255,34/255,34/255])
xlabel('Net worth','interpreter','latex')
ylabel('Semi-Elasticity','interpreter','latex')
xlim([0.48 4])
yyaxis right
plot(vCashGridDist,mDistributionProduction(ceil(2*nProd/3),:) / sum(mDistributionProduction(ceil(2*nProd/3),:)),'linewidth',1.5,'linestyle','--','color',[.5,.24,.5])
ylabel('Mass of firms')
title('High Productivity','interpreter','latex','fontsize',14)
set(gcf,'color','w')
grid on
hold off

print('../Results/decomposition_elasticities_byprod.eps','-depsc')



