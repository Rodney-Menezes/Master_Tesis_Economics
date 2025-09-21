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

currentCapitalInc = nansum(currentCapitalLeverage,2);
nextCapitalInc    = nansum(nextCapitalLeverage,2);

denomInc = max(currentCapitalInc,eps);
denomLeverage = max(currentCapitalLeverage,eps);
denomDefault  = max(currentCapitalDefault,eps);

aggregateIKInc        = (nextCapitalInc - (1 - ddelta) * EOmegaTerm2 .* currentCapitalInc) ./ denomInc;
mIKByLeverage         = (nextCapitalLeverage - (1 - ddelta) * EOmegaTerm2 * currentCapitalLeverage) ./ denomLeverage;
mIKByDefaultDistance  = (nextCapitalDefault - (1 - ddelta) * EOmegaTerm2 * currentCapitalDefault) ./ denomDefault;

aggregateIKInc(~isfinite(aggregateIKInc))                       = NaN;
mIKByLeverage(~isfinite(mIKByLeverage))                       = NaN;
mIKByDefaultDistance(~isfinite(mIKByDefaultDistance))         = NaN;

zeroCapitalInc      = (denomInc      <= eps);
zeroCapitalLeverage = (denomLeverage <= eps);
zeroCapitalDefault  = (denomDefault  <= eps);

aggregateIKInc(zeroCapitalInc)                     = NaN;
mIKByLeverage(zeroCapitalLeverage)           = NaN;
mIKByDefaultDistance(zeroCapitalDefault)     = NaN;

ikBaseline = 1 - (1 - ddelta) * EOmegaTerm2;

ikDeviationLeverage        = 100 * (mIKByLeverage - ikBaseline);
ikDeviationDefaultDistance = 100 * (mIKByDefaultDistance - ikBaseline);

ikDeviationLeverage(~isfinite(ikDeviationLeverage))                       = 0;
ikDeviationDefaultDistance(~isfinite(ikDeviationDefaultDistance))         = 0;

irfIKLowLev    = 100 * (mIKByLeverage(:,1) ./ aggregateIKInc - 1);
irfIKHighLev   = 100 * (mIKByLeverage(:,2) ./ aggregateIKInc - 1);
irfIKCloseDef  = 100 * (mIKByDefaultDistance(:,1) ./ aggregateIKInc - 1);
irfIKFarDef    = 100 * (mIKByDefaultDistance(:,2) ./ aggregateIKInc - 1);

%----------------------------------------------------------------
% Interaction between investment and firm characteristics
%----------------------------------------------------------------

% Build series that include steady state benchmarks
leverageSeriesFull          = [vLeverageByLeverageSS; mLeverageByLeverage(1:T,:)];
distanceSeriesFull          = [vDistanceByDefaultDistanceSS; mDistanceByDefaultDistance(1:T,:)];
massSeriesLeverageFull      = massSeriesLeverage;
massSeriesDefaultFull       = massSeriesDefault;

% Use the steady state distribution when weighting the interactions so that
% the model counterpart mirrors the empirical regressions that rely on
% predetermined (pre-shock) firm characteristics.
massBaselineLeverage        = massSeriesLeverageFull(1,:);
massBaselineDefault         = massSeriesDefaultFull(1,:);
leverageBaseline            = leverageSeriesFull(1,:);
distanceBaseline            = distanceSeriesFull(1,:);

interactionLevLevel         = NaN(T,1);
interactionDistLevel        = NaN(T,1);

for t = 1 : T

        ikRowLev    = mIKByLeverage(t,:);
        validLev    = isfinite(massBaselineLeverage) & isfinite(leverageBaseline) & isfinite(ikRowLev);
        weightLev   = sum(massBaselineLeverage(validLev));
        if weightLev > 0
                interactionLevLevel(t,1) = sum(massBaselineLeverage(validLev) .* leverageBaseline(validLev) .* ikRowLev(validLev)) / weightLev;
        end

        ikRowDist    = mIKByDefaultDistance(t,:);
        validDist    = isfinite(massBaselineDefault) & isfinite(distanceBaseline) & isfinite(ikRowDist);
        weightDist   = sum(massBaselineDefault(validDist));
        if weightDist > 0
                interactionDistLevel(t,1) = sum(massBaselineDefault(validDist) .* distanceBaseline(validDist) .* ikRowDist(validDist)) / weightDist;
        end

end

validBaselineLeverage       = isfinite(massBaselineLeverage) & isfinite(leverageBaseline);
weightBaselineLeverage      = sum(massBaselineLeverage(validBaselineLeverage));
if weightBaselineLeverage > 0
        baselineInteractionLeverage = ikBaseline * sum(massBaselineLeverage(validBaselineLeverage) .* leverageBaseline(validBaselineLeverage)) / ...
                                      weightBaselineLeverage;
else
        baselineInteractionLeverage = NaN;
end

validBaselineDefault        = isfinite(massBaselineDefault) & isfinite(distanceBaseline);
weightBaselineDefault       = sum(massBaselineDefault(validBaselineDefault));
if weightBaselineDefault > 0
        baselineInteractionDefault = ikBaseline * sum(massBaselineDefault(validBaselineDefault) .* distanceBaseline(validBaselineDefault)) / ...
                                      weightBaselineDefault;
else
        baselineInteractionDefault = NaN;
end

interactionLevIRF           = NaN(T,1);
interactionDistIRF          = NaN(T,1);

if isfinite(baselineInteractionLeverage) && abs(baselineInteractionLeverage) > eps
        interactionLevIRF   = 100 * (interactionLevLevel / baselineInteractionLeverage - 1);
end
if isfinite(baselineInteractionDefault) && abs(baselineInteractionDefault) > eps
        interactionDistIRF  = 100 * (interactionDistLevel / baselineInteractionDefault - 1);
end

%----------------------------------------------------------------
% Group-to-group comparisons for leverage and default distance
%----------------------------------------------------------------

% Leverage groups: Low vs High
lowLevSeries          = mIKByLeverage(:,1);
highLevSeries         = mIKByLeverage(:,2);
validLevMask          = isfinite(lowLevSeries) & isfinite(highLevSeries) & (lowLevSeries > 0) & (highLevSeries > 0);

rel_low_vs_high       = NaN(size(lowLevSeries));
rel_high_vs_low       = NaN(size(lowLevSeries));
gap_sym_lev           = NaN(size(lowLevSeries));

rel_low_vs_high(validLevMask) = 100 * (lowLevSeries(validLevMask) ./ highLevSeries(validLevMask) - 1);
rel_high_vs_low(validLevMask) = 100 * (highLevSeries(validLevMask) ./ lowLevSeries(validLevMask) - 1);
gap_sym_lev(validLevMask)     = 100 * (highLevSeries(validLevMask) - lowLevSeries(validLevMask)) ./ ...
                                        ((highLevSeries(validLevMask) + lowLevSeries(validLevMask)) / 2);

% Default distance groups: Far vs Near
farDefSeries          = mIKByDefaultDistance(:,1);
nearDefSeries         = mIKByDefaultDistance(:,2);
validDefMask          = isfinite(farDefSeries) & isfinite(nearDefSeries) & (farDefSeries > 0) & (nearDefSeries > 0);

rel_far_vs_near       = NaN(size(farDefSeries));
rel_near_vs_far       = NaN(size(farDefSeries));
gap_sym_def           = NaN(size(farDefSeries));

rel_far_vs_near(validDefMask) = 100 * (farDefSeries(validDefMask) ./ nearDefSeries(validDefMask) - 1);
rel_near_vs_far(validDefMask) = 100 * (nearDefSeries(validDefMask) ./ farDefSeries(validDefMask) - 1);
gap_sym_def(validDefMask)     = 100 * (nearDefSeries(validDefMask) - farDefSeries(validDefMask)) ./ ...
                                        ((nearDefSeries(validDefMask) + farDefSeries(validDefMask)) / 2);

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
% Local projections for transition panels (MATLAB implementation)
%----------------------------------------------------------------

try
        horizons_local_projection   = 0:12;
        shock_parameters_lp = struct('length', 12, 'size', -0.0025, ...
                                     'decay', 0.5, 'scale', -4);

        panel_variable_names = who('mTransitionPanel_*');
        if isempty(panel_variable_names)
                fprintf('No transition panel data available for local projection analysis.\n');
        else
                fprintf('Processing %d transition panel(s) for local projections.\n', numel(panel_variable_names));

                panel_tables_lp = {};
                for iPanel = 1:numel(panel_variable_names)
                        panel_name_lp = panel_variable_names{iPanel};
                        panel_data_lp = eval(panel_name_lp);

                        t_pre_value = NaN;
                        numeric_tokens = regexp(panel_name_lp,'\d+','match');
                        if ~isempty(numeric_tokens)
                                t_pre_value = str2double(numeric_tokens{end});
                        end

                        fprintf('  %s %s (t_{pre} = %d)\n', char(8226), panel_name_lp, t_pre_value);

                        panel_results_table = compute_local_projections_for_panel(panel_data_lp, t_pre_value, ...
                                horizons_local_projection, shock_parameters_lp, panel_name_lp);
                        if ~isempty(panel_results_table)
                                panel_tables_lp{end+1,1} = panel_results_table; %#ok<AGROW>
                        end
                end

                if ~isempty(panel_tables_lp)
                        local_projection_results = vertcat(panel_tables_lp{:});

                        results_directory = fullfile('..','Results');
                        if ~exist(results_directory,'dir')
                                mkdir(results_directory);
                        end

                        results_output_path = fullfile(results_directory,'transition_panel_local_projections_simple.csv');
                        writetable(local_projection_results, results_output_path);
                        fprintf('Local projection estimates saved to: %s\n', results_output_path);

                        summary_table_lp = summarize_local_projection_results(local_projection_results, horizons_local_projection);
                        if isempty(summary_table_lp)
                                fprintf('No valid regression results were produced, skipping the plot.\n');
                        else
                                plot_local_projection_summary(summary_table_lp, horizons_local_projection, results_directory);
                        end
                else
                        fprintf('No valid regression results were produced for the available panels.\n');
                end
        end
catch ME
        warning('Local projection analysis failed: %s', ME.message);
end


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


%% ------------------------------------------------------------------------
%% Local function definitions for the local projection analysis
%% ------------------------------------------------------------------------

function panel_table = compute_local_projections_for_panel(panel_matrix, t_pre, horizons, shock_params, panel_name)

panel_table = table();
if isempty(panel_matrix)
        return
end

in_sample_mask = panel_matrix(:,3) == 1;
panel_matrix = panel_matrix(in_sample_mask,:);
if isempty(panel_matrix)
        return
end

panel_matrix = sortrows(panel_matrix,[1 2]);

firm_id     = panel_matrix(:,1);
quarter_id  = panel_matrix(:,2);
capital     = panel_matrix(:,8);
cash        = panel_matrix(:,9);
debt        = panel_matrix(:,10);

log_capital = NaN(size(capital));
valid_capital = capital > 0;
log_capital(valid_capital) = log(capital(valid_capital));

ratio_mask = capital ~= 0;
leverage = NaN(size(capital));
leverage(ratio_mask) = debt(ratio_mask) ./ capital(ratio_mask);

cash_to_capital = NaN(size(capital));
cash_to_capital(ratio_mask) = cash(ratio_mask) ./ capital(ratio_mask);

leverage_dm = demean_by_group(leverage, firm_id);
cash_to_capital_dm = demean_by_group(cash_to_capital, firm_id);

n_obs = numel(firm_id);
lag_log_capital = NaN(n_obs,1);
lag_leverage_dm = NaN(n_obs,1);
lag_dd_dm       = NaN(n_obs,1);
if n_obs > 1
        previous_idx = (1:(n_obs-1))';
        next_idx = previous_idx + 1;
        same_firm_previous = firm_id(next_idx) == firm_id(previous_idx);
        sequential_previous = (quarter_id(next_idx) - quarter_id(previous_idx)) == 1;
        valid_lag = same_firm_previous & sequential_previous;
        valid_positions = next_idx(valid_lag);

        lag_log_capital(valid_positions) = log_capital(previous_idx(valid_lag));
        lag_leverage_dm(valid_positions) = leverage_dm(previous_idx(valid_lag));
        lag_dd_dm(valid_positions)       = cash_to_capital_dm(previous_idx(valid_lag));
end

max_quarter = max(quarter_id);
shock_series = build_shock_series(max_quarter, t_pre, shock_params);
shock_values = shock_series(quarter_id);

lev_shock = lag_leverage_dm .* shock_values;
dd_shock  = lag_dd_dm .* shock_values;

max_h = max(horizons);
lead_matrix = compute_lead_matrix(log_capital, firm_id, quarter_id, max_h);

measure_names = {'leverage','default_distance'};
interaction_series = {lev_shock, dd_shock};

results_struct = struct('panel_name',{},'t_pre',{},'measure',{},'horizon',{},'coefficient',{},'std_error',{},'n_obs',{});

for m = 1:numel(measure_names)
        current_interaction = interaction_series{m};
        for h = horizons
                response_series = lead_matrix(:,h+1) - lag_log_capital;
                [coef, se, n_reg_obs] = estimate_local_projection(firm_id, quarter_id, response_series, current_interaction);

                results_struct(end+1,1) = struct( ...
                        'panel_name', panel_name, ...
                        't_pre', t_pre, ...
                        'measure', measure_names{m}, ...
                        'horizon', h, ...
                        'coefficient', coef, ...
                        'std_error', se, ...
                        'n_obs', n_reg_obs); %#ok<AGROW>
        end
end

if ~isempty(results_struct)
        panel_table = struct2table(results_struct);
end

end


function shock_series = build_shock_series(max_quarter, t_pre, params)

shock_series = zeros(max_quarter,1);
if isnan(t_pre) || t_pre + 1 > max_quarter
        return
end

for h = 1:params.length
        period = t_pre + h;
        if period > max_quarter
                break
        end
        shock_series(period,1) = params.size * (params.decay .^ (h - 1));
end

shock_series = shock_series * params.scale;

end


function lead_matrix = compute_lead_matrix(log_capital, firm_id, quarter_id, max_h)

n_obs = numel(log_capital);
lead_matrix = NaN(n_obs, max_h + 1);
lead_matrix(:,1) = log_capital;

if max_h == 0 || n_obs == 0
        return
end

for h = 1:max_h
        if h >= n_obs
                break
        end

        idx_base = (1:(n_obs - h))';
        same_firm = firm_id(idx_base) == firm_id(idx_base + h);
        sequential_period = (quarter_id(idx_base + h) - quarter_id(idx_base)) == h;
        valid_entries = same_firm & sequential_period;
        valid_idx = idx_base(valid_entries);

        lead_matrix(valid_idx, h + 1) = log_capital(valid_idx + h);
end

end


function [beta, se, n_obs] = estimate_local_projection(firm_id, quarter_id, response, interaction)

mask = isfinite(response) & isfinite(interaction);
firm_id = firm_id(mask);
quarter_id = quarter_id(mask);
y = response(mask);
x = interaction(mask);

n_obs = numel(y);
if n_obs == 0
        beta = NaN;
        se = NaN;
        return
end

n_firm = max(firm_id);
n_time = max(quarter_id);

n_firm_effective = max(n_firm - 1, 0);
n_time_effective = max(n_time - 1, 0);
n_columns = n_firm_effective + n_time_effective;

if n_columns > 0
        obs_index = (1:n_obs)';
        rows = [];
        cols = [];

        if n_firm_effective > 0
                firm_mask = firm_id > 1;
                rows = [rows; obs_index(firm_mask)]; %#ok<AGROW>
                cols = [cols; firm_id(firm_mask) - 1]; %#ok<AGROW>
        end

        if n_time_effective > 0
                time_mask = quarter_id > 1;
                rows = [rows; obs_index(time_mask)]; %#ok<AGROW>
                cols = [cols; n_firm_effective + quarter_id(time_mask) - 1]; %#ok<AGROW>
        end

        values = ones(numel(rows),1);
        B = sparse(rows, cols, values, n_obs, n_columns);

        tol = 1e-10;
        max_iterations = 1000;

        if isempty(B)
                y_tilde = y - mean(y);
                x_tilde = x - mean(x);
        else
                [coeff_y, flag_y] = lsqr(B, y, tol, max_iterations);
                if flag_y > 1
                        warning('LSQR did not converge for the response variable (flag %d).', flag_y);
                end
                y_tilde = y - B * coeff_y;

                [coeff_x, flag_x] = lsqr(B, x, tol, max_iterations);
                if flag_x > 1
                        warning('LSQR did not converge for the regressor (flag %d).', flag_x);
                end
                x_tilde = x - B * coeff_x;
        end
else
        y_tilde = y - mean(y);
        x_tilde = x - mean(x);
end

denominator = sum(x_tilde .^ 2);
if denominator <= eps || ~isfinite(denominator)
        beta = NaN;
        se = NaN;
        return
end

beta = (x_tilde' * y_tilde) / denominator;
residuals = y_tilde - x_tilde * beta;

cluster_sums = accumarray(firm_id, x_tilde .* residuals, [n_firm, 1], @sum, 0);
variance_beta = (cluster_sums' * cluster_sums) / (denominator ^ 2);

cluster_counts = accumarray(firm_id, 1, [n_firm, 1], @sum, 0);
G = sum(cluster_counts > 0);
K = 1;

if G > 1
        variance_beta = variance_beta * (G / (G - 1));
end
if n_obs > K
        variance_beta = variance_beta * ((n_obs - 1) / (n_obs - K));
end

se = sqrt(variance_beta);

end


function summary_table = summarize_local_projection_results(results_table, horizons)

summary_table = table();
if isempty(results_table)
        return
end

measures = unique(results_table.measure);
summary_entries = struct('measure',{},'horizon',{},'coefficient',{},'std_error',{},'n_obs',{},'n_panels',{},'ci_low',{},'ci_high',{});

for iMeasure = 1:numel(measures)
        measure_key = measures{iMeasure};
        for h = horizons
                mask = strcmp(results_table.measure, measure_key) & results_table.horizon == h & isfinite(results_table.coefficient);
                if ~any(mask)
                        continue
                end

                coeff = weighted_mean_safe(results_table.coefficient(mask), results_table.n_obs(mask));
                se = sqrt(weighted_mean_safe(results_table.std_error(mask) .^ 2, results_table.n_obs(mask)));
                n_obs = sum(results_table.n_obs(mask));
                n_panels = sum(mask);
                ci_low = coeff - 1.645 * se;
                ci_high = coeff + 1.645 * se;

                summary_entries(end+1,1) = struct( ...
                        'measure', measure_key, ...
                        'horizon', h, ...
                        'coefficient', coeff, ...
                        'std_error', se, ...
                        'n_obs', n_obs, ...
                        'n_panels', n_panels, ...
                        'ci_low', ci_low, ...
                        'ci_high', ci_high); %#ok<AGROW>
        end
end

if ~isempty(summary_entries)
        summary_table = struct2table(summary_entries);
end

end


function value = weighted_mean_safe(values, weights)

valid_entries = isfinite(values) & isfinite(weights) & (weights > 0);
if ~any(valid_entries)
        value = NaN;
        return
end

values = values(valid_entries);
weights = weights(valid_entries);

value = sum(values .* weights) / sum(weights);

end


function demeaned_series = demean_by_group(series, group_ids)

if isempty(series)
        demeaned_series = series;
        return
end

if ~isequal(size(series), size(group_ids))
        error('Series and group identifiers must have the same dimensions.');
end

valid_entries = isfinite(series) & isfinite(group_ids);
demeaned_series = NaN(size(series));

if ~any(valid_entries)
        return
end

max_group_id = max(group_ids(valid_entries));
sum_by_group = accumarray(group_ids(valid_entries), series(valid_entries), [max_group_id, 1], @sum, 0);
count_by_group = accumarray(group_ids(valid_entries), 1, [max_group_id, 1], @sum, 0);

group_means = NaN(max_group_id, 1);
positive_counts = count_by_group > 0;
group_means(positive_counts) = sum_by_group(positive_counts) ./ count_by_group(positive_counts);

demeaned_series(valid_entries) = series(valid_entries) - group_means(group_ids(valid_entries));

end


function plot_local_projection_summary(summary_table, horizons, results_directory)

if isempty(summary_table)
        return
end

measure_order = {'leverage','default_distance'};
panel_titles = {'Panel (a): Heterogeneidad por apalancamiento', ...
                'Panel (b): Heterogeneidad por distancia al default'};
panel_colours = [178, 34, 34; 70, 130, 180] / 255;

fig = figure('Color','w','Position',[100 100 960 420]);

for iPanel = 1:numel(measure_order)
        subplot(1, numel(measure_order), iPanel);
        hold on

        data_subset = summary_table(strcmp(summary_table.measure, measure_order{iPanel}), :);
        if isempty(data_subset)
                text(0.5,0.5,'Sin resultados','HorizontalAlignment','center');
                axis off
                hold off
                continue
        end

        [h_values, sort_idx] = sort(data_subset.horizon);
        coeff_values = data_subset.coefficient(sort_idx);
        low_values = data_subset.ci_low(sort_idx);
        high_values = data_subset.ci_high(sort_idx);

        fill([h_values; flipud(h_values)], [low_values; flipud(high_values)], panel_colours(iPanel,:), ...
                'FaceAlpha',0.2,'EdgeColor','none');
        plot(h_values, coeff_values, '-', 'Color', panel_colours(iPanel,:), 'LineWidth', 1.5);
        plot(h_values, coeff_values, 'o', 'Color', panel_colours(iPanel,:), ...
                'MarkerFaceColor', panel_colours(iPanel,:), 'LineWidth', 1.5);

        xlim([min(horizons), max(horizons)]);
        xticks(horizons);
        xlabel('Trimestres','Interpreter','latex');
        ylabel('Efecto acumulado en la inversi\''on','Interpreter','latex');
        title(panel_titles{iPanel},'Interpreter','latex','FontSize',13);
        grid on
        set(gca,'FontSize',12);
        hold off
end

sgtitle({'FHeterogeneidad financiera en la dinamica de la inversion ante un shock monetario expansivo'...
         'con datos simulados con el modelo teorico'}, 'Interpreter','latex','FontSize',14);

plot_path = fullfile(results_directory,'transition_panel_local_projections_simple.png');
print(fig, plot_path, '-dpng', '-r300');
fprintf('Plot saved to: %s\n', plot_path);

end

