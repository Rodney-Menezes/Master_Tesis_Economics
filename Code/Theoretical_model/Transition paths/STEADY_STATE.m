% Computes and analyzes steady state
%
% Pablo Ottonello and Thomas Winberry
% First draft: September 29nd, 2017
% This draft: June 27th, 2020

% Housekeeping
clear all
close all

%----------------------------------------------------------------
% Set options for what is to be done
%----------------------------------------------------------------

global option_continuous_optimization option_continuous_optimization_VFI option_calibration

option_continuous_optimization			= 1;			% use continuous optimizer once discrete VFI has converged
option_continuous_optimization_VFI		= 0;
option_calibration 						= 0;


%----------------------------------------------------------------
% Set parameters
%----------------------------------------------------------------

set_parameters_model;
set_parameters_numerical;


%----------------------------------------------------------------
% Compute steady state real wage
%----------------------------------------------------------------

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


%%%
% Compute steady state objects
%%%

% Compute objects
core_steady_state;
steady_state_aggregates;

% Calibrate disutility of labor supply
cchi     = wageSS * (aggregateConsumption ^ (-ssigma));


%%%
% Calibration targets and other analysis
%%%

% Compute quarterly lifecycle dynamics
steady_state_lifecycle;

% Compute calibration targets using annual panel simulation
steady_state_calibration_stats;


%----------------------------------------------------------------
% Replicate results from paper
%----------------------------------------------------------------

%%%
% Calibrated parameters
%%%

rows             = {'SD TFP','Operating cost','k0','Entrants mean shift','Recovery rate','SD capital quality','Exit shock'};
columns          = {'Calibrated_value'};
column1          = [ssigmaProd;ppsi_0;k0;mmuEnt / (ssigmaProd / sqrt(1 - rrhoProd));aalpha;ssigmaOmega;ppiExit];
parameters_table = table(column1,'RowNames',rows,'VariableNames',columns)
writetable(parameters_table,'../Results/calibrated_parameters.xls','WriteRowNames',true)


%%%
% Calibration targets %% 333333333   Corregir la edad aqui.3333333
%%%

rows 			  = {'SD(i/k)','Mean gross leverage ratio','Annualized default rate','Emp. share of age 1 firms',...
					'Emp. share of age 2-10 firms','Emp. share of age >= 10 firms','Share of age 1 firms',...
					'Share of age 2 firms','Frac with positive debt','Mean exit rate'};
columns			  = {'Model';'Data'};
column1 		  = [sdIKAnnual;averageGrossLeverage;400*meanDefault;vEmploymentShareAnnual(1,1);sum(vEmploymentShareAnnual(2:5,1));...
					1 - sum(vEmploymentShareAnnual(1:5,1));vExtensiveMarginAnnual(1,1);vExtensiveMarginAnnual(2,1);...
					fracPosDebt;4*massEntrants];
column2			  = [0.02;0.01;0.08;0.04;0.0133;0.0828;0.0178;0.088;0.021;0.06];
calibration_targets_table = table(column1,column2,'RowNames',rows,'VariableNames',columns)
writetable(calibration_targets_table,'../Results/calibration_targets.xls','WriteRowNames',true)


%%%
% Untargetted statistics %% 333333333   Corregir la edad aqui.3333333
%%%

rows 			  = {'Annualized credit spread','Mean net leverage ratio','aggregate TFP','Fraction unconstrained',...
					'Fraction Risky Constrained','Fraction Risk-Free Constrained'};
columns			  = {'Model';'Data'};
column1 		  = [400*meanSpread;averageNetLeverage;aggregateTFP;fractionUnconstrained;fractionRiskyConstrained;...
					1 - fractionUnconstrained - fractionRiskyConstrained];
column2			  = [4.12;0.517;0.00;0.00;1.00;0.000422];
calibration_untargetted_table = table(column1,column2,'RowNames',rows,'VariableNames',columns)
writetable(calibration_untargetted_table,'../Results/calibration_untargetted_stats.xls','WriteRowNames',true)


%%%
% Compustat sample
%%%

rows 			  = {'Frac with positive debt','Mean gross leverage ratio','Mean net leverage ratio','Employment ratio','Age ratio'};
columns			  = {'Model';'Data'};
column1			  = [fracPosDebtCompustat;avgGrossLeverageCompustat;avgNetLeverageCompustat;employmentRatioCompustat;ageRatioCompustat];
column2			  = [0.968;0.341;0.231;1;5];
compustat_table   = table(column1,column2,'RowNames',rows,'VariableNames',columns)
writetable(compustat_table,'../Results/compustat_sample.xls','WriteRowNames',true)


%%%
% Lifecycle dynamics
%%%

% In the data
vMeanGrowthData 	= [.075;.165;.096;.153;.156;.099;.167;.075];

% Aggregate model and data categories
vCategoriesAggregated = [-2;-.2;-.05;.05;.2;2];
vMeanGrowthAggregated	= [vMeanGrowth(1,1);vMeanGrowth(2,1);sum(vMeanGrowth(3:6));vMeanGrowth(7,1);vMeanGrowth(8,1)];
vMeanGrowthDataAggregated	= [vMeanGrowthData(1,1);vMeanGrowthData(2,1);sum(vMeanGrowthData(3:6));vMeanGrowthData(7,1);vMeanGrowthData(8,1)];
vBar = [vMeanGrowthAggregated(1,1) vMeanGrowthDataAggregated(1,1)];
for iCat = 2:5
	vBar = [vBar; vMeanGrowthAggregated(iCat,1) vMeanGrowthDataAggregated(iCat,1)];
end

figure

h               = gcf;
h.PaperUnits    = 'inches';
h.PaperPosition = [0 0 13 4];

subplot(1,2,1)
b = bar(vBar);
b(1).FaceColor = [8/255,62/255,118/255];
b(2).FaceColor = [178/255,34/255,34/255];
set(gca, 'XTickLabel', {'(-2.0,-0.2]','(-0.2,-0.05]','(-0.05,0.05)','(0.05,0.2]','(0.2,2.0)'})
set(gcf,'color','w')
legend('Model','Data')
title('Distribution of growth rates, model vs. data','interpreter','latex','fontsize',14)
grid on
hold off

vTime = linspace(1,20,20)';

subplot(1,2,2)
hold on
plot(vTime(2:end),vGrowthRates(2:end) + 0.1,'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime(2:end),.1*ones(size(vGrowthRates(2:end))),'linewidth',1.5,'linestyle','--','color','k')
xlabel('Age in years','interpreter','latex')
ylabel('Average growth rate','interpreter','latex')
xlim([2 20])
set(gcf,'color','w')
grid on
title('Age-growth profile','interpreter','latex','fontsize',14)
hold off

print('../Results/lifecycle_dynamics.eps','-depsc')


%%%
% Steady state decision rules
%%%

mBorrowingDist		= reshape(interpn(mProdGrid,mCashGrid,reshape(vDebtPriceOptimal,nProd,nCash) .* ...
					   mDebtPrime,mStateGridDist(:,1),mStateGridDist(:,2)),nProd,nCashDist);

fig = figure;
left_color = [8/255,62/255,118/255];
right_color = [.5 .24 .5];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

h               = gcf;
h.PaperUnits    = 'inches';
h.PaperPosition = [0 0 9 3];

subplot(1,2,1)
hold on
yyaxis left
plot(vCashGridDist,mCapitalPrimeDist(ceil(nProd/3),:),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vCashGridDist,mBorrowingDist(ceil(nProd/3),:),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vCashGridDist,mDividendsDist(ceil(nProd/3),:),'linewidth',1.5,'linestyle','-','color',[.3,.78,.5])
plot(vCashGridDist,zeros(nCashDist,1),'linewidth',1.5,'linestyle','--','color','k')
xlabel('Net worth','interpreter','latex')
ylabel('Decisions','interpreter','latex')
xlim([cashMin 10])
yyaxis right
plot(vCashGridDist,mDistributionProduction(ceil(nProd/3),:) / sum(mDistributionProduction(ceil(nProd/3),:)),'linewidth',1.5,'linestyle','--','color',[.5,.24,.5])
ylabel('Mass of firms')
title('Low productivity','interpreter','latex','fontsize',14)
set(gcf,'color','w')
h	 = legend('Capital','Borrowing','Dividends');
set(h,'interpreter','latex','location','northeast','fontsize',10)
grid on
hold off

subplot(1,2,2)
hold on
yyaxis left
plot(vCashGridDist,mCapitalPrimeDist(ceil(2*nProd/3),:),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vCashGridDist,mBorrowingDist(ceil(2*nProd/3),:),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vCashGridDist,mDividendsDist(ceil(2*nProd/3),:),'linewidth',1.5,'linestyle','-','color',[.1,.78,.5])
plot(vCashGridDist,zeros(nCashDist,1),'linewidth',1.5,'linestyle','--','color','k')
xlabel('Net worth','interpreter','latex')
ylabel('Decisions','interpreter','latex')
xlim([cashMin 10])
yyaxis right
plot(vCashGridDist,mDistributionProduction(ceil(2*nProd/3),:) / sum(mDistributionProduction(ceil(2*nProd/3),:)),'linewidth',1.5,'linestyle','--','color',[.5,.24,.5])
ylabel('Mass of firms')
title('High productivity','interpreter','latex','fontsize',14)
set(gcf,'color','w')
grid on
hold off

print('../Results/decision_rules.png','-dpng')


%%%
% Marginal propensity to invest
%%%

% Capital elasticity
mCapitalElasticity   = (mCapitalPrimeDist(:,2:end) - mCapitalPrimeDist(:,1:end-1)) ./ ...
												(mCashGridDist(:,2:end) - mCashGridDist(:,1:end-1));
mCapitalElasticity  = [mCapitalElasticity mCapitalElasticity(:,end-1)];

% Compute debt elasticity
mDebtElasticity     = (mDebtPrimeDist(:,2:end) - mDebtPrimeDist(:,1:end-1)) ./ ...
												(mCashGridDist(:,2:end) - mCashGridDist(:,1:end-1));
mDebtElasticity(mDebtElasticity == -Inf) = 0;           % sometimes b' = 0
mDebtElasticity     = [mDebtElasticity mDebtElasticity(:,end-1)];

% Compute average elasticity w.r.t. productivity
vCapitalElasticityAverage   = sum(repmat(vProdErgodic,[1 nCashDist]) .* mCapitalElasticity,1)';
vDebtElasticityAverage      = sum(repmat(vProdErgodic,[1 nCashDist]) .* mDebtElasticity,1)';
vDistributionAverage        = sum(mDistributionProduction,1)';

% Average w.r.t. productivity
fig = figure;
left_color = [8/255,62/255,118/255];
right_color = [.5 .24 .5];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

h               = gcf;
h.PaperUnits    = 'inches';
h.PaperPosition = [0 0 8 5];

hold on
yyaxis left
plot(vCashGridDist,vCapitalElasticityAverage,'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vCashGridDist,zeros(nCashDist,1),'linewidth',1.5,'linestyle','--','color','k')
xlabel('Net worth','interpreter','latex')
ylabel('Marginal Propensity to Invest','interpreter','latex')
xlim([vCashGridDist(2) 5])
yyaxis right
plot(vCashGridDist,vDistributionAverage,'linewidth',1.5,'linestyle','--','color',[.5,.24,.5])
ylabel('Mass of firms')
title('Average Marginal Propensity to Invest','interpreter','latex','fontsize',14)
set(gcf,'color','w')
grid on
hold off

print('../Results/marginal_propensity_invest.eps','-depsc')


%%%
% Leverage vs. net worth
%%%

% Compute leverage and credit spreads over whole distribution
mLeverageDist				 = reshape(vDebtPrimeDist ./ max(k0,vCapitalPrimeDist),nProd,nCashDist);
mSpreadDist					 = reshape((1 ./ vDebtPriceProduction) - (1 / bbeta),nProd,nCashDist);

% Average over productivity
vLeverageDistAvg		     = sum(repmat(vProdErgodic,[1 nCashDist]) .* mLeverageDist,1)';
vSpreadDistAvg			     = sum(repmat(vProdErgodic,[1 nCashDist]) .* mSpreadDist,1)';
vDistributionAverage	     = sum(mDistributionProduction,1)';

fig = figure;
left_color = [8/255,62/255,118/255];
right_color = [.5 .24 .5];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

h               = gcf;
h.PaperUnits    = 'inches';
h.PaperPosition = [0 0 8 5];

hold on
yyaxis left
plot(vCashGridDist,vLeverageDistAvg,'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vCashGridDist,zeros(nCashDist,1),'linewidth',1.5,'linestyle','--','color','k')
xlabel('Net worth','interpreter','latex')
ylabel('Leverage $\frac{b^{\prime}}{k^{\prime}}$','interpreter','latex')
xlim([0.075 5])
yyaxis right
plot(vCashGridDist,vDistributionAverage,'linewidth',1.5,'linestyle','--','color',[.5,.24,.5])
ylabel('Mass of firms')
set(gcf,'color','w')
grid on
hold off

print('../Results/leverage_vs_networth.eps','-depsc')


%%%
% Average lifecycle dynamics
%%%

vTime		  = linspace(1,T_lifecycle,T_lifecycle)';
vAvgSpread    = 400 * ((1 ./ vAvgDebtPrice) - (1 / bbeta));

figure

h               = gcf;
h.PaperUnits    = 'inches';
h.PaperPosition = [0 0 9 6];

subplot(2,3,1)
hold on
plot(vTime,vCapitalAggregate,'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
xlim([1 50])
xlabel('Age in quarters','interpreter','latex')
ylabel('Capital')
set(gcf,'color','w')
grid on
title('Capital','interpreter','latex','fontsize',14)
hold off

subplot(2,3,2)
hold on
plot(vTime,vDebtAggregate,'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
xlim([1 50])
xlabel('Age in quarters','interpreter','latex')
ylabel('Debt')
set(gcf,'color','w')
grid on
title('Debt','interpreter','latex','fontsize',14)
hold off

subplot(2,3,3)
hold on
plot(vTime,vLeverageAggregate,'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
xlim([1 50])
xlabel('Age in quarters','interpreter','latex')
ylabel('Leverage')
set(gcf,'color','w')
grid on
title('Leverage','interpreter','latex','fontsize',14)
hold off

subplot(2,3,4)
hold on
plot(vTime,vAvgProductivityAggregate,'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
xlim([1 50])
xlabel('Age in quarters','interpreter','latex')
ylabel('Productivity')
set(gcf,'color','w')
grid on
title('Average Productivity','interpreter','latex','fontsize',14)
hold off

subplot(2,3,5)
hold on
plot(vTime,vLaborAggregate,'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
xlim([1 50])
xlabel('Age in quarters','interpreter','latex')
ylabel('Employment')
set(gcf,'color','w')
grid on
title('Employment','interpreter','latex','fontsize',14)
hold off

subplot(2,3,6)
hold on
plot(vTime,vAvgSpread,'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
xlim([1 50])
xlabel('Age in quarters','interpreter','latex')
ylabel('Annualized p.p.')
set(gcf,'color','w')
grid on
title('Credit Spread','interpreter','latex','fontsize',14)
hold off

print('../Results/lifecycle_dynamics_2.eps','-depsc')


%%%
% Summary statistics of model's simulated panel
% (NB: other entries are in steady_state_panel_analysis.do)
%%%

% Autocorrelation
vAutoCorrIKAnnual               = zeros(tAnnual,1);
vAutoCorrIKAnnual_full 			= zeros(tAnnual,1);
for t = 2 : tAnnual
    vAutocorrIKAnnual(t,1)      = corr(mInvestmentRateAnnual(mBalancedPanelIndicator(:,tAnnual),t),...
                                    mInvestmentRateAnnual(mBalancedPanelIndicator(:,tAnnual),t-1));
	vAutocorrIKAnnualFull(t,1)  = corr(mInvestmentRateAnnual(logical(mInSampleAnnual(:,t)),t),...
								    mInvestmentRateAnnual(logical(mInSampleAnnual(:,t)),t-1));
end
autocorrIKAnnual                = mean(vAutocorrIKAnnual);
autocorrIKAnnualFull 			= mean(vAutocorrIKAnnualFull);

% Mean and standard deviation
meanIKAnnual   					= mean(mInvestmentRateAnnual(logical(mBalancedPanelIndicator)));
meanIKAnnualFull   				= mean(mInvestmentRateAnnual(logical(mInSampleAnnual)));
sdIKAnnual						= std(mInvestmentRateAnnual(logical(mBalancedPanelIndicator)));
sdIKAnnualFull					= std(mInvestmentRateAnnual(logical(mInSampleAnnual)));

rows 			  = {'E(i/k)','SD(i/k)','autocorr(i/k)'};
columns			  = {'Data';'Model_selected';'Model_full'};
column1 		  = [0.194;1.79;-0.00281];
column2 		  = [meanIKAnnual;sdIKAnnual;autocorrIKAnnual];
column3 		  = [meanIKAnnualFull;sdIKAnnualFull;autocorrIKAnnualFull];
panel_stats_table = table(column1,column2,column3,'RowNames',rows,'VariableNames',columns)
writetable(panel_stats_table,'../Results/panel_stats.xls','WriteRowNames',true)


%----------------------------------------------------------------
% Compute the heatmaps to assess identification
%----------------------------------------------------------------

%%%
% Compute the Jacobian of moments w.r.t. parameters
%%%

% Preallocate
nMoments              = 8;
nParms                = 6;
mJacobian             = zeros(nMoments,nParms);

% Load in calibrated parameters
load parameters_calibration.mat
vParmsInit                = vParms;

% Compute function value at initial guess
vG                    = calibration_targets_individual(vParmsInit);

% Compute Jacobian by finite differences
hh                    = 1e-3;
for iParms = 1 : nParms
  vParmsDiff            = vParmsInit;
  vParmsDiff(iParms)    = vParms(iParms) + hh;
  vGDiff                = calibration_targets_individual(vParmsDiff);
  mJacobian(:,iParms)   = (vGDiff - vG) / hh;
end

%%%
% Compute the sensitivity measure
%%%

% Assume identifiy weighting matrix
mWeighting              = eye(nMoments);
mLambda                 = inv(mJacobian' * mWeighting * mJacobian) * (mJacobian' * mWeighting);

% Rescale so that changes are 1% change in each moment
mMoments                = repmat([33.7;34.4;3;80.22;2.8;21.1;10.5;8.1]',[nParms 1]);
mParms                  = repmat(vParmsInit,[1 nMoments]);
mLambdaScaled           = mLambda .* (mMoments ./ mParms);

% Rescale jacobian
mJacobianScaled         = mJacobian .* (mParms' ./ mMoments');


%%%
% Visualize the heat map
%%%

xvalues     = {'$\sigma$','$\xi$','$k_{0}$','$m$','$\alpha$','$\sigma_{\omega}$'};
yvalues     = {'$\sigma(\frac{i}{k})$','E[default rate]','Frac($b > 0$)','E[gross leverage]','$N_{1}/N$','$N_{1-10}/N$','$M_{1}/M$','$M_{2}/M$'};

% Colormap
n                 = 500;
indexValue        = 0;
bottomcolor       =[30/255 80/255 178/255];
indexColor        =[1 1 1];
topColor          =[178/255 60/255 60/255];

%%% BELOW IS PASTED FROM ONLINE
% Colormap for Jacobian
scaleJacobian     = max(abs(mJacobianScaled(:)));
largest           = scaleJacobian; smallest = -scaleJacobian;
index             = knnsearch(linspace(smallest,largest,n)',0);

% Create color map ranging from bottom color to index color
% Multipling number of points by 100 adds more resolution
customCMap1       = [linspace(bottomcolor(1),indexColor(1),index)',...
            			linspace(bottomcolor(2),indexColor(2),index)',...
            			linspace(bottomcolor(3),indexColor(3),index)'];

% Create color map ranging from index color to top color
% Multipling number of points by 100 adds more resolution
customCMap2       = [linspace(indexColor(1),topColor(1),(n-index))',...
            		linspace(indexColor(2),topColor(2),(n-index))',...
            		linspace(indexColor(3),topColor(3),(n-index))'];
colormap_jacobian = [customCMap1;customCMap2]; % Combine colormaps

% Colormap for Lambda
scaleLambda       = max(abs(mLambdaScaled(:)));
largest           = scaleLambda; smallest = -scaleLambda;
index             = knnsearch(linspace(smallest,largest,n)',0);

% Create color map ranging from bottom color to index color
% Multipling number of points by 100 adds more resolution
customCMap1       = [linspace(bottomcolor(1),indexColor(1),index)',...
            linspace(bottomcolor(2),indexColor(2),index)',...
            linspace(bottomcolor(3),indexColor(3),index)'];

% Create color map ranging from index color to top color
% Multipling number of points by 100 adds more resolution
customCMap2       = [linspace(indexColor(1),topColor(1),(n-index))',...
            linspace(indexColor(2),topColor(2),(n-index))',...
            linspace(indexColor(3),topColor(3),(n-index))'];
colormap_lambda   = [customCMap1;customCMap2]; % Combine colormaps

figure

h               = gcf;
h.PaperUnits    = 'inches';
h.PaperPosition = [0 0 9 6];

subplot(2,1,1)
h1 = heatmap(mJacobianScaled,xvalues,yvalues, '%0.2f','TickAngle',45,'Colormap',colormap_jacobian,'Colorbar','true','textcolor',[.1 .1 .1],'MinColorValue',-3,'MaxColorValue',3,'ShowAllTicks',true,'gridlines',':');
set(gcf,'color','w')
title('Elasticity of moments w.r.t. parameters','interpreter','latex','fontsize',14)

subplot(2,1,2)
h2 = heatmap(mLambdaScaled',xvalues,yvalues, '%0.2f','TickAngle',45,'Colormap',colormap_lambda,'Colorbar','true','textcolor',[.1 .1 .1],'MinColorValue',-5,'MaxColorValue',5,'ShowAllTicks',true,'gridlines',':');
set(gcf,'color','w')
title('Elasticity of parameters w.r.t. moments','interpreter','latex','fontsize',14)

print('../Results/identification_heatmaps.eps','-depsc')
