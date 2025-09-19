function irfResults = impulse_response_heterogeneity(mInvestmentPanel,mCapitalPanel,mDebtPanel,mCashPanel,mDefaultCutoffPanel,mInSample,tPre,varargin)
%IMPULSE_RESPONSE_HETEROGENEITY Compute IRFs for leverage and default-distance heterogeneity.
%   irfResults = IMPULSE_RESPONSE_HETEROGENEITY(mInvestmentPanel, mCapitalPanel,
%   mDebtPanel, mCashPanel, mDefaultCutoffPanel, mInSample, tPre) computes
%   the cumulative impulse responses of investment rates along two
%   heterogeneity channels (leverage and distance-to-default) over a default
%   horizon of 12 quarters using the simulated panel outputs from
%   TRANSITION_PATH_PANEL.
%
%   Additional name-value pair arguments:
%       'Horizon'      - Number of post-shock quarters to report (default 12).
%       'ResultsDir'   - Directory where Excel outputs are stored
%                        (default ../Results relative to this file).
%       'Quantiles'    - Two-element vector with lower/upper quantiles used
%                        to split the sample (default [0.25 0.75]).
%
%   The function returns a structure with the computed impulse responses
%   and a handle to the generated figure.

    narginchk(7, inf);

    panelSize = size(mInvestmentPanel);
    if ~isequal(size(mCapitalPanel), panelSize) || ~isequal(size(mDebtPanel), panelSize) || ...
            ~isequal(size(mCashPanel), panelSize) || ~isequal(size(mDefaultCutoffPanel), panelSize) || ...
            ~isequal(size(mInSample), panelSize)
        error('All panel inputs must share the same dimensions.');
    end

    parser = inputParser;
    parser.FunctionName = mfilename;
    addParameter(parser,'Horizon',12,@(x)isnumeric(x) && isscalar(x) && x > 0);
    addParameter(parser,'ResultsDir','',@(x)ischar(x) || isstring(x));
    addParameter(parser,'Quantiles',[0.25 0.75],@(x)isnumeric(x) && numel(x)==2 && all(x>=0) && all(x<=1) && x(1) < x(2));
    parse(parser,varargin{:});

    horizonRequested = round(parser.Results.Horizon);
    quantiles = sort(parser.Results.Quantiles(:)');

    if isempty(parser.Results.ResultsDir)
        thisFile = mfilename('fullpath');
        thisDir = fileparts(thisFile);
        resultsDir = fullfile(thisDir,'..','Results');
    else
        resultsDir = char(parser.Results.ResultsDir);
    end

    if ~exist(resultsDir,'dir')
        mkdir(resultsDir);
    end

    totalPeriods = panelSize(2);
    if tPre >= totalPeriods
        error('tPre must be smaller than the number of simulated quarters.');
    end

    horizon = min(horizonRequested, totalPeriods - tPre);
    if horizon < 1
        error('Not enough post-shock periods available to compute the IRFs.');
    end

    classificationIdx = min(tPre + 1, totalPeriods);
    baselineIdx = tPre;

    denom = max(mCapitalPanel, 1e-8);
    investmentRate = mInvestmentPanel ./ denom;
    investmentRate(~isfinite(investmentRate)) = NaN;
    investmentRate(mInSample == 0) = NaN;

    validLeverage = (mInSample(:,classificationIdx) == 1) & isfinite(investmentRate(:,classificationIdx));
    if ~any(validLeverage)
        error('No valid observations to compute leverage groups at the classification period.');
    end
    leverageMeasure = mDebtPanel(:,classificationIdx) ./ max(mCapitalPanel(:,classificationIdx),1e-8);
    leverageMeasure(~validLeverage) = NaN;
    leverageQuantiles = quantile(leverageMeasure(validLeverage), quantiles);
    groupLowLeverage = leverageMeasure <= leverageQuantiles(1);
    groupHighLeverage = leverageMeasure >= leverageQuantiles(2);

    distanceMeasure = mCashPanel(:,classificationIdx) - mDefaultCutoffPanel(:,classificationIdx);
    validDistance = (mInSample(:,classificationIdx) == 1) & ~isnan(distanceMeasure);
    if ~any(validDistance)
        error('No valid observations to compute distance-to-default groups at the classification period.');
    end
    distanceQuantiles = quantile(distanceMeasure(validDistance), quantiles);
    groupCloseDefault = distanceMeasure <= distanceQuantiles(1);
    groupFarDefault = distanceMeasure >= distanceQuantiles(2);

    baselineMask = (mInSample(:,baselineIdx) == 1);
    baselineLowLeverage = mean(investmentRate(groupLowLeverage & baselineMask, baselineIdx), 'omitnan');
    baselineHighLeverage = mean(investmentRate(groupHighLeverage & baselineMask, baselineIdx), 'omitnan');
    baselineCloseDefault = mean(investmentRate(groupCloseDefault & baselineMask, baselineIdx), 'omitnan');
    baselineFarDefault = mean(investmentRate(groupFarDefault & baselineMask, baselineIdx), 'omitnan');

    if any(isnan([baselineLowLeverage, baselineHighLeverage, baselineCloseDefault, baselineFarDefault]))
        error('Unable to compute baseline investment rates for one or more groups.');
    end

    quarters = (1:horizon)';
    irfLeverageLow = NaN(horizon,1);
    irfLeverageHigh = NaN(horizon,1);
    irfDistanceClose = NaN(horizon,1);
    irfDistanceFar = NaN(horizon,1);

    for h = 1:horizon
        idx = tPre + h;
        leverageLowMask = groupLowLeverage & (mInSample(:,idx) == 1);
        leverageHighMask = groupHighLeverage & (mInSample(:,idx) == 1);
        distanceCloseMask = groupCloseDefault & (mInSample(:,idx) == 1);
        distanceFarMask = groupFarDefault & (mInSample(:,idx) == 1);

        irfLeverageLow(h,1) = mean(investmentRate(leverageLowMask, idx), 'omitnan') - baselineLowLeverage;
        irfLeverageHigh(h,1) = mean(investmentRate(leverageHighMask, idx), 'omitnan') - baselineHighLeverage;
        irfDistanceClose(h,1) = mean(investmentRate(distanceCloseMask, idx), 'omitnan') - baselineCloseDefault;
        irfDistanceFar(h,1) = mean(investmentRate(distanceFarMask, idx), 'omitnan') - baselineFarDefault;
    end

    irfLeverageLow = cumsum(irfLeverageLow);
    irfLeverageHigh = cumsum(irfLeverageHigh);
    irfDistanceClose = cumsum(irfDistanceClose);
    irfDistanceFar = cumsum(irfDistanceFar);

    irfLeverageLow = 100 * irfLeverageLow;
    irfLeverageHigh = 100 * irfLeverageHigh;
    irfDistanceClose = 100 * irfDistanceClose;
    irfDistanceFar = 100 * irfDistanceFar;

    tblLeverage = table(quarters, irfLeverageLow, irfLeverageHigh, ...
        'VariableNames', {'Quarter','BajoApalancamientoAcum','AltoApalancamientoAcum'});
    tblDistance = table(quarters, irfDistanceClose, irfDistanceFar, ...
        'VariableNames', {'Quarter','CercaDefaultAcum','LejosDefaultAcum'});

    writetable(tblLeverage, fullfile(resultsDir, 'IRF_Leverage.xlsx'));
    writetable(tblDistance, fullfile(resultsDir, 'IRF_dd.xlsx'));

    fig = figure('Color','w');
    subplot(1,2,1);
    hold on;
    plot(quarters, irfLeverageLow, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', [0.0 0.45 0.74]);
    plot(quarters, irfLeverageHigh, 'LineWidth', 1.5, 'LineStyle', '--', 'Color', [0.85 0.33 0.10]);
    yline(0,'--','Color',[0.2 0.2 0.2]);
    grid on;
    xlabel('Trimestres');
    ylabel('p.p. acumulados respecto al estado estacionario');
    title('Canal de apalancamiento');
    legend({'Bajo apalancamiento','Alto apalancamiento'},'Location','best');
    hold off;

    subplot(1,2,2);
    hold on;
    plot(quarters, irfDistanceClose, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', [0.0 0.45 0.74]);
    plot(quarters, irfDistanceFar, 'LineWidth', 1.5, 'LineStyle', '--', 'Color', [0.85 0.33 0.10]);
    yline(0,'--','Color',[0.2 0.2 0.2]);
    grid on;
    xlabel('Trimestres');
    ylabel('p.p. acumulados respecto al estado estacionario');
    title('Canal distancia al default');
    legend({'Cerca del default','Lejos del default'},'Location','best');
    hold off;

    irfResults = struct();
    irfResults.quarters = quarters;
    irfResults.leverage.low = irfLeverageLow;
    irfResults.leverage.high = irfLeverageHigh;
    irfResults.distance.close = irfDistanceClose;
    irfResults.distance.far = irfDistanceFar;
    irfResults.figure = fig;
    irfResults.resultsDir = resultsDir;
    irfResults.quantiles = quantiles;
end