function irfResults = impulse_response_heterogeneity(mCapitalPanelShock,mCapitalPanelBase,mDebtPanel,mCashPanel,mDefaultCutoffPanel,mInSampleShock,mInSampleBase,tPre,varargin)
%IMPULSE_RESPONSE_HETEROGENEITY Compute log-capital IRFs for leverage and default-distance heterogeneity.
%   irfResults = IMPULSE_RESPONSE_HETEROGENEITY(mCapitalPanelShock, mCapitalPanelBase,
%   mDebtPanel, mCashPanel, mDefaultCutoffPanel, mInSampleShock, mInSampleBase, tPre)
%   computes impulse-response functions for the log capital stock (expressed
%   in percentage-point deviations relative to a no-shock baseline) along two
%   heterogeneity channels (leverage and distance-to-default). The baseline
%   corresponds to the same economy and initial distribution simulated with
%   the monetary shock turned off. The function reports the response for a
%   horizon of up to 12 quarters (or as many quarters as are available) and
%   returns both the series for each heterogeneity group and a handle to the
%   generated figure.
%%
%   Additional name-value pair arguments:
%       'Horizon'      - Number of post-shock quarters to report (default 12).
%       'ResultsDir'   - Directory where Excel outputs are stored
%                        (default ../Results relative to this file).
%       'Quantiles'    - Two-element vector with lower/upper quantiles used
%                        to split the sample (default [0.25 0.75]).
%
%   The function returns a structure with the computed impulse responses
%   and a handle to the generated figure.

    narginchk(8, inf);

    panelSize = size(mCapitalPanelShock);
    if ~isequal(size(mCapitalPanelBase), panelSize) || ~isequal(size(mDebtPanel), panelSize) || ...
            ~isequal(size(mCashPanel), panelSize) || ~isequal(size(mDefaultCutoffPanel), panelSize)
        error('All panel inputs must share the same dimensions.');
    end
    if ~isequal(size(mInSampleShock), panelSize) || ~isequal(size(mInSampleBase), panelSize)
        error('In-sample indicators must share the same dimensions as the panels.');
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
    epsilon = 1e-12;

    timeIndices = tPre + (0:horizon);
    logCapitalShock = log(max(mCapitalPanelShock(:, timeIndices), epsilon));
    logCapitalBase = log(max(mCapitalPanelBase(:, timeIndices), epsilon));

    inSampleShockWindow = (mInSampleShock(:, timeIndices) == 1);
    inSampleBaseWindow = (mInSampleBase(:, timeIndices) == 1);
    validWindow = inSampleShockWindow & inSampleBaseWindow;

    logCapitalShock(~validWindow) = NaN;
    logCapitalBase(~validWindow) = NaN;

    baselineShock = logCapitalShock(:,1);
    baselineBase = logCapitalBase(:,1);

    deltaShock = logCapitalShock - baselineShock;
    deltaBase = logCapitalBase - baselineBase;
    diffLogCapital = deltaShock - deltaBase;

    quarters = (0:horizon)';
    irfLeverageLow = NaN(numel(quarters),1);
    irfLeverageHigh = NaN(numel(quarters),1);
    irfDistanceClose = NaN(numel(quarters),1);
    irfDistanceFar = NaN(numel(quarters),1);

    validLeverage = (mInSampleShock(:,classificationIdx) == 1) & (mInSampleBase(:,classificationIdx) == 1);
    validLeverage = validLeverage & isfinite(mDebtPanel(:,classificationIdx)) & (mCapitalPanelShock(:,classificationIdx) > 0);
    if ~any(validLeverage)
        error('No valid observations to compute leverage groups at the classification period.');
    end
    leverageMeasure = mDebtPanel(:,classificationIdx) ./ max(mCapitalPanelShock(:,classificationIdx),1e-8);
    leverageMeasure(~validLeverage) = NaN;
    leverageQuantiles = quantile(leverageMeasure(validLeverage), quantiles);
    groupLowLeverage = leverageMeasure <= leverageQuantiles(1);
    groupHighLeverage = leverageMeasure >= leverageQuantiles(2);
    groupLowLeverage(~validLeverage) = false;
    groupHighLeverage(~validLeverage) = false;

    distanceMeasure = mCashPanel(:,classificationIdx) - mDefaultCutoffPanel(:,classificationIdx);
    validDistance = (mInSampleShock(:,classificationIdx) == 1) & (mInSampleBase(:,classificationIdx) == 1) & ~isnan(distanceMeasure);
    if ~any(validDistance)
        error('No valid observations to compute distance-to-default groups at the classification period.');
    end
    distanceQuantiles = quantile(distanceMeasure(validDistance), quantiles);
    groupCloseDefault = distanceMeasure <= distanceQuantiles(1);
    groupFarDefault = distanceMeasure >= distanceQuantiles(2);
    groupCloseDefault(~validDistance) = false;
    groupFarDefault(~validDistance) = false;

    for h = 0:horizon
        colIdx = h + 1;
        timeIdx = tPre + h;

        leverageLowMask = groupLowLeverage & (mInSampleShock(:,timeIdx) == 1) & (mInSampleBase(:,timeIdx) == 1);
        leverageHighMask = groupHighLeverage & (mInSampleShock(:,timeIdx) == 1) & (mInSampleBase(:,timeIdx) == 1);
        distanceCloseMask = groupCloseDefault & (mInSampleShock(:,timeIdx) == 1) & (mInSampleBase(:,timeIdx) == 1);
        distanceFarMask = groupFarDefault & (mInSampleShock(:,timeIdx) == 1) & (mInSampleBase(:,timeIdx) == 1);

        if any(leverageLowMask)
            irfLeverageLow(colIdx,1) = mean(diffLogCapital(leverageLowMask, colIdx), 'omitnan');
        end
        if any(leverageHighMask)
            irfLeverageHigh(colIdx,1) = mean(diffLogCapital(leverageHighMask, colIdx), 'omitnan');
        end
        if any(distanceCloseMask)
            irfDistanceClose(colIdx,1) = mean(diffLogCapital(distanceCloseMask, colIdx), 'omitnan');
        end
        if any(distanceFarMask)
            irfDistanceFar(colIdx,1) = mean(diffLogCapital(distanceFarMask, colIdx), 'omitnan');
        end
    end

    irfLeverageLow = 100 * irfLeverageLow;
    irfLeverageHigh = 100 * irfLeverageHigh;
    irfDistanceClose = 100 * irfDistanceClose;
    irfDistanceFar = 100 * irfDistanceFar;

    tblLeverage = table(quarters, irfLeverageLow, irfLeverageHigh, ...
        'VariableNames', {'Quarter','BajoApalancamientoLogK','AltoApalancamientoLogK'});
    tblDistance = table(quarters, irfDistanceClose, irfDistanceFar, ...
        'VariableNames', {'Quarter','CercaDefaultLogK','LejosDefaultLogK'});

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
    yylabel({'Variación log capital', ...
        '(p.p. vs. base sin shock)'});
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
    ylabel({'Variación log capital', ...
        '(p.p. vs. base sin shock)'});
    title('Canal distancia al default');
    legend({'Cerca del default','Lejos del default'},'Location','best');
    hold off;

    irfResults = struct();
    irfResults.quarters = quarters;
    irfResults.leverage.low = irfLeverageLow;
    irfResults.leverage.high = irfLeverageHigh;
    irfResults.distance.close = irfDistanceClose;
    irfResults.distance.far = irfDistanceFar;
    irfResults.meanDifference = 100 * mean(diffLogCapital, 1, 'omitnan')';
    irfResults.figure = fig;
    irfResults.resultsDir = resultsDir;
    irfResults.quantiles = quantiles;
end