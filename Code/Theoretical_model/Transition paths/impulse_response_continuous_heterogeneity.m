function irfResults = impulse_response_continuous_heterogeneity(mCapitalPanelShock,mCapitalPanelBase,mDebtPanel,mCashPanel,mDefaultCutoffPanel,mInSampleShock,mInSampleBase,tPre,varargin)
%IMPULSE_RESPONSE_CONTINUOUS_HETEROGENEITY Compute log-capital IRFs using continuous heterogeneity measures.
%   irfResults = IMPULSE_RESPONSE_CONTINUOUS_HETEROGENEITY(mCapitalPanelShock,
%   mCapitalPanelBase,mDebtPanel,mCashPanel,mDefaultCutoffPanel,mInSampleShock,
%   mInSampleBase,tPre) calculates the impulse-response functions of the log
%   capital stock (expressed in percentage-point deviations relative to the
%   shocked-economy baseline) when heterogeneity is summarized by continuous
%   (demeaned) leverage and distance-to-default measures. The baseline
%   corresponds to the pre-shock period in the economy facing the monetary
%   shock, so no counterfactual path is subtracted.
%
%   Additional name-value pair arguments:
%       'Horizon'      - Number of post-shock quarters to report (default 12).
%       'ResultsDir'   - Directory where Excel outputs are stored
%                        (default ../Results relative to this file).
%
%   The function returns a structure with the computed impulse responses,
%   including the regression slope coefficients, the implied effect of a
%   one-unit increase in the heterogeneity variable, and the effect of a
%   one-standard-deviation increase (using the cross-sectional dispersion at
%   the classification period), as well as a handle to the generated figure.

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
    parse(parser,varargin{:});

    horizonRequested = round(parser.Results.Horizon);

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

    inSampleShockWindow = (mInSampleShock(:, timeIndices) == 1);
    inSampleBaseWindow = (mInSampleBase(:, timeIndices) == 1);
    validWindow = inSampleShockWindow & inSampleBaseWindow;

    logCapitalShock(~validWindow) = NaN;

    baselineShock = logCapitalShock(:,1);

    deltaLogCapital = logCapitalShock - baselineShock;

    quarters = (0:horizon)';
    slopeLeverage = NaN(numel(quarters),1);
    slopeDistance = NaN(numel(quarters),1);
    effectLeveragePerUnit = NaN(numel(quarters),1);
    effectDistancePerUnit = NaN(numel(quarters),1);

    denom = max(mCapitalPanelShock, 1e-8);
    leverageMeasure = mDebtPanel ./ denom;
    leverageMeasure(~isfinite(leverageMeasure)) = NaN;
    leverageMeasure((mInSampleShock == 0) | (mInSampleBase == 0)) = NaN;
    leverageMean = mean(leverageMeasure, 2, 'omitnan');
    leverageDemeaned = leverageMeasure - leverageMean;

    distanceMeasure = mCashPanel - mDefaultCutoffPanel;
    distanceMeasure(~isfinite(distanceMeasure)) = NaN;
    distanceMeasure((mInSampleShock == 0) | (mInSampleBase == 0)) = NaN;
    distanceMean = mean(distanceMeasure, 2, 'omitnan');
    distanceDemeaned = distanceMeasure - distanceMean;

    xLeverage = leverageDemeaned(:,classificationIdx);
    xDistance = distanceDemeaned(:,classificationIdx);
    inSampleClass = (mInSampleShock(:,classificationIdx) == 1) & (mInSampleBase(:,classificationIdx) == 1);
    xLeverage(~inSampleClass) = NaN;
    xDistance(~inSampleClass) = NaN;

    leverageStd = std(xLeverage,'omitnan');
    distanceStd = std(xDistance,'omitnan');

    for h = 0:horizon
        colIdx = h + 1;
        yCurrent = deltaLogCapital(:,colIdx);

        validLeverage = ~isnan(yCurrent) & ~isnan(xLeverage);
        if any(validLeverage)
            yL = yCurrent(validLeverage);
            xL = xLeverage(validLeverage);
            xL = xL - mean(xL,'omitnan');
            yL = yL - mean(yL,'omitnan');
            denomL = sum(xL.^2,'omitnan');
            if denomL > 0
                slopeLeverage(colIdx) = sum(xL .* yL,'omitnan') / denomL;
                effectLeveragePerUnit(colIdx) = slopeLeverage(colIdx);
            end
        end

        validDistance = ~isnan(yCurrent) & ~isnan(xDistance);
        if any(validDistance)
            yD = yCurrent(validDistance);
            xD = xDistance(validDistance);
            xD = xD - mean(xD,'omitnan');
            yD = yD - mean(yD,'omitnan');
            denomD = sum(xD.^2,'omitnan');
            if denomD > 0
                slopeDistance(colIdx) = sum(xD .* yD,'omitnan') / denomD;
                effectDistancePerUnit(colIdx) = slopeDistance(colIdx);
            end
        end
    end

    slopeLeverage = 100 * slopeLeverage;
    slopeDistance = 100 * slopeDistance;
    effectLeveragePerUnit = 100 * effectLeveragePerUnit;
    effectDistancePerUnit = 100 * effectDistancePerUnit;

    effectLeveragePerStd = leverageStd .* effectLeveragePerUnit;
    effectDistancePerStd = distanceStd .* effectDistancePerUnit;

    tblLeverage = table(quarters, slopeLeverage, effectLeveragePerUnit, effectLeveragePerStd, ...
        'VariableNames', {'Quarter','Slope','EffectPerUnit','EffectPerStd'});
    tblDistance = table(quarters, slopeDistance, effectDistancePerUnit, effectDistancePerStd, ...
        'VariableNames', {'Quarter','Slope','EffectPerUnit','EffectPerStd'});

    writetable(tblLeverage, fullfile(resultsDir, 'IRF_Leverage_continuous.xlsx'));
    writetable(tblDistance, fullfile(resultsDir, 'IRF_dd_continuous.xlsx'));

    fig = figure('Color','w');

    subplot(1,2,1);
    hold on;
    plot(quarters, effectLeveragePerStd, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', [0.0 0.45 0.74]);
    yline(0,'--','Color',[0.2 0.2 0.2]);
    grid on;
    xlabel('Trimestres');
    ylabel({'Variación log capital', ...
        '(p.p. vs. línea base, por 1 desv. estándar)'});
    title('Apalancamiento (centrado por firma)');
    hold off;

    subplot(1,2,2);
    hold on;
    plot(quarters, effectDistancePerStd, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', [0.85 0.33 0.10]);
    yline(0,'--','Color',[0.2 0.2 0.2]);
    grid on;
    xlabel('Trimestres');
    ylabel({'Variación log capital', ...
        '(p.p. vs. línea base, por 1 desv. estándar)'});
    hold off;

    irfResults = struct();
    irfResults.quarters = quarters;
    irfResults.leverage.slope = slopeLeverage;
    irfResults.leverage.effectPerUnit = effectLeveragePerUnit;
    irfResults.distance.slope = slopeDistance;
    irfResults.leverage.effectPerStd = effectLeveragePerStd;
    irfResults.leverage.stdAtClassification = leverageStd;
    irfResults.distance.effectPerUnit = effectDistancePerUnit;
    irfResults.distance.effectPerStd = effectDistancePerStd;
    irfResults.distance.stdAtClassification = distanceStd;
    irfResults.meanDifference = 100 * mean(deltaLogCapital, 1, 'omitnan')';
    irfResults.figure = fig;

end