function irfResults = impulse_response_continuous_heterogeneity(mInvestmentPanel,mCapitalPanel,mDebtPanel,mCashPanel,mDefaultCutoffPanel,mInSample,tPre,varargin)
%IMPULSE_RESPONSE_CONTINUOUS_HETEROGENEITY Compute IRFs using continuous heterogeneity measures.
%   irfResults = IMPULSE_RESPONSE_CONTINUOUS_HETEROGENEITY(mInvestmentPanel,
%   mCapitalPanel,mDebtPanel,mCashPanel,mDefaultCutoffPanel,mInSample,tPre)
%   calculates the cumulative impulse-response functions of investment
%   rates (investment-to-capital) expressed as accumulated percentage-point
%   deviations from a post-shock baseline ("línea base") when heterogeneity
%   is summarized by
%   continuous (demeaned) leverage and distance-to-default measures rather
%   than discrete quantile groups.
%
%   Additional name-value pair arguments:
%       'Horizon'      - Number of post-shock quarters to report (default 12).
%       'ResultsDir'   - Directory where Excel outputs are stored
%                        (default ../Results relative to this file).
%
%   The function returns a structure with the computed impulse responses,
%   including the regression slope coefficients and the implied effect of a
%   one-unit increase in the heterogeneity variable, as well as a handle to
%   the generated figure.

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

    denom = max(mCapitalPanel, 1e-8);
    investmentRate = mInvestmentPanel ./ denom;
    investmentRate(~isfinite(investmentRate)) = NaN;
    investmentRate(mInSample == 0) = NaN;

    leverageMeasure = mDebtPanel ./ denom;
    leverageMeasure(~isfinite(leverageMeasure)) = NaN;
    leverageMeasure(mInSample == 0) = NaN;
    leverageMean = mean(leverageMeasure, 2, 'omitnan');
    leverageDemeaned = leverageMeasure - leverageMean;

    distanceMeasure = mCashPanel - mDefaultCutoffPanel;
    distanceMeasure(~isfinite(distanceMeasure)) = NaN;
    distanceMeasure(mInSample == 0) = NaN;
    distanceMean = mean(distanceMeasure, 2, 'omitnan');
    distanceDemeaned = distanceMeasure - distanceMean;

    xLeverage = leverageDemeaned(:,classificationIdx);
    xDistance = distanceDemeaned(:,classificationIdx);

    postShockIdx = tPre + (1:horizon);
    baselinePath = mean(investmentRate(:, postShockIdx), 1, 'omitnan');
    if any(isnan(baselinePath))
        error('Unable to compute the aggregate baseline investment path after the shock.');
    end

    investmentDeviation = bsxfun(@minus, investmentRate(:, postShockIdx), baselinePath);
    cumulativeDeviation = cumsum(investmentDeviation, 2);

    quarters = (1:horizon)';
    slopeLeverage = NaN(horizon,1);
    slopeDistance = NaN(horizon,1);
    effectLeveragePerUnit = NaN(horizon,1);
    effectDistancePerUnit = NaN(horizon,1);

    for h = 1:horizon
        yCurrent = cumulativeDeviation(:,h);

        validLeverage = ~isnan(yCurrent) & ~isnan(xLeverage);
        if any(validLeverage)
            yL = yCurrent(validLeverage);
            xL = xLeverage(validLeverage);
            xL = xL - mean(xL,'omitnan');
            yL = yL - mean(yL,'omitnan');
            denomL = sum(xL.^2,'omitnan');
            if denomL > 0
                slopeLeverage(h) = sum(xL .* yL,'omitnan') / denomL;
                effectLeveragePerUnit(h) = slopeLeverage(h);
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
                slopeDistance(h) = sum(xD .* yD,'omitnan') / denomD;
                effectDistancePerUnit(h) = slopeDistance(h);
            end
        end
    end

    slopeLeverage = 100 * slopeLeverage;
    slopeDistance = 100 * slopeDistance;
    effectLeveragePerUnit = 100 * effectLeveragePerUnit;
    effectDistancePerUnit = 100 * effectDistancePerUnit;

    tblLeverage = table(quarters, slopeLeverage, effectLeveragePerUnit, ...
        'VariableNames', {'Quarter','Slope','EffectPerUnit'});
    tblDistance = table(quarters, slopeDistance, effectDistancePerUnit, ...
        'VariableNames', {'Quarter','Slope','EffectPerUnit'});

    writetable(tblLeverage, fullfile(resultsDir, 'IRF_Leverage_continuous.xlsx'));
    writetable(tblDistance, fullfile(resultsDir, 'IRF_dd_continuous.xlsx'));

    fig = figure('Color','w');

    subplot(1,2,1);
    hold on;
    plot(quarters, effectLeveragePerUnit, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', [0.0 0.45 0.74]);
    yline(0,'--','Color',[0.2 0.2 0.2]);
    grid on;
    xlabel('Trimestres');
    ylabel({'Variación acumulada inv./capital', ...
        '(p.p. vs. línea base, por 1 unidad)'});
    title('Apalancamiento (centrado por firma)');
    hold off;

    subplot(1,2,2);
    hold on;
    plot(quarters, effectDistancePerUnit, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', [0.85 0.33 0.10]);
    yline(0,'--','Color',[0.2 0.2 0.2]);
    grid on;
    xlabel('Trimestres');
    ylabel({'Variación acumulada inv./capital', ...
        '(p.p. vs. línea base, por 1 unidad)'});
    title('Distancia al default (centrada por firma)');
    hold off;

    irfResults = struct();
    irfResults.quarters = quarters;
    irfResults.leverage.slope = slopeLeverage;
    irfResults.leverage.effectPerUnit = effectLeveragePerUnit;
    irfResults.distance.slope = slopeDistance;
    irfResults.distance.effectPerUnit = effectDistancePerUnit;
    irfResults.baselinePath = baselinePath(:);
    irfResults.figure = fig;

end