%%%
% Aggregates next period (K_{t+1} and B_{t+1})
%%%

aggregateCapital		= (1 - ppiExit) * sum(vCapitalPrimeDist .* vDistContContinue .* mDistribution(:)) + ...
								massEntrants * sum(vCapitalPrimeDist .* vDistContContinue .* vDistributionEntrants);
aggregateDebt			= (1 - ppiExit) * sum(vDebtPrimeDist .* vDistContContinue .* mDistribution(:)) + ...
								massEntrants * massSurviveEntrants * b0;

%%%
% Average leverage
%%%

averageLeverage			= (1 - ppiExit) * sum((vDebtPrimeDist ./ max(0.1,vCapitalPrimeDist)) .* ...
								vDistContContinue .* mDistribution(:)) + ...
								massEntrants * massSurviveEntrants * (b0 / k0);
averageGrossLeverage	= (1 - ppiExit) * sum((max(0,vDebtPrimeDist) ./ (max(.1,vCapitalPrimeDist) + max(0,-vDebtPrimeDist))) .* ...
								vDistContContinue .* mDistribution(:)) + ...
								massEntrants * massSurviveEntrants * (b0 / k0);


%%%
% Aggregate capital by category
%%%

% Unconstrained cutoff
mDistUnconstrained 	= (mCashGridDist >= repmat(vUnconstrainedCutoff,[1 nCashDist]));
vDistUnconstrained	= mDistUnconstrained(:);

% For unconstrained firms
aggregateCapitalUnconstrained		= (1 - ppiExit) * sum(vCapitalPrimeDist .* vDistUnconstrained .* mDistribution(:));

% For constrained firms
aggregateCapitalConstrained			= (1 - ppiExit) * sum(vCapitalPrimeDist .* (1 - vDistUnconstrained) .* vDistContContinue .* mDistribution(:));


%%%
% Aggregates this period (N_{t}, Y_{t}, Mass_{t}, I_{t}, C_{t}, TFP_{t}, fraction unconstrained_{t})
%%%

% Default indicators
mContinueIncumbent			= ppiExit * mContinueExitDist + (1 - ppiExit) * mContinueContDist;
vContinueEntrant			= ppiExit * vContinueExitEnt + (1 - ppiExit) * vContinueContEnt;

%%%
% Employment
%%%

% Incumbents
aLaborIntegrandIncumbent	= reshape((((pSS * nnu * exp(mProdPrimeDist) .* ((exp(mOmegaPrimeDist) .* ...
								 mCapitalPrimePrimeDist) .^ ttheta)) / wage) .^ (1 / (1 - nnu))) .* ...
								 mContinueIncumbent,nShocks,nProd,nCashDist);

mLaborIntegrandIncumbent  = zeros(nProd,nCashDist);
for iProd = 1 : nProd
    mLaborIntegrandIncumbent(iProd,:)   = mShocksTransition(iProd,:) * squeeze(aLaborIntegrandIncumbent(:,iProd,:));
end

% New entrants
vLaborIntegrandEntrant		= (((pSS * nnu * exp(mShocksGrid(:,1)) .* ((exp(mShocksGrid(:,2)) * k0) .^ ttheta)) / ...
								wage) .^ (1 / (1 - nnu))) .* vContinueEntrant;

% Aggregate
aggregateLabor           	= (1 - ppiExit) * sum(mLaborIntegrandIncumbent(:) .* mDistribution(:) .* vDistContContinue) + ... % incumbents
                          		massEntrants * sum(vLaborIntegrandEntrant .* vDistEnt);               % new entrant

%%%
% Output
%%%

% Incumbents
aOutputIntegrandIncumbent	= reshape(exp(mProdPrimeDist) .* ((exp(mOmegaPrimeDist) .* mCapitalPrimePrimeDist) .^ ttheta) .* ...
								(reshape(aLaborIntegrandIncumbent,nShocks,nStateDist) .^ nnu) .* mContinueIncumbent,nShocks,nProd,nCashDist);

mOutputIntegrandIncumbent = zeros(nProd,nCashDist);
for iProd = 1 : nProd
    mOutputIntegrandIncumbent(iProd,:) = mShocksTransition(iProd,:) * squeeze(aOutputIntegrandIncumbent(:,iProd,:));
end

% New entrants
vOutputIntegrandEntrant		= exp(mShocksGrid(:,1)) .* ((exp(mShocksGrid(:,2)) * k0) .^ ttheta) .* ...
								(vLaborIntegrandEntrant .^ nnu) .* vContinueEntrant;

% Aggregate
aggregateOutput           = (1 - ppiExit) * sum(mOutputIntegrandIncumbent(:) .* mDistribution(:) .* vDistContContinue) + ... % incumbents
                            	massEntrants * sum(vOutputIntegrandEntrant .* vDistEnt);               % new entrant

%%%
% The rest
%%%

% Mass of firms in production
aggregateMass					= sum(mDistribution(:) .* vContinue);

% Investment
aggregateInvestment				= (1 - (1 - ddelta) * EOmegaTerm2) * (1 + (massEntrants * k0 / aggregateCapital)) * aggregateCapital;

% Consumption
aggregateConsumption			= aggregateOutput - aggregateInvestment - ppsi_0 * aggregateMass;

% TFP
aggregateTFP					= aggregateOutput / ((aggregateCapital ^ ttheta) * (aggregateLabor ^ nnu));


%%%
% Mass of firms around thresholds (unconstrained and defaults)
%%%

% Mass of firms unconstrained
fractionUnconstrained	= sum(sum((mCashGridDist >= repmat(vUnconstrainedCutoff,[1 nCashDist])) .* ...
							mDistribution));

% Mass of firms that default
fractionDefault			= sum(sum((1 - vContinue) .* mDistribution(:)));
meanDefault				= fractionDefault;

% Mass of firms that exit
otherExits				= ppiExit * sum((mStateGridDist(:,2) >= 0) .* mDistribution(:));
meanExit				= meanDefault + otherExits;

% Credit spread
vDebtPriceProduction = interpn(mProdGrid,mCashGrid,reshape(vDebtPriceOptimal,nProd,nCash),...
                        mStateGridDist(:,1),mStateGridDist(:,2));
meanPrice			 = (1 - ppiExit) * sum(vDebtPriceProduction .* vDistContContinue .* mDistribution(:)) ./ ...
						((1 - ppiExit) * sum(vDistContContinue .* mDistribution(:)));
meanSpread           = (1 / meanPrice) - (1 / bbeta);

% Fraction of risky constrained firms
fractionRiskyConstrained= sum((vDebtPriceProduction < bbeta) .* mDistribution(:));


%%%
% Misc. calibration statistics (checking distribution of leverage)
%%%

% Compute necessary integrands
vGrossLeverageDist		= max(0,vDebtPrimeDist) ./ (max(.1,vCapitalPrimeDist) + max(0,-vDebtPrimeDist));
vNetLeverageDist 		= vDebtPrimeDist ./ (max(.1,vCapitalPrimeDist) + max(0,-vDebtPrimeDist));

% Compute useful statistics
fracPosDebt 			= (1 - ppiExit) * sum((vDebtPrimeDist > 0) .* ...
								vDistContContinue .* mDistribution(:));
averageNetLeverage 		= (1 - ppiExit) * sum(vNetLeverageDist .* ...
								vDistContContinue .* mDistribution(:)) + ...
								massEntrants * massSurviveEntrants * (b0 / k0);
