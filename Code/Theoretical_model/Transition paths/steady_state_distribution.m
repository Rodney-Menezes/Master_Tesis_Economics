%----------------------------------------------------------------
% Compute transition matrix mapping firms this period
% into firms next period.
% CONVENTION: rows = distribution state variables next period;
% columns = distribution state variables this period
%----------------------------------------------------------------

% In general, the implied net worth n'(z',z,x) is not on the grid of cash on hand n.
% We find the nearest two gridpoints and assign of fraction of firms with n'(z',z,n),
% weighted by distance.

% Compute weighting matrices for the two nearest gridpoints
[vIndicesAbove,vIndicesBelow,vWeightAbove,vWeightBelow] = ...
	compute_weighting_matrices(vCashPrimeDist,vCashGridDist);

% Preallocate some matrices
mTransitionAbove    = zeros(nCashDist,nShocks * nStateDist);
mTransitionBelow    = zeros(nCashDist,nShocks * nStateDist);

% Construct the transition matrices for each of the points
for x = 1 : nCashDist
	mTransitionBelow(x,vIndicesBelow == x)	= vWeightBelow(vIndicesBelow == x);
	mTransitionAbove(x,vIndicesAbove == x)	= vWeightAbove(vIndicesAbove == x);
end


%%%
% Integrate out the capital quality shocks
%%%

aTransition				= reshape(mTransitionBelow + mTransitionAbove,nCashDist,nOmega,nProd * nStateDist);
mTransition				= reshape(sum(repmat(reshape(vOmegaWeights,1,nOmega,1),[nCashDist 1 ...
							nProd * nStateDist]) .* aTransition,2),nCashDist,nProd * nStateDist);


%%%
% Combine the above transition matrix with the transition matrix for productivity shocks
%%%

% Permute transition matrix for cash on hand to have the correct dimensionality
mTransitionCash			= reshape(permute(reshape(mTransition,nCashDist,nProd,nStateDist),...
							[2 1 3]),nStateDist,nStateDist);

% Combine with transition matrix for productivity
mTransitionProd			= reshape(repmat(reshape(mProdTransition',nProd,1,nProd,1),[1 nCashDist 1 nCashDist]),...
												nStateDist,nStateDist);
mTransitionIncumbents	= sparse(mTransitionCash .* mTransitionProd);

%%%
% Indicator variable for whether (z,n) will default for continuing firms
%%%

mDistContContinue		= (mCashGridDist >= repmat(vDefaultCutoff,[1 nCashDist]));
vDistContContinue		= mDistContContinue(:);


%----------------------------------------------------------------
% Compute vector mapping new entrants next period
% into firms next period.
% CONVENTION: rows = distribution state variables next period;
% columns = productivity draw upon entry next period
%----------------------------------------------------------------

% Compute weighting matrices for nearest point in the state space
[vIndicesAbove,vIndicesBelow,vWeightAbove,vWeightBelow] = ...
	compute_weighting_matrices(vCashImpliedEnt,vCashGridDist);

% Preallocate some matrices
mTransitionAbove    = zeros(nCashDist,nShocks);
mTransitionBelow    = zeros(nCashDist,nShocks);

% Construct the transition matrix
for x = 1 : nCashDist
    mTransitionBelow(x,vIndicesBelow == x)  =  vWeightBelow(vIndicesBelow == x);
	mTransitionAbove(x,vIndicesAbove == x)  =  vWeightAbove(vIndicesAbove == x);
end

% Integrate out the capital quality shocks
aTransition			     = reshape(mTransitionBelow + mTransitionAbove,nCashDist,nOmega,nProd);
mTransition			     = reshape(sum(repmat(reshape(vOmegaWeights,1,nOmega,1),[nCashDist 1 nProd]) .* ...
							aTransition,2),nCashDist,nProd);

% Combine into vector
mTransitionEntrants		 = mTransition';		% want productivity in first dimension
vDistributionEntrants	 = reshape(mTransitionEntrants .* repmat(vProdDistEnt,[1 nCashDist]),nStateDist,1);


%----------------------------------------------------------------
% Compute stationary distribution
%----------------------------------------------------------------

% Initialize iteration counters
err                     = 100;
iteration               = 0;
vDistribution			= ones(nStateDist,1) / nStateDist;

% Do the iteration
while err > 1e-16 && iteration < 1e3

	% Choose mass of entrants to ensure mass of firms in production = 1
	vContinue				= ppiExit * (mStateGridDist(:,2) >= 0) + (1 - ppiExit) * vDistContContinue;
	massSurviveIncumbents	= (1 - ppiExit) * sum((mTransitionIncumbents * (vDistContContinue .* vDistribution)) .* vContinue);
	massSurviveEntrants		= sum(vContinue .* vDistributionEntrants);
	massEntrants			= (1 - massSurviveIncumbents) / massSurviveEntrants;

	% Update distribution
	vDistributionNew		= (1 - ppiExit) * mTransitionIncumbents * (vDistContContinue .* vDistribution) + ... % incumbent firms
								massEntrants * vDistributionEntrants;											 % new entrants

	% Update iteration
	err                 	= max(abs(vDistributionNew - vDistribution));
	iteration           	= iteration + 1;
	vDistribution			= vDistributionNew;

end


%----------------------------------------------------------------
% Save output from the iteration
%----------------------------------------------------------------

% Distribution of firms at the beginning of the period
mDistribution				= reshape(vDistribution,nProd,nCashDist);

% Distribution of firms in production
mDistributionProduction		= reshape(vContinue .* vDistribution,nProd,nCashDist);
