function [vProdGrid,mProdTransition,vProdErgodic] = compute_tauchen_transition(rrhoProd,ssigmaProd,nProd);

% Declare global variables
global prodMin prodMax

% Compute grid of productivity shocks
vProdGrid       = linspace(prodMin,prodMax,nProd)';

% Compute Markov transition matrix approximating AR(1) using Tauchen (1986) method
vPlus 					= zeros(1,nProd);
    vPlus(1,nProd) 			= 1e9;
    vPlus(1,1:nProd-1) 		= vProdGrid(2:nProd)';
vMinus					= zeros(1,nProd);
    vMinus(1,1) 			= -1e9;
    vMinus(1,2:nProd) 		= vProdGrid(1:nProd-1)';
mPlusCutoff 			= .5 * ones(nProd,1) * (vProdGrid' + vPlus) - rrhoProd * vProdGrid * ones(1,nProd);
    mMinusCutoff 			= .5 * ones(nProd,1) * (vProdGrid' + vMinus) - rrhoProd * vProdGrid * ones(1,nProd);
mProdTransition	 		= normcdf(mPlusCutoff,0,ssigmaProd) - normcdf(mMinusCutoff,0,ssigmaProd);

% Compute implied ergodic distribution
vProdErgodic            = mProdTransition ^ 1e4;    % iterate 1e4 times on transition matrix
vProdErgodic            = vProdErgodic(1,:)';       % extract relevant row
