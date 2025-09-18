function [vOmegaGrid,vOmegaWeights]	= compute_quadrature_rule(ssigmaOmega,nOmega)

global omegaMin omegaMax

%%%
% Compute evenly-spaced grid
%%%

vOmegaGrid       = linspace(omegaMin,omegaMax,nOmega)';


%%%
% Compute weights over that grid
%%%

vPlus 					= zeros(nOmega,1);
    vPlus(nOmega,1) 		= 1e9;
    vPlus(1:nOmega-1,1) 	= vOmegaGrid(2:nOmega);
vMinus					= zeros(nOmega,1);
    vMinus(1,1) 			= -1e9;
    vMinus(2:nOmega,1) 		= vOmegaGrid(1:nOmega-1,1);

vPlusCutoff				= .5 * (vOmegaGrid + vPlus);
vMinusCutoff			= .5 * (vOmegaGrid + vMinus);

vOmegaWeights			= normcdf(vPlusCutoff,0,ssigmaOmega) - normcdf(vMinusCutoff,0,ssigmaOmega);
