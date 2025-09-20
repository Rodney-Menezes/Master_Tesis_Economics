%%%
% Declare global variables
%%%

global ttheta nnu ddelta rrhoProd ssigmaProd aalpha ppsi_0 ppiExit kInitial b0 ...
  pphiPrice ggamma pphiInflation pphiOutput pphiCapital qSS bbeta nSS cchi ssigmaM rrhoM ...
  mmuEnt ssigmaEnt ssigma pSS kRepSS yRepSS k0 A tthetaHat cRepSS muRepSS wRepSS inflationSS ...
  cchi ssigmaOmega pphiPrice pphiS ssigmaS rrhoS

%%%
% Load in calibrated parameters
%%%

if option_calibration == 1  % use parameters from calibration routine (used to create heatmap)

    ssigmaProd      = vParms(1);
    ppsi_0          = vParms(2);
    kInitial        = vParms(3);
    mmuEnt          = vParms(4);
    aalpha          = vParms(5);
    ssigmaOmega     = vParms(6);

else

    load parameters_calibration.mat 
    ssigmaProd      = vParms(1);
    ppsi_0          = vParms(2);
    kInitial        = vParms(3);
    mmuEnt          = vParms(4);
    aalpha          = vParms(5);
    ssigmaOmega     = vParms(6);
    
end


% Parametros exogenos del riesgo pais:
pphiS = 0.4;
ssigmaS = 0.02;
rrhoS   = 0.8;

% Exogenously fix persistence of productivity shocks
rrhoProd    = 0.90;


%%%
% Fixed parameters
%%%

% Heterogeneous production firms
ttheta 			 = .33;								% capital coefficient
nnu 			 = .67;								% labor coefficient
ddelta           = .025;							% capital depreciation rate

% Firm lifecycle
ppiExit          = (.1) / (4);			 % exogenous exit rate (quarterly); hit 8.7% annual exit rate given a 3% annual default rate. 
b0               = 0;
mmuEnt			 = -mmuEnt * ssigmaProd / sqrt(1 - rrhoProd); % mean of entrant distribution
ssigmaEnt		 = ssigmaProd / sqrt(1 - rrhoProd);	% SD of entrant distribution

% New Keynesian block
ggamma          = 10;							    % elasticity of substitution
pphiPrice		= 90;							    % slope of NKPC = pphi / (ggamma - 1)
pphiInflation   = 1.25;						        % coefficient on inflation in Taylor rule
pphiOutput      = .5;							    % coefficient on output in Taylor rule

% Capital good producers
pphiCapital     = 21;							    % coefficient in capital adjustment cost; 1 / pphi = price elasticity w.r.t. I/K
qSS             = 1;						        % steady state price of capital = 1

% Representative household
bbeta		    = 0.99;								% discount factor (quarterly)
nSS				= 0.6;								% steady state employment rate.
cchi            = 1.72;                                % disuility of labor (will be recalibrated)
ssigma          = 1;                                % coefficient of relative risk aversion

% NUEVO: Sensibilidad al riesgo soberano
pphiS           = 0.10;    % valor inicial (ajustar según calibración)

% Si modelas shock AR(1) de s_t:
rrhoS           = 0.8;     % persistencia aproximada
ssigmaS         = 0.02;    % desviación estándar

lambda0         = 5;        % coeficiente base de apalancamiento  
alpha           = 2;        % sensibilidad del apalancamiento al riesgo s_t  


% Aggregate shocks
ssigmaM 	    = .0025;							% SD of monetary policy shock
rrhoM			= .5;							    % persistence of monetary policy shock


%%%
% Parameters that depend on representative agent steady state
%%%

% Solve for representative agent steady state
pSS				 = (ggamma - 1) / ggamma; % relative price of production firms' output
kRepSS			 = ((pSS * ttheta * (nSS ^ nnu)) / (qSS * ((1 / bbeta) - (1 - ddelta)))) ^ ...
                    (1 / (1 - ttheta)); % capital stock
yRepSS			 = (kRepSS ^ ttheta) * (nSS ^ nnu); % output

% Set parameters
ppsi_0			 = ppsi_0 * yRepSS;					% rescale fixed operating cost
k0               = kInitial * kRepSS; % initial capital stock of new entrants

% Some useful quantities in the computation
A				 = (pSS ^ (1 / (1 - nnu))) * (nnu ^ (nnu / (1 - nnu)) -...
					   nnu ^ (1 / (1 - nnu))); 		 % common term in profit function
tthetaHat		 = ttheta / (1 - nnu);               % curvature in profit function

% Additional quantities from representative agent steady state
cRepSS           = yRepSS - ddelta * kRepSS - ppsi_0;% consumption
muRepSS			 = cRepSS ^ (-ssigma);				 % marginal utility
wRepSS			 = nnu * pSS * yRepSS / nSS;         % real wage
inflationSS	     = 1;                                % steady state inflation
cchi             = wRepSS * (cRepSS ^(-1));          % disuility of labor supply in rep agent model
