%----------------------------------------------------------------
% Get inflation and real interest rate from Taylor rule and EE for bonds
%----------------------------------------------------------------

for t = T:-1:1

	vAggregateInflation(t,1)	= exp((1 / pphiInflation) * (log((vAggregateMarginalUtility(t,1) ./ ...
									vAggregateMarginalUtility(t+1,1)) * vAggregateInflation(t+1,1)) - vEpsilon_m(t,1)));	% Combines Taylor rule and EE for bonds

	vAggregateR_nom(t,1)        = exp(log(1 / bbeta) + pphiInflation * log(vAggregateInflation(t,1)) + ...
									vEpsilon_m(t,1));     % Taylor rule
end


%----------------------------------------------------------------
% Wage and real interest rate
%----------------------------------------------------------------

% Labor-leisure FOC
vAggregateWage		                 = cchi ./ vAggregateMarginalUtility;

% Fisher equation
vAggregateR_real		             = vAggregateR_nom(1:T) ./ vAggregateInflation(2:T+1);


%----------------------------------------------------------------
% Relative price of output
%----------------------------------------------------------------

vAggregatePriceIntermediate(1:T)	 = exp((pphiPrice / (ggamma - 1)) * (log(vAggregateInflation(1:T)) - ...
										bbeta * log(vAggregateInflation(2:T+1)))) * pSS;


%----------------------------------------------------------------
% Prices relevant for decision rules
%----------------------------------------------------------------

% Create vector versions of the usual objects (over the whole transition path)
vP                           = vAggregatePriceIntermediate;
vQ                           = vAggregatePriceCapital;
vWage                        = vAggregateWage;
vA                           = (vP .^ (1 / (1 - nnu))) * (nnu ^ (nnu / (1 - nnu)) -...
									nnu ^ (1 / (1 - nnu)));

vSDF			             = vAggregateOutput;
vSDF(1:end) 	             = bbeta;
vSDF(1:T,1) 	             = bbeta * (vAggregateMarginalUtility(2:end) ./ vAggregateMarginalUtility(1:T));

% Last period: use no change (could also use steady state here)
vP(T+1,1)					 = vP(T,1);
vQ(T+1,1)					 = vQ(T,1);
vWage(T+1,1)			     = vWage(T,1);
vA(T+1,1)					 = vA(T,1);
vSDF(T+1,1)				     = vSDF(T,1);
vAggregateInflation(T+1,1)	 = vAggregateInflation(T,1);
