%----------------------------------------------------------------
% Solve unconstrained problem
%----------------------------------------------------------------

steady_state_unconstrained;

%----------------------------------------------------------------
% Compute default threshold
%----------------------------------------------------------------

% Use results from unconstrained decisions to make capital and debt choice grids
compute_choice_grids;

% Comptute cutoff
steady_state_default_threshold;

%----------------------------------------------------------------
% Compute value function and decision rules
%----------------------------------------------------------------

% Use results from unconstrained decisions and default thresholds to compute grids for cash on hand
compute_cash_grid;

% Compute full value function using value function iteration
steady_state_value_function;

% Compute decisions using continuous optimizer given the converged value function
steady_state_continuous_decisions;
mCapitalPrime       = mCapitalPrimeContinuous;
mDebtPrime          = mDebtPrimeContinuous;
mDividends          = mDividendsContinuous;
vDebtPriceOptimal   = mDebtPriceContinuous;


%----------------------------------------------------------------
% Compute stationary distribution
%----------------------------------------------------------------

compute_distribution_grids;

steady_state_distribution;
