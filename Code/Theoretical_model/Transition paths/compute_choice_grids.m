% Declare global variables
global nChoices mCapitalGrid mDebtGrid mChoicesGrid capitalMin capitalMax vCapitalGrid ...
    debtMin debtMax vDebtGrid

% Nonlinear grid for capital
capitalMin               = capitalMinMult * min(vCapitalUnconstrained);
capitalMax               = capitalMaxMult * max(vCapitalUnconstrained);
vCapitalRaw			     = linspace(0,1,nCapital)';
vCapitalGrid		     = capitalMin + (capitalMax - capitalMin) * (vCapitalRaw .^ (1 / capitalPower));

% Nonlinear grid for debt (putting more points for high levels of debt)
debtMin                  = min(b0,debtMinMult * min(vDebtUnconstrained));
debtMax				     = debtMaxMult * max(vCapitalUnconstrained);
vDebtRaw			     = linspace(1,0,nDebt)';
vDebtGrid			     = debtMax - (debtMax - debtMin) * (vDebtRaw .^ (1 / debtPower));

% Compute tensor product grid of choices (useful for computations later on)
nChoices                 = nCapital * nDebt;
[mCapitalGrid,mDebtGrid] = ndgrid(vCapitalGrid,vDebtGrid);
mChoicesGrid             = [mCapitalGrid(:) mDebtGrid(:)];
