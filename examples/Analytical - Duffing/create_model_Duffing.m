%% --- Linear Parameters ---

% Mass
m = 1; 

% Stiffness
k = 1; 

% Damping
cv = 2*0.02;

% Forcing position
Lf = 1; 

%% --- Nonlinear parameters ----

% Cubic spring
mu(1)      = 0.1;                  % Cubic coefficient
fNL{1}     = @(y, yd) y.^3;
dFNL{1}    = @(y, yd) 3*y.^2;
dFNLvel{1} = @(y, yd) 0;
Lnl(1)     = 1;

% % (optional) Quadratic spring
% mu(2)      = 0.1;                  % Quadratic coefficient
% fNL{2}     = @(y, yd) y.^2;
% dFNL{2}    = @(y, yd) 2*y;
% dFNLvel{2} = @(y, yd) 0;
% Lnl(2)     = 1;

%% --- Build Extended State-Space ---
[A, Be, C, De, model] = hbss_build_extended_ss(m, cv, k, Lf, Lnl, mu);

clear m c k mu Lf Lnl;