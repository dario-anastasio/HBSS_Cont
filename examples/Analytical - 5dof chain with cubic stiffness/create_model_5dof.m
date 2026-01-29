%% --- Linear parameters  ---

%  Number of DOF
nDof = 5;

% Mass matrix
m   = 1;
M = diag(ones(nDof,1)*m);

% Stiffness matrix 
k = 1e4;
K = [2*k   -k   0    0    0;
     -k  2*k  -k    0    0;
      0   -k  2*k  -k    0;
      0    0   -k  2*k  -k;
      0    0    0   -k  2*k];

% Damping matrix
c  = 5;
Cv = K*c/k;

% External force location vector
Lf = zeros(nDof,1); Lf(1)=1;

%% --- Nonlinear parameters ---

J = 5;          % Number of nonlinear terms

% Nonlinear coefficients
mu = 1e7*ones(1,5);

% Nonlinear basis functions
fNL = cell(J,1);
fNL{1} = @(y,yd) (y(1,:)-y(2,:)).^3;  % g1 
fNL{2} = @(y,yd) (y(2,:)-y(3,:)).^3;  % g2
fNL{3} = @(y,yd) (y(3,:)-y(4,:)).^3;  % g3
fNL{4} = @(y,yd) (y(4,:)-y(5,:)).^3;  % g4
fNL{5} = @(y,yd) (y(5,:)).^3;         % g5

% Location vectors
Lnl = zeros(nDof, J);
Lnl(:,1) = [ 1; -1;  0;  0;  0];   % 1-2
Lnl(:,2) = [ 0;  1; -1;  0;  0];   % 2-3
Lnl(:,3) = [ 0;  0;  1; -1;  0];   % 3-4
Lnl(:,4) = [ 0;  0;  0;  1; -1];   % 4-5
Lnl(:,5) = [ 0;  0;  0;  0;  1];   % 5-ground

% Analytic derivatives
dFNL    = cell(J, nDof);
dFNLvel = cell(J, nDof);

% Initialize all to zero
for inl = 1:J
    for j = 1:nDof
        dFNL{inl,j}    = @(y,yd) 0*y(1,:);
        dFNLvel{inl,j} = @(y,yd) 0*y(1,:);
    end
end

% g1 = (y1-y2)^3
dFNL{1,1} = @(y,yd)  3*(y(1,:)-y(2,:)).^2;
dFNL{1,2} = @(y,yd) -3*(y(1,:)-y(2,:)).^2;

% g2 = (y2-y3)^3
dFNL{2,2} = @(y,yd)  3*(y(2,:)-y(3,:)).^2;
dFNL{2,3} = @(y,yd) -3*(y(2,:)-y(3,:)).^2;

% g3 = (y3-y4)^3
dFNL{3,3} = @(y,yd)  3*(y(3,:)-y(4,:)).^2;
dFNL{3,4} = @(y,yd) -3*(y(3,:)-y(4,:)).^2;

% g4 = (y4-y5)^3
dFNL{4,4} = @(y,yd)  3*(y(4,:)-y(5,:)).^2;
dFNL{4,5} = @(y,yd) -3*(y(4,:)-y(5,:)).^2;

% g5 = (y5)^3
dFNL{5,5} = @(y,yd)  3*(y(5,:)).^2;

%% --- Build Extended State-Space ---
[A, Be, C, De, model] = hbss_build_extended_ss(M, Cv, K, Lf, Lnl, mu);

clear M Cv K m c k mu nDof Lf Lnl;
