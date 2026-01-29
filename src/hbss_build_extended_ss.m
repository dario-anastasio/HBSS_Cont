function [A, Be, C, De, model] = hbss_build_extended_ss(M, Cv, K, Lf, Lnl, mu)
% HBSS_BUILD_EXTENDED_SS  Build extended state-space matrices from an analytical model
%
% Builds the continuous-time extended state-space model:
%   xdot = A x + Be u_e
%   y    = C x + De u_e
%
% with states x = [y; ydot], outputs y, and extended inputs:
%   u_e = [ u ; g1 ; ... ; gJ ]
%
% Assumptions / conventions:
% - The external forcing input u is a scalar and it is the FIRST channel of u_e.
% - Nonlinear channels are scalar inputs g_j (defined by the user as g_j = -f_j(y, ydot)).
% - mu(j) scales the j-th nonlinear channel inside Be.
% - Outputs are always displacements y.
%
% Inputs:
%   M       (n×n)  mass matrix
%   Cv      (n×n)  viscous damping matrix
%   K       (n×n)  stiffness matrix
%   Lf      (n×1)  force distribution vector for the scalar input u
%                  (e.g., Lf = [1; 0; 0] means forcing applied at DOF 1 out of 3)
%   Lnl     (n×J)  position of the nonlinear basis functions
%   mu      (1×J) or (J×1) nonlinear gains multiplying each column of Lnl
%
%
% Output:
%   extended state-space matrices: A, Be, C, De
%   model struct 
%
% -------------------------------------------------------------------------
% Author: Dario Anastasio (Politecnico di Torino)
% Part of: HBSS_Cont toolbox
% -------------------------------------------------------------------------

n = size(M,1);
assert(ismatrix(M) && size(M,2)==n, 'M must be square.');
assert(ismatrix(Cv) && all(size(Cv)==[n n]), 'Cq must be n×n.');
assert(ismatrix(K)  && all(size(K)==[n n]),  'K must be n×n.');
assert(isvector(Lf) && numel(Lf)==n, 'Lf must be n×1.');
Lf = Lf(:);

if isempty(Lnl)
    Lnl = zeros(n,0);
end
assert(size(Lnl,1)==n, 'Lnl must have n rows.');

J = size(Lnl,2);
if isempty(mu)
    mu = zeros(J,1);
end
mu = mu(:);
assert(numel(mu)==J, 'mu must have length equal to number of nonlinear channels (size(Lnl,2)).');

% -------------------- Continuous-time matrices --------------------
Minv = M \ eye(n);

A = [zeros(n) eye(n);
     -Minv*K   -Minv*Cv];

% Build q-level input directions: [ external ; nonlinear channels ]
% External input is scalar: direction inForce
% Nonlinear channels: mu(j)*Lnl(:,j)
Bq = zeros(n, 1+J);
Bq(:,1) = Lf;
for j = 1:J
    Bq(:,1+j) = mu(j) * Lnl(:,j);
end

Be = [zeros(n, 1+J);
      Minv * Bq];

% Outputs are displacements y = q
C  = [eye(n) zeros(n)];
De = zeros(n, 1+J);

% -------------------- Structure --------------------
model = struct();
model.A = A;
model.Be = Be;
model.C = C;
model.De = De;

model.nDof = n;
model.nx = size(A,1);
model.ny = size(C,1);
model.nu = size(Be,2);

model.nForc = 1;
model.nNL   = J;
end
