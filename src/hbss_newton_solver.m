function [z, info] = hbss_newton_solver(fun, z0, opts)
% HBSS_NEWTON_SOLVER Newton solver
%
% INPUTS:
%   fun   : function handle returning residual vector
%   z0    : initial guess
%   opts  : struct with optional fields
%       .tol        : convergence tolerance on norm(residual)
%       .maxIter    : maximum number of Newton iterations
%       .fdStep     : finite-difference step (relative)
%
% OUTPUTS:
%   z     : converged solution
%   info  : struct with convergence information
%
% -------------------------------------------------------------------------
% Author: Dario Anastasio (Politecnico di Torino)
% Part of: HBSS_Cont toolbox
% -------------------------------------------------------------------------

% Defaults
if ~isfield(opts,'tol'),      opts.tol = 1e-8;        end
if ~isfield(opts,'maxIter'),  opts.maxIter = 20;      end
if ~isfield(opts,'fdStep'),   opts.fdStep = 1e-6;     end

z = z0(:);
n = numel(z);

for it = 1:opts.maxIter

    R = fun(z);
    nR = norm(R);

    if nR < opts.tol
        info.converged = true;
        info.iter = it;
        info.resNorm = nR;
        return;
    end

    % --- Finite-difference Jacobian ---
    J = zeros(numel(R), n);
    h = opts.fdStep * max(1, abs(z));

    for j = 1:n
        zj = z;
        zj(j) = zj(j) + h(j);
        Rj = fun(zj);
        J(:,j) = (Rj - R) / h(j);
    end

    % --- Newton step ---
    dz = - J \ R;

    % --- Basic damping (optional) ---
    alpha = 1.0;
    z = z + alpha * dz;
end

info.converged = false;
info.iter = opts.maxIter;
info.resNorm = nR;
warning('hbss_newton_solver:NoConvergence', ...
        'Newton solver did not converge.');
end
