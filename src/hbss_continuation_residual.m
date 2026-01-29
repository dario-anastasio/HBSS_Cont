function [r, aux] = hbss_continuation_residual(contPoint, contPointPred, contTangentPrev, HB)
% HBSS_CONTINUATION_RESIDUAL  residual for corrector step
%
% contPoint       = [Y_scaled; Omega] current unknown in the corrector
% contPointPred   = predictor point (previous contPoint + h * contTangent)
% contTangentPrev = tangent vector used for pseudo-arc-length constraint
%
% Pseudo-arc-length constraint:
%   contTangentPrev' * (contPoint - contPointPred) = 0
%
% -------------------------------------------------------------------------
% Author: Dario Anastasio (Politecnico di Torino)
% Part of: HBSS_Cont toolbox
% -------------------------------------------------------------------------

Y_scaled = contPoint(1:end-1);
Omega    = contPoint(end);

[rY, auxHB] = hbss_residual(Y_scaled, Omega, HB);

% Pseudo-arc-length constraint
rArc = contTangentPrev.' * (contPoint - contPointPred);

r = [rY; rArc];
aux = auxHB;
end
