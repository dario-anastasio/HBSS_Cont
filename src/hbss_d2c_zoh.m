function sysc = hbss_d2c_zoh(sysd, Ts)
% HBSS_D2C_ZOH  Discrete-to-continuous conversion with zoh
%
%   sysc = hbss_d2c_zoh(sysd, Ts)
%
% Inputs:
%   sysd : struct with fields
%          .A  (nx x nx)  discrete-time state matrix
%          .Be (nx x nu)  discrete-time input matrix (extended inputs)
%          .C  (ny x nx)  output matrix
%          .De (ny x nu)  feedthrough matrix
%   Ts  : sample time [s], scalar > 0
%
% Output:
%   sysc : struct with fields .A, .Be, .C, .De (continuous-time equivalent),
%          plus .Ts=0
%
% -------------------------------------------------------------------------
% Author: Dario Anastasio (Politecnico di Torino)
% Part of: HBSS_Cont tooolbox
% -------------------------------------------------------------------------

assert(isstruct(sysd) && all(isfield(sysd,{'A','Be','C','De'})), ...
    'sysd must be a struct with fields A, Be, C, De.');
assert(isscalar(Ts) && Ts > 0, 'Ts must be a positive scalar.');

Ad = sysd.A; Bd = sysd.Be; Cc = sysd.C; Dc = sysd.De;

[nx1,nx2] = size(Ad);
assert(nx1==nx2, 'A must be square.');

[nxB,nu] = size(Bd);
assert(nxB==nx1, 'Be must have nx rows.');

% Build matrix for ZOH
Md = [Ad, Bd;
      zeros(nu,nx1), eye(nu)];

% Conditioning check
rc = rcond(Md);
if rc < 1e-10
    warning(['hbss_d2c_zoh: ill conditioned matrix ', ...
             '(rcond = %.2e). ', ...
             'Discrete-to-continuous conversion via logm may be inaccurate.'], rc);
end

% Compute matrix logarithm
[Lm, logmWarn] = logm(Md);
if logmWarn ~= 0
    warning(['hbss_d2c_zoh: ', ...
             'logm returned a warning flag (%d). ', ...
             'Continuous equivalent may be unreliable.'], logmWarn);
end

L   = Lm / Ts;
Ac  = L(1:nx1, 1:nx1);
Bec = L(1:nx1, nx1+1:nx1+nu);

% Output
sysc = struct();
sysc.A  = Ac;
sysc.Be = Bec;
sysc.C  = Cc;
sysc.De = Dc;
sysc.Ts = 0;
end
