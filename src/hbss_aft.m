function [Fnl, yt, fnl, ytd] = hbss_aft(Yh, Omega, HB)
% HBSS_AFT  Alternating Frequency-Time (AFT) for nonlinear forces.
%
% Inputs:
%   Yh     (nH+1)×ny complex, one-sided Fourier coefficients of outputs
%   Omega  scalar rad/s
%   HB     struct
%
% Outputs:
%   Fnl    nNL×(nH+1) complex Fourier coefficients (one-sided)
%   yt     ny×Nstep real time histories over one period
%   fnl    nNL×Nstep real nonlinear forces in time
%   ytd    ny×Nstep real time derivative of yt
%
% -------------------------------------------------------------------------
% Author: Dario Anastasio (Politecnico di Torino)
% Part of: HBSS_Cont toolbox
% -------------------------------------------------------------------------

assert(isscalar(Omega) && Omega > 0, 'Omega must be a positive scalar.');

Nstep = HB.nSamples;
fNL   = HB.NL.fNL;
nNL   = numel(fNL);

[nH1, ny] = size(Yh);
nH = nH1 - 1;

% IFFT to time domain (periodic, one period)
yt = zeros(ny, Nstep);
for iy = 1:ny
    Yr = zeros(Nstep, 1);
    Yr(1) = Yh(1, iy);
    Yr(2:nH+1) = Yh(2:nH+1, iy) / 2;
    Yr(Nstep-nH+1:Nstep) = conj(Yh(nH+1:-1:2, iy)) / 2;
    yt(iy, :) = ifft(Nstep * Yr).';
end
yt = real(yt);

% Time step (one period)
T  = 2*pi / Omega;
dt = T / Nstep;

% Periodic finite-difference derivative
ytd = zeros(ny, Nstep);
ytd(:, 2:Nstep-1) = (yt(:, 3:Nstep) - yt(:, 1:Nstep-2)) / (2*dt);
ytd(:, 1)         = (yt(:, 2) - yt(:, end)) / (2*dt);
ytd(:, end)       = (yt(:, 1) - yt(:, end-1)) / (2*dt);

% Nonlinear forces in time
fnl = zeros(nNL, Nstep);
for inl = 1:nNL
    fnl(inl, :) = fNL{inl}(yt, ytd);
end
fnl = real(fnl);

% FFT back to one-sided coefficients
Fnl = zeros(nNL, nH+1);
for inl = 1:nNL
    Fnlt = fft(fnl(inl, :)) / Nstep;
    Fnlt = Fnlt(1:nH+1);
    Fnlt(2:end) = 2 * Fnlt(2:end);
    Fnl(inl, :) = Fnlt;
end
end
