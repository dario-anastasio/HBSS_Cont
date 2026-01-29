function [ymax, ymin] = hbss_maxmin_from_Yh(Yh, Harmonics, Omega, nSamples)
% Time-domain max/min for each output channel from harmonics
%
% Yh: (nH+1)×ny complex one-sided coefficients (DC + positive harmonics)
% Harmonics: [0 1 ... nH]
% Omega: rad/s
% nSamples: number of time samples over one period
%
% -------------------------------------------------------------------------
% Author: Dario Anastasio (Politecnico di Torino)
% Part of: HBSS_Cont toolbox
% -------------------------------------------------------------------------

Harmonics = Harmonics(:).';
nH = numel(Harmonics)-1;
ny = size(Yh,2);

T = 2*pi/abs(Omega);
t = linspace(0, T, nSamples);

yt = zeros(ny, nSamples);
for h = 0:nH
    Yh1 = Yh(h+1,:).'; % ny×1
    if h==0
        yt = yt + real(Yh1) * ones(1,nSamples);
    else
        yt = yt + real(Yh1)*cos(h*abs(Omega)*t) - imag(Yh1)*sin(h*abs(Omega)*t);
    end
end

ymax = max(yt,[],2);
ymin = min(yt,[],2);
end
