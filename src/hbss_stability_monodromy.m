function stab = hbss_stability_monodromy(y, HB)
% HBSS_STABILITY_MONODROMY  Floquet multipliers via monodromy
%
% Inputs:
%   y  = [Y_scaled; Omega], 
%
% Output struct 'stab':
%   stab.sigma   : Floquet multipliers
%
% -------------------------------------------------------------------------
% Author: Dario Anastasio (Politecnico di Torino)
% Part of: HBSS_Cont toolbox
% -------------------------------------------------------------------------

Y_scaled = y(1:end-1);
Omega    = abs(y(end));
Harmonics = HB.Harmonics(:).';

A  = HB.A;  Be = HB.Be;  C = HB.C;  De = HB.De;
nx = HB.nx; ny = HB.ny; nu = HB.nu;
nNL = HB.nNL; nForc = HB.nForc;

% Period and discretization
T = 2*pi/Omega;
N = HB.nSamples;
dt = T / N;

Y_phys = Y_scaled(:) / HB.scaleY;
Yh = hbss_Y_to_Yh(Y_phys, ny, Harmonics);   % (nH+1)×ny

% AFT
[~, yt, ~, ytd] = hbss_aft(Yh, Omega, HB);

% Matrices to store partial derivatives d(fnl)/dy and d(fnl)/dyd
DFy  = zeros(nNL, ny, N);
DFyd = zeros(nNL, ny, N);

% Check for analytical Jacobians
hasAnalyticY  = isfield(HB,'NL') && isfield(HB.NL,'dFNL') && iscell(HB.NL.dFNL) ...
                && all(size(HB.NL.dFNL) == [nNL, ny]);
hasAnalyticYd = isfield(HB,'NL') && isfield(HB.NL,'dFNLvel') && iscell(HB.NL.dFNLvel) ...
                && all(size(HB.NL.dFNLvel) == [nNL, ny]);
if ~hasAnalyticY
    epsFD = 1e-7;
end
if ~hasAnalyticYd
    epsFD = 1e-7;
end

% Go
for k = 1:N
    xk  = yt(:,k);
    xdk = ytd(:,k);
    for inl = 1:nNL
        for j = 1:ny
            if hasAnalyticY
                DFy(inl,j,k) = HB.NL.dFNL{inl,j}(xk, xdk);
            else
                xp = xk; xp(j) = xp(j) + epsFD;
                DFy(inl,j,k) = (HB.NL.fNL{inl}(xp, xdk) - HB.NL.fNL{inl}(xk, xdk)) / epsFD;
            end
            
            if hasAnalyticYd
                DFyd(inl,j,k) = HB.NL.dFNLvel{inl,j}(xk, xdk);
            else
                xdp = xdk; xdp(j) = xdp(j) + epsFD;
                DFyd(inl,j,k) = (HB.NL.fNL{inl}(xk, xdp) - HB.NL.fNL{inl}(xk, xdk)) / epsFD;
            end
        end
    end
end

% Nonlinear input is g = -fnl, so dg/dy = -DFy and dg/dyd = -DFyd
dGdY  = -DFy;
dGdYd = -DFyd;

% Build Jacobian Jk over time 
J = zeros(nx, nx, N);
for k = 1:N
    DUeY  = zeros(nu, ny);
    DUeYd = zeros(nu, ny);
    DUeY(nForc+1:end, :)  = dGdY(:,:,k);
    DUeYd(nForc+1:end, :) = dGdYd(:,:,k);
        
    M = (eye(nu) - DUeY * De - DUeYd * (C * Be));
    
    J(:,:,k) = A + Be * (M \ (DUeY * C + DUeYd * (C * A)));
end


% --- Monodromy matrix ---
H = eye(nx);

% % Matrix exponential (slow)
% for k = 1:N
%     H = expm(J(:,:,k) * dt) * H;
% end

% Second order approximation (faster)
I = eye(size(J,1));
for k = 1:N
    Jdt = J(:,:,k) * dt;
    H = (I + Jdt + 0.5 * (Jdt * Jdt)) * H;
end

% Floquet multipliers
sigma = eig(H);
sigma = sigma(:);

stab = struct();
stab.sigma = sigma;

end