function [rY, aux] = hbss_residual(Y_scaled, Omega, HB)
% HBSS_RESIDUAL  HB residual
%
% -------------------------------------------------------------------------
% Author: Dario Anastasio (Politecnico di Torino)
% Part of: HBSS_Cont toolbox
% -------------------------------------------------------------------------

Harmonics = HB.Harmonics(:).';
nH = numel(Harmonics)-1;
ny = HB.ny;
scaleY = HB.scaleY;
Y_phys = Y_scaled(:) / scaleY;
Yh_unknown = hbss_Y_to_Yh(Y_phys, ny, Harmonics);  

% AFT nonlinear forces
[Fnl, yt, fnl, ytd] = hbss_aft(Yh_unknown, Omega, HB);   
GNL = -Fnl(:,1:nH+1);

% Target from extended transfer
Yh_target = zeros(nH+1, ny);

A  = HB.A;  Be = HB.Be;  C = HB.C;  De = HB.De;

% Current force amplitude
Fh = hbss_eval_forcing(HB, Omega);

for ih = 1:(nH+1)
    h = Harmonics(ih);

    Fext = [Fh(:,ih); GNL(:,ih)];
    Tmat = 1i*h*Omega*eye(size(A)) - A;
   
    [Lt,Ut] = lu(Tmat);
    He = De + ((C/Ut)*(Lt\Be));

    y_h = He * Fext;          
    Yh_target(ih,:) = y_h.';   
end

% Target
Y_target = hbss_Yh_to_Y(Yh_target, ny, Harmonics);

% Scaled residual
rY = Y_scaled(:) - scaleY * Y_target(:);

aux = struct();
aux.Yh_unknown = Yh_unknown;
aux.Yh_target  = Yh_target;
aux.Fnl       = Fnl;
aux.FNL       = GNL;
aux.yt        = yt;
aux.ytd       = ytd;
aux.fnl_time  = fnl;
end
