function He = hbss_extended_FRF(HB, h, Omega)
% HBSS_EXTENDED-FRF extended FRF matrix He
%
% -------------------------------------------------------------------------
% Author: Dario Anastasio (Politecnico di Torino)
% Part of: HBSS_Cont tooolbox
% -------------------------------------------------------------------------

A  = HB.A;
Be = HB.Be;
C  = HB.C;
De = HB.De;
nx = size(A,1);

if HB.ssType == 'c'
    Tmat = 1i*h*Omega*eye(nx) - A;
elseif HB.ssType == 'd'
    z_var = exp(1i*h*Omega*HB.Ts);
    Tmat = z_var*eye(nx) - A;
else
    error('Unknown HB.ssType: must be ''c'' or ''d''.');
end

[Lt,Ut] = lu(Tmat);
He = De + ((C/Ut) * (Lt\Be));
end