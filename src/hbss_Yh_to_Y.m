function Y = hbss_Yh_to_Y(U, ny, Harmonics)
% Convert complex coefficients to [re;im] stack
%
% -------------------------------------------------------------------------
% Author: Dario Anastasio (Politecnico di Torino)
% Part of: HBSS_Cont toolbox
% -------------------------------------------------------------------------


Harmonics = Harmonics(:).';
assert(Harmonics(1)==0, 'Harmonics must start with 0.');
nH = numel(Harmonics)-1;

assert(size(U,1)==nH+1 && size(U,2)==ny, 'U must be (nH+1)×ny.');

Y = zeros(ny*(2*nH+1),1);

% DC
Y(1:ny) = real(U(1,:)).';

k = ny;
for h = 1:nH
    Yc = real(U(h+1,:)).';
    Ys = -imag(U(h+1,:)).';
    Y(k+1:k+ny) = Yc; k = k+ny;
    Y(k+1:k+ny) = Ys; k = k+ny;
end
end
