function U = hbss_Y_to_Yh(Y, ny, Harmonics)
% Convert [re;im] stack to complex coefficients
%
% -------------------------------------------------------------------------
% Author: Dario Anastasio (Politecnico di Torino)
% Part of: HBSS_Cont toolbox
% -------------------------------------------------------------------------

Harmonics = Harmonics(:).';
assert(Harmonics(1)==0, 'Harmonics must start with 0.');
nH = numel(Harmonics)-1;

Y = Y(:);
assert(numel(Y) == ny*(2*nH+1), 'Size mismatch: Y must be ny*(2*nH+1).');

U = zeros(nH+1, ny);

% DC
U(1,:) = Y(1:ny).';

k = ny;
for h = 1:nH
    Yc = Y(k+1:k+ny); k = k+ny;
    Ys = Y(k+1:k+ny); k = k+ny;
    U(h+1,:) = (Yc - 1i*Ys).';
end
end
