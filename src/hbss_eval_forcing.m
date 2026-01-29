function Fh = hbss_eval_forcing(HB, Omega)
% HBSS_EVAL_FORCING Returns forcing Fourier coefficients at frequency Omega

Harmonics = HB.Harmonics(:).';
nH = numel(Harmonics) - 1;
nForc = HB.nForc;

if isnumeric(HB.F)
    Fh = zeros(nForc, nH+1);
    ih1 = find(Harmonics==1,1);
    Fh(:,ih1) = HB.F;
    return;
end

if isa(HB.F,'function_handle')
    amp1 = HB.F(Omega);
    amp1 = amp1(:);
    assert(numel(amp1)==nForc, 'HB.F(Omega) must return nForc x 1.');

    Fh = zeros(nForc, nH+1);
    ih1 = find(Harmonics==1,1);
    assert(~isempty(ih1), 'HB.Harmonics must include harmonic 1.');
    Fh(:,ih1) = amp1;
    return;
end

error('HB.F must be numeric or a function handle.');
end
