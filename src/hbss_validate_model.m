function HB = hbss_validate_model(HB)
% HBSS_VALIDATE_MODEL  Validate and normalize HB model fields.

req = {'A','Be','C','De'};
for k=1:numel(req)
    assert(isfield(HB,req{k}), 'HB.%s is required.', req{k});
end
if ~isfield(HB,'Ts'), HB.Ts = 0; end

A  = HB.A;  Be = HB.Be;  C  = HB.C;  De = HB.De;

% Basic sizes
nx = size(A,1);
assert(size(A,2)==nx, 'A must be square.');
assert(size(Be,1)==nx, 'Be must have nx rows.');
ny = size(C,1);
assert(size(C,2)==nx, 'C must have nx columns.');
assert(size(De,1)==ny, 'De must have ny rows.');
assert(size(De,2)==size(Be,2), 'De must have same number of columns as Be.');

% Time
if HB.Ts==0
    HB.ssType = 'c';
    HB.fs = 0;
else
    assert(HB.Ts>0, 'Ts must be >0 for discrete models.');
    HB.ssType = 'd';
    HB.fs = 1/HB.Ts;
end

% Harmonics
assert(isfield(HB,'Harmonics') && HB.Harmonics(1)==0, 'HB.Harmonics must start with 0.');

% NL fields
assert(isfield(HB,'NL') && isfield(HB.NL,'fNL'), 'HB.NL.fNL is required.');
nNL = numel(HB.NL.fNL);

% Forcing (numeric matrix or function handle)
assert(isfield(HB,'F'), 'HB.F is required.');
nu = size(Be,2);
nForc = nu - nNL;    HB.nForc = nForc;
assert(nForc >= 1, 'Extended input size mismatch: nu=%d, nNL=%d -> nForc=%d.', nu, nNL, nForc);
if isnumeric(HB.F)
    assert(size(HB.F,2)==1, 'HB.F must have 1 column (fundamental harmonic component).');
    assert(size(HB.F,1)==nForc, 'HB.F must have nForc=%d rows.', nForc);
elseif isa(HB.F,'function_handle')
    % check
    test = HB.F(1.0);
    assert(isvector(test) && numel(test)==nForc, ...
        'If HB.F is a function handle, it must return nForc x 1 (fundamental amplitude).');
    assert(any(HB.Harmonics==1), 'HB.Harmonics must include harmonic 1.');
else
    error('HB.F must be numeric or a function handle.');
end

% Defaults
if ~isfield(HB,'nSamples'), HB.nSamples = 512; end
if ~isfield(HB,'solver'), HB.solver = struct(); end

if ~isfield(HB.solver,'name'), HB.solver.name = 'newton'; end
if ~isfield(HB.solver,'opts')
    switch HB.solver.name
        case 'fsolve'
            HB.solver.opts = optimoptions('fsolve','Display','off',...
                'TolFun',1e-10,'TolX',1e-10,'MaxIter',50);
            % No diplay
            HB.solver.opts = optimoptions(HB.solver.opts,'Display','off');
        case 'newton'
            HB.solver.opts.tol=1e-8;            % Tolerance
            HB.solver.opts.maxIter=50;          % Maximum number of iterations
            HB.solver.opts.fdStep=1e-6;         % Finite difference step
    end
end


if ~isfield(HB,'scaleY'), HB.scaleY = 'auto'; end
if ~isfield(HB,'cont'), HB.cont = struct(); end
if ~isfield(HB.cont,'hmax'), HB.cont.hmax = 1; end
if ~isfield(HB.cont,'hmin'), HB.cont.hmin = 1e-5; end
if ~isfield(HB.cont,'maxIterCont'), HB.cont.maxIterCont = 10; end

if ~isfield(HB,'plot'), HB.plot = struct(); end
if ~isfield(HB.plot,'enable'), HB.plot.enable = true; end
if ~isfield(HB.plot,'chOut'), HB.plot.chOut = 1; end
if ~isfield(HB.plot,'updateEvery'), HB.plot.updateEvery = 10; end

% Store dimensions
HB.nx = nx; HB.ny = ny; HB.nu = nu; HB.nNL = nNL; HB.nForc = nForc;
end
