function HB = hbss_validate_model(HB)
% HBSS_VALIDATE_MODEL  Validate HB model fields.
%
% -------------------------------------------------------------------------
% Author: Dario Anastasio (Politecnico di Torino)
% Part of: HBSS_Cont toolbox
% -------------------------------------------------------------------------

%% MINIMUM REQUIREMENTS AND STATE DIMENSIONS
req = {'A','Be','C','De'};
for k = 1:numel(req)
    assert(isfield(HB,req{k}), 'HB.%s is required.', req{k});
end

A  = HB.A;  Be = HB.Be;  C  = HB.C;  De = HB.De;
nx = size(A,1);
ny = size(C,1);
nu = size(Be,2);

assert(size(A,2)==nx, 'A must be square.');
assert(size(Be,1)==nx, 'Be must have nx rows.');
assert(size(C,2)==nx, 'C must have nx columns.');
assert(size(De,1)==ny, 'De must have ny rows.');
assert(size(De,2)==nu, 'De must have same number of columns as Be.');

%% TIME AND SAMPLING
if ~isfield(HB,'Ts'), HB.Ts = 0; end

if HB.Ts == 0
    HB.ssType = 'c';
    HB.fs = 0;
else
    assert(HB.Ts > 0, 'Ts must be >0 for discrete models.');
    HB.ssType = 'd';
    HB.fs = 1/HB.Ts;

    % Build continuous surrogate model for stability analisys
    sysd = struct();
    sysd.A  = HB.A;
    sysd.Be = HB.Be;
    sysd.C  = HB.C;
    sysd.De = HB.De;

    sysc = hbss_d2c_zoh(sysd, HB.Ts);
    HB.sys_c = struct();
    HB.sys_c.A  = sysc.A;
    HB.sys_c.Be = sysc.Be;
    HB.sys_c.C  = sysc.C;
    HB.sys_c.De = sysc.De;
end

if ~isfield(HB,'nSamples'), HB.nSamples = 512; end

%% HARMONICS AND NONLINEARITIES
assert(isfield(HB,'Harmonics') && HB.Harmonics(1)==0, 'HB.Harmonics must start with 0.');

assert(isfield(HB,'NL') && isfield(HB.NL,'fNL'), 'HB.NL.fNL is required.');
nNL = numel(HB.NL.fNL);

%% FORCING TERM (Numeric Matrix or Function Handle)
assert(isfield(HB,'F'), 'HB.F is required.');
nForc = nu - nNL;
assert(nForc >= 1, 'Extended input size mismatch: nu=%d, nNL=%d -> nForc=%d.', nu, nNL, nForc);

if isnumeric(HB.F)
    assert(size(HB.F,2)==1, 'HB.F must have 1 column (fundamental harmonic component).');
    assert(size(HB.F,1)==nForc, 'HB.F must have nForc=%d rows.', nForc);
elseif isa(HB.F,'function_handle')
    test = HB.F(1.0);
    assert(isvector(test) && numel(test)==nForc, ...
        'If HB.F is a function handle, it must return nForc x 1 (fundamental amplitude).');
    assert(any(HB.Harmonics==1), 'HB.Harmonics must include harmonic 1.');
else
    error('HB.F must be numeric or a function handle.');
end

%% SOLVER OPTIONS (Newton / fsolve)
if ~isfield(HB,'solver'), HB.solver = struct(); end
if ~isfield(HB.solver,'name'), HB.solver.name = 'newton'; end

if ~isfield(HB.solver,'opts')
    switch HB.solver.name
        case 'fsolve'
            HB.solver.opts = optimoptions('fsolve','Display','off',...
                'TolFun',1e-10,'TolX',1e-10,'MaxIter',50);
        case 'newton'
            HB.solver.opts.tol = 1e-8;
            HB.solver.opts.maxIter = 50;
            HB.solver.opts.fdStep = 1e-6;
    end
end

%% CONTINUATION AND SCALING OPTIONS
if ~isfield(HB,'scaleY'), HB.scaleY = 'auto'; end
if ~isfield(HB,'cont'), HB.cont = struct(); end

if ~isfield(HB.cont,'hmax'),         HB.cont.hmax = 1; end
if ~isfield(HB.cont,'hmin'),         HB.cont.hmin = 1e-5; end
if ~isfield(HB.cont,'maxIterCont'),  HB.cont.maxIterCont = 10; end

%% BIFURCATION OPTIONS (Post-Processing Analysis)
if ~isfield(HB,'bifOpts'), HB.bifOpts = struct(); end

if ~isfield(HB.bifOpts,'tolRad'),    HB.bifOpts.tolRad = 0.1;    end
if ~isfield(HB.bifOpts,'nPersist'),  HB.bifOpts.nPersist = 3;    end
if ~isfield(HB.bifOpts,'tolCross'),  HB.bifOpts.tolCross = 1e-4; end
if ~isfield(HB.bifOpts,'tolConj'),   HB.bifOpts.tolConj = 1e-3;  end

%% REAL-TIME PLOTTING OPTIONS
if ~isfield(HB,'plot'), HB.plot = struct(); end

if ~isfield(HB.plot,'enable'),       HB.plot.enable = true; end
if ~isfield(HB.plot,'chOut'),        HB.plot.chOut = 1; end
if ~isfield(HB.plot,'updateEvery'),  HB.plot.updateEvery = 10; end

%% STORE DIMENSIONS AND RETURN
HB.nx = nx; 
HB.ny = ny; 
HB.nu = nu; 
HB.nNL = nNL; 
HB.nForc = nForc;

end