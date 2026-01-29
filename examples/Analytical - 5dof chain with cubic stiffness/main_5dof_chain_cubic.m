%% examples/main_5dof_chain_cubic.m
% HBSS_Cont - Example: 5DOF chain with cubic nonlinearities

% Numerical example taken from:
% Karaagaçlı T., Özgüven H. N., MSSP 2021
% https://doi.org/10.1016/j.ymssp.2020.107023
 
clearvars; close all; clc;

%% ---- Setup path (required!) ----
thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);          % .../examples/...
examplesDir = fileparts(thisDir);        % .../examples
rootDir     = fileparts(examplesDir);    % .../HBSS_Cont
srcDir = fullfile(rootDir, 'src');
addpath(srcDir);

%% ---- User defined parameters ----

% =========================
% Harmonic Balance
% =========================
Harmonics = 0:5;           % Harmonic indices [0 1 ... H]
nSamples  = 512;           % Time samples for AFT

% =========================
% Frequency range
% =========================
fLim = [5 35];             % Frequency range [Hz]

% =========================
% Forcing
% =========================
F0 = 50;                 % Fundamental amplitude of harmonic forcing

% =========================
% Scaling (important!)
% =========================
% Scaling factor for harmonic coefficients (improves numerical conditioning)
scaleY = 'auto';            % Automatic selection
% scaleY = 1e3;             

% =========================
% Continuation steps
% =========================
% hmax: maximum allowed change per continuation step
% hmin: minimum step size before stopping
hmax = 2;                  
hmin = 1e-2;               

% =========================
% Plotting
% =========================
chOut = 1;                 % Output channel to plot (index)
updateEvery = 10;          % Update real-time plot every x iterations

% =========================
% Nonlinear solver
% =========================
% These parameters control the local Newton solver used
% inside each continuation corrector step
solverTol = 1e-10;         % Solver tolerance
maxIter = 20;               % Max Newton iterations

%% ---- State-space from analytical model ----

% Continuous-time model from theoretical M,C,K matrices
create_model_5dof;

% Number of outputs
ny = size(C,1);

%% ---- Build HB struct ----

HB = struct();
HB.A  = A; HB.Be = Be;
HB.C  = C; HB.De = De;
HB.Harmonics = Harmonics;
HB.fLim      = fLim;
HB.scaleY    = scaleY;
HB.nSamples  = nSamples;
HB.F = F0;
HB.cont = struct('hmax',hmax,'hmin',hmin);
HB.plot = struct('enable',true,'chOut',chOut,'updateEvery',updateEvery);

HB.solver = struct();
HB.solver.name = 'newton';
HB.solver.opts = struct('tol',solverTol,'maxIter',maxIter);

HB.NL = struct();
HB.NL.fNL  = fNL;
try
    HB.NL.dFNL = dFNL;    % <- analytic Jacobian of NL (recommended)
catch
    % If no dFNL is defined -> numerical Jacobian
end
try
    HB.NL.dFNLvel = dFNLvel;    % <- analytic Jacobian of NL (recommended)
catch
    % If no dFNL is defined -> numerical Jacobian
end

%% ---- Run ----
out = hbss_continuation(HB);

%% ---- Post plot ----

% Fundamental harmonic response of all output channels
chPlot = 1:5;    % all outputs     
hPlot  = 1;      % fundamental
hbss_plot_nfrc(out, chPlot, hPlot);
