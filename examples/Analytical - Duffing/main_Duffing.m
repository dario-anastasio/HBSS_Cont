%% examples/main_Duffing.m
% HBSS_Cont - Example: Duffing oscillator
% State space model is created in create_model_Duffing script. 
% Quadratic stiffness can be added by de-commenting lines 24-28. 

clearvars; close all; clc;

%% Setup path (HBSS_Cont)
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
Harmonics = 0:7;           % Harmonic indices [0 1 ... H]
nSamples  = 1024;          % Time samples for AFT

% =========================
% Frequency range
% =========================
fLim = [0.1 5]*1e-1;      % Frequency range [Hz]

% =========================
% Forcing
% =========================
F0 = 0.2;                 % Amplitude of harmonic forcing

% =========================
% Scaling (important!)
% =========================
% Scaling factor for harmonic coefficients (improves numerical conditioning)
scaleY = 'auto';           % Automatic selection

% =========================
% Pseudo-Arc-Length Continuation
% =========================
% hmax: maximum allowed change per continuation step
% hmin: minimum step size before stopping
hmax = 0.1;               
hmin = 1e-3;              

% =========================
% Plotting
% =========================
chOut = 1;                 % Output channel to plot (index)
updateEvery = 20;          % Update real-time plot every x iterations

% =========================
% Nonlinear solver
% =========================
% These parameters control the local Newton solver used
% inside each continuation corrector step
solverTol = 1e-10;         % Solver tolerance
maxIter = 20;               % Max Newton iterations

%% ---- State-space from analytical model ----

% Continuous-time model from theoretical M,C,K matrices
create_model_Duffing;

% Number of outputs
ny = size(C,1);

%% ---- Build HB struct ---

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

% First 5 harmonic contributions
chPlot = 1;       % output channel
hPlot  = 1:5;     % maximum
hbss_plot_nfrc(out, chPlot, hPlot);

% Maximum amplitude
chPlot = 1;       % output channel
hPlot  = "max";   % maximum
hbss_plot_nfrc(out, chPlot, hPlot);