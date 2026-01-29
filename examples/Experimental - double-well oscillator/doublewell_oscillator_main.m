%% examples/doublewell_oscillator_main.m
% HBSS_Cont - Example: double-well oscillator
% State-space model obtained experimentally with NSI method

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
nSamples  = 512;           % Time samples for AFT

% =========================
% Frequency range
% =========================
fLim = [1 16];             % Frequency range [Hz]

% =========================
% Forcing
% =========================
% Forcing amplitude as a function of omega from experimental testing
F0 = @(omega) 0.0014*omega.^2 - 0.0127*omega + 0.4;

% =========================
% Scaling (important!)
% =========================
% Scaling factor for harmonic coefficients (improves numerical conditioning)
scaleY = 'auto';              % Automatic selection

% =========================
% Pseudo-Arc-Length Continuation
% =========================
% hmax: maximum allowed change per continuation step
% hmin: minimum step size before stopping
hmax = 0.2;               
hmin = 1e-3;              

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

%% ---- Load state-space model ----
load ssmodel;

% Discrete-time model
Ts=1/fs;

% Number of output channels
ny = size(C,1);

%% ---- Nonlinearities ----

% Nonlinear basis functions
J = size(Be,2)-1;                 % number of nonlinear terms

fNL    = cell(J,1);
fNL{1} = @(x,xd) x(1,:).^3;
fNL{2} = @(x,xd) x(1,:).^2;

% Analytic derivatives 

dFNL = cell(J,1);
dFNL{1,1} = @(y,yd)  3*(y(1,:)).^2;
dFNL{2,1} = @(y,yd)  2*(y(1,:)).^1;

dFNLvel = cell(J,1);  
dFNLvel{1,1} = @(y,yd)  0*(y(1,:));
dFNLvel{2,1} = @(y,yd)  0*(y(1,:));

%% ---- Build HB struct ---

HB = struct();
HB.A  = A; HB.Be = Be;
HB.C  = C; HB.De = De;
HB.Ts = Ts;
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

% Single harmonic contributions
chPlot = 1;     % output channel
hPlot  = 0:3;   % harmonic contribution
hbss_plot_nfrc(out, chPlot, hPlot);

% Maximum response amplitude %
chPlot = 1;       % output channel
hPlot  = "max";   % harmonic contribution
hbss_plot_nfrc(out, chPlot, hPlot);