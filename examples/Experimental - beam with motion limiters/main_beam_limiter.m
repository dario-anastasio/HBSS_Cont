%% examples/main_beam_limiter.m
% HBSS_Cont - Example: Cantilever beam with motion limiters
% State-space model obtained experimetnally with NFR-ID method

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
Harmonics = 0:5;           % Harmonic indices [0 1 ... H]
nSamples  = 2048;          % Time samples for AFT

% =========================
% Frequency range
% =========================
fLim = [8 16];             % Frequency range [Hz]

% =========================
% Forcing
% =========================
F0 = 0.015;                 % Amplitude of harmonic forcing

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
load ssmodel_NFRID;

% Discrete-time model
Ts=1/fs;

% Number of output channels
ny = size(C,1);

%% ---- Nonlinearities ----

% Nonlinear basis functions: vector of clearances
gap_vector = (0.5:0.2:1.2)*1e-3;        % (m)
J = length(gap_vector);                 % number of nonlinear terms
fNL = cell(J,1);
for ij = 1:J
    fNL{ij} = @(x,xd) (x(1,:)-gap_vector(ij)*sign(x(1,:))).*(abs(x(1,:))>gap_vector(ij));
end

% Analytic derivatives

dFNL    = cell(J,1);   
dFNLvel = cell(J,1);   

for ij = 1:J
    g = gap_vector(ij);
    dFNL{ij,1}    = @(x,xd) (abs(x(1,:)) > g);
    dFNLvel{ij,1} = @(x,xd) 0.*x(1,:);
end

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

% Fundamental harmonic response
chPlot = 1;     % output channel
hPlot  = 1;     % harmonic contribution
hbss_plot_nfrc(out, chPlot, hPlot);

% Maximum
chPlot = 1;     % output channel
hPlot  = "max"; % maximum
hbss_plot_nfrc(out, chPlot, hPlot);