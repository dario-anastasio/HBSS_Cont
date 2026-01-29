function bif = hbss_classify_bifurcation_crossing(sigma_prev, sigma_cur, opts)
% HBSS_CLASSIFY_BIFURCATION_CROSSING  Bifurcation detection using Floquet multipliers
%
% -------------------------------------------------------------------------
% Author: Dario Anastasio (Politecnico di Torino)
% Part of: HBSS_Cont toolbox
% -------------------------------------------------------------------------

if nargin < 3, opts = struct(); end
if ~isfield(opts,'tolCross'), opts.tolCross = 1e-3; end % Minimum jump to confirm crossing
% if ~isfield(opts,'tolIm'),        opts.tolIm = 1e-3;        end % Threshold for real part
if ~isfield(opts,'tolRad'),       opts.tolRad = 0.1;        end % Threshold for phase angle detection

tolCross = opts.tolCross;
% tolIm    = opts.tolIm;
tolRad   = opts.tolRad;

sigma_prev = sigma_prev(:);
sigma_cur  = sigma_cur(:);

% Overall stability test (current point)
stable = all(abs(sigma_cur) <= 1 + 1e-5); 

% Initialize bifurcation flags
Fold = false; PD = false; NS = false;
label = "none";

% --- Case 1: No previous step -> no bifurcation tracking ---
if isempty(sigma_prev)
    bif = struct('stable',1,'Fold',false,'PD',false,'NS',false,'label',categorical("none",["none","Fold","PD","NS"]));
    return;
end

% --- Case 2: Crossing-based Classification ---
nCur = numel(sigma_cur);
matchPrevIdx = nan(nCur,1);
idxPrevAvailable = 1:numel(sigma_prev);

% Nearest-neighbor (track multipliers trajectory)
for i = 1:nCur
    [~, m_idx] = min(abs(sigma_cur(i) - sigma_prev(idxPrevAvailable)));
    matchPrevIdx(i) = idxPrevAvailable(m_idx);
    idxPrevAvailable(m_idx) = [];
    if isempty(idxPrevAvailable), break; end
end

% Detect crossings
crossIdx = false(nCur,1);
for i = 1:nCur
    j = matchPrevIdx(i);
    if isnan(j), continue; end
    
    g_prev = abs(sigma_prev(j)) - 1;
    g_cur  = abs(sigma_cur(i))  - 1;
    
    % If they have different signs (up to tolerance) -> crossing
    if (g_prev * g_cur < 0) && (abs(g_cur - g_prev) > tolCross)
        crossIdx(i) = true;
    end
end

% Return if no crossing
if ~any(crossIdx)
    bif = struct('stable',stable,'Fold',false,'PD',false,'NS',false,'label',categorical("none",["none","Fold","PD","NS"]));
    return;
end

% --- Bifurcation Classification (check phase-angle) ---

crossed_list = find(crossIdx);

for i = crossed_list.'
    s_cur = sigma_cur(i);
    
    % Calculate the absolute phase angle (from 0 to pi)
    % A multiplier near +1 has an angle near 0.
    % A multiplier near -1 has an angle near pi.
    phase_angle = abs(angle(s_cur)); 
    
    % 1. FOLD: Crossing occurs near the positive real axis (+1)
    % tolRad rad tolerance
    if phase_angle < tolRad
        Fold = true; 
        label = "Fold";
        break; % Priority found
        
    % 2. PERIOD DOUBLING (PD): Crossing occurs near the negative real axis (-1)
    elseif phase_angle > (pi - tolRad)
        PD = true; 
        label = "PD";
        break;
        
    % 3. NEIMARK-SACKER (NS): Crossing occurs elsewhere (complex pair)
    else
        NS = true; 
        label = "NS";
        break;
    end
end

% Create output structure
bif = struct();
bif.stable = stable;
bif.Fold   = Fold;
bif.PD     = PD;
bif.NS     = NS;
bif.label  = categorical(label, ["none","Fold","PD","NS"]);

end