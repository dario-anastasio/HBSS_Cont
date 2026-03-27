function out = hbss_identify_bifurcations_post(out, opts)
% HBSS_IDENTIFY_BIFURCATIONS_POST Bifurcation analysis 
%
% Steps:
% 1. Sort multipliers to track their evolution
% 2. Detect and classify bifurcations
%
% -------------------------------------------------------------------------
% Author: Dario Anastasio (Politecnico di Torino)
% Part of: HBSS_Cont toolbox
% -------------------------------------------------------------------------

nPoints = length(out.sigma_unsorted);
if nPoints < 2, return; end 

nMult = length(out.sigma_unsorted{1});
rawSigma = cell2mat(cellfun(@(x) x(:), out.sigma_unsorted, 'UniformOutput', false));

% Delete unsorted sigma from out struct
out = rmfield(out,"sigma_unsorted");

% --- STEP 1: MULTIPLIER TRACKING ---
sortedSigma = zeros(nMult, nPoints);
sortedSigma(:, 1) = rawSigma(:, 1);
for i = 2:nPoints
    prevS = sortedSigma(:, i-1);
    currS = rawSigma(:, i);
    distMat = abs(currS - prevS.'); 
    newOrder = zeros(nMult, 1);
    assigned = false(nMult, 1);
    for j = 1:nMult
        dists = distMat(:, j);
        dists(assigned) = Inf;
        [~, bestMatchIdx] = min(dists);
        newOrder(j) = bestMatchIdx;
        assigned(bestMatchIdx) = true;
    end
    sortedSigma(:, i) = currS(newOrder);
end
out.sigma_sorted = sortedSigma;

% --- STEP 2: DETECTION ---
out.bif.Fold  = false(1, nPoints);
out.bif.PD    = false(1, nPoints);
out.bif.NS    = false(1, nPoints);
out.bif.label = categorical(repmat("none", 1, nPoints), ["none","Fold","PD","NS"]);

R = 1.0; 
K = opts.nPersist;
branchType = strings(nMult, 1); % Memory for each branch ("Fold", "PD", "NS")

for i = 2:(nPoints - K)
    for idx = 1:nMult
        val = sortedSigma(idx, i);
        m_prev = abs(sortedSigma(idx, i-1));
        m_curr = abs(val);
        
        % CASE A: STABILITY LOSS (Crossing OUT)
        if (m_prev <= R) && (m_curr > R + opts.tolCross)
            if all(abs(sortedSigma(idx, i:(i+K-1))) > R)
                
                type = classifyBif(val, sortedSigma(:, i), opts);
                branchType(idx) = type; % Save for re-entry
                applyLabel(i, type);
            end
            
        % CASE B: STABILITY RECOVERY (Crossing IN)
        elseif (m_prev > R) && (m_curr <= R - opts.tolCross)
            if all(abs(sortedSigma(idx, i:(i+K-1))) <= R)
                
                % Use memory or classify
                if branchType(idx) ~= ""
                    type = branchType(idx);
                else
                    type = classifyBif(val, sortedSigma(:, i), opts);
                end
                applyLabel(i, type);
                branchType(idx) = ""; % Reset memory
            end
        end
    end
end

% --- Helper function --- %
function type = classifyBif(s, all_s, opts)
    ang = abs(angle(s));
    if ang < opts.tolRad
        % Fold
        type = "Fold";
    elseif ang > (pi - opts.tolRad)
        % Period doubling
        type = "PD";
    else
        % NS: check for a conjugate pair in the current multipliers
        target = conj(s);
        distToConj = abs(all_s - target);
        [minDist, ~] = min(distToConj);
        
        if minDist < opts.tolConj
            type = "NS";
        else
            type = "Unclear"; % do nothing
        end
    end
end

function applyLabel(pos, type)
    if type == "Fold", out.bif.Fold(pos) = true; end
    if type == "PD",   out.bif.PD(pos)   = true; end
    if type == "NS",   out.bif.NS(pos)   = true; end
    out.bif.label(pos) = type;
end
end