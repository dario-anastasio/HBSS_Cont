function hbss_plot_nfrc(out, chOut, hPlot, varargin)
% HBSS_PLOT_NFRC Plot NFRC with stability/bifurcation markers.
%
%   hbss_plot_nfrc(out, chOut, hPlot)
%   hbss_plot_nfrc(out, chOut, hPlot, lineSpec, color, lineWidth)
%
% Inputs
%   out     : output struct from hbss_continuation
%   chOut   : output channel index/indices to plot (e.g. 1 or [1 2 3])
%   hPlot   : harmonic index/indices to plot (e.g. 1 or [1 3 5]), OR:
%             - "max" : plot out.Ymax(chOut,:) (time-domain maximum)
%             - "min" : plot out.Ymin(chOut,:) (time-domain minimum)
%
% Behavior
%   - Multiple channels: one figure per channel.
%   - Multiple harmonics: amplitude only (no phase), overlaid with different colors.
%     Unstable points use marker '.' with a lighter shade of the corresponding curve.
%     Fold/PD/NS markers are drawn for each curve, but appear only once in legend.
%   - If hPlot is "max" or "min": harmonic selection is disabled, but multiple
%     channels are allowed.
%
% Optional inputs (single-harmonic or max/min)
%   lineSpec  : line style (default '-')
%   color     : line color (default 'k') [ignored for multi-harmonic]
%   lineWidth : line width (default 1.5)
%
% Markers:
%   - Unstable: '.' with lighter shade of the NFRC curve color
%   - Fold: magenta square (ms)
%   - PD: blue triangle (b^)
%   - NS: green diamond (gd)

%% ---------------- Default parameters ----------------
lineSpec  = '-';
lineColor = [0 0 0];
lineWidth = 1.5;

if nargin >= 4 && ~isempty(varargin{1}), lineSpec  = varargin{1}; end
if nargin >= 5 && ~isempty(varargin{2}), lineColor = varargin{2}; end
if nargin >= 6 && ~isempty(varargin{3}), lineWidth = varargin{3}; end

% Utility: unstable branches use lighter shade of the corresponding curve
lightenColor   = @(c,alpha) c + alpha*(1-c);
unstableAlpha  = 0.6;     % 0=no change, 1=white
unstableMSize  = 8;

%% ---------------- Checks ----------------
assert(isfield(out,'freqHz') && ~isempty(out.freqHz), 'out must contain freqHz.');
assert(isfield(out,'bif') && isfield(out.bif,'label'), ...
    'out must contain out.bif.label (strings: Fold/PD/NS/none).');
assert(isfield(out,'isStable'), 'out must contain out.isStable (0/1).');

f = out.freqHz(:);
lab = string(out.bif.label(:));

isFold = (lab == "Fold");
isPD   = (lab == "PD");
isNS   = (lab == "NS");
isUn   = (out.isStable(:) == 0);

% allow multiple channels
chVec = chOut(:).';
assert(all(chVec >= 1) && all(mod(chVec,1)==0), 'chOut must contain positive integers.');

%% ---------------- Interpret hPlot ----------------
isMaxMin = false;
modeStr  = "";

isMultiHarm = false;

if ischar(hPlot) || isstring(hPlot)
    hp = lower(string(hPlot));
    if hp == "max"
        isMaxMin = true;
        modeStr = "max";
        assert(isfield(out,'Ymax') && ~isempty(out.Ymax), 'Requested hPlot="max" but out.Ymax is missing/empty.');
        titStr = 'HBSS\_Cont – NFRC (time-domain maximum)';
    elseif hp == "min"
        isMaxMin = true;
        modeStr = "min";
        assert(isfield(out,'Ymin') && ~isempty(out.Ymin), 'Requested hPlot="min" but out.Ymin is missing/empty.');
        titStr = 'HBSS\_Cont – NFRC (time-domain minimum)';
    else
        error('hPlot must be a harmonic index (numeric), a vector of harmonics, or "max"/"min".');
    end
else
    % numeric harmonic(s)
    assert(isfield(out,'Yh') && ~isempty(out.Yh), 'Harmonic plot requested but out.Yh is missing/empty.');
    hVec = hPlot(:).';
    assert(all(mod(hVec,1)==0) && all(hVec >= 0), 'hPlot must be nonnegative integer harmonic index/indices.');
    assert(max(hVec)+1 <= size(out.Yh,1), 'Requested harmonic index exceeds available harmonics.');
    modeStr = "harmonic";
    isMultiHarm = (numel(hVec) > 1);
end

%% ---------------- Plot loop over channels ----------------
for ch = chVec

    if isMaxMin
        % ---------- Single plot (max/min, no phase) ----------
        if modeStr == "max"
            assert(ch <= size(out.Ymax,1), 'Requested channel exceeds number of outputs in out.Ymax.');
            yPlot = out.Ymax(ch,:).';
            yLab  = sprintf('Y_{max,%d}', ch);
        else
            assert(ch <= size(out.Ymin,1), 'Requested channel exceeds number of outputs in out.Ymin.');
            yPlot = out.Ymin(ch,:).';
            yLab  = sprintf('Y_{min,%d}', ch);
        end

        figure('Color','w');
        ax1 = axes(); hold(ax1,'on'); grid(ax1,'on');

        hNFRC = plot(ax1, f, yPlot, ...
            'LineStyle', lineSpec, ...
            'Color', lineColor, ...
            'LineWidth', lineWidth);

        % unstable points: lighter shade of lineColor
        cUn = lightenColor(lineColor, unstableAlpha);
        hUn = plot(ax1, f(isUn), yPlot(isUn), ...
            '.', 'Color', cUn, 'MarkerSize', unstableMSize, 'LineWidth', 1.0);

        hFold = plot(ax1, f(isFold), yPlot(isFold), 'ms', 'MarkerSize', 6, 'LineWidth', 1.2);
        hPD   = plot(ax1, f(isPD),   yPlot(isPD),   'b^', 'MarkerSize', 6, 'LineWidth', 1.2);
        hNS   = plot(ax1, f(isNS),   yPlot(isNS),   'gd', 'MarkerSize', 6, 'LineWidth', 1.2);

        xlabel(ax1,'Frequency (Hz)');
        ylabel(ax1, yLab);
        title(ax1, titStr);

        % Dynamic legend
        hLeg = hNFRC;
        legTxt = {'NFRC'};

        if any(isFold), hLeg(end+1) = hFold; legTxt{end+1} = 'Fold'; end
        if any(isPD),   hLeg(end+1) = hPD;   legTxt{end+1} = 'PD';   end
        if any(isNS),   hLeg(end+1) = hNS;   legTxt{end+1} = 'NS';   end
        if any(isUn),   hLeg(end+1) = hUn;   legTxt{end+1} = 'Unstable'; end

        legend(ax1, hLeg, legTxt, 'Location','best');

    else
        % ---------- Harmonic plotting ----------
        if ~isMultiHarm
            % ----- Single harmonic: amplitude + phase -----
            hh = hVec(1);
            Yh = squeeze(out.Yh(hh+1, ch, :));  % complex coeff
            amp = abs(Yh);
            ph  = angle(Yh);

            figure('Color','w');

            % Amplitude
            ax1 = subplot(2,1,1); hold(ax1,'on'); grid(ax1,'on');

            hNFRC = plot(ax1, f, amp, ...
                'LineStyle', lineSpec, ...
                'Color', lineColor, ...
                'LineWidth', lineWidth);

            % unstable points: lighter shade of lineColor
            cUn = lightenColor(lineColor, unstableAlpha);
            hUn = plot(ax1, f(isUn), amp(isUn), ...
                '.', 'Color', cUn, 'MarkerSize', unstableMSize, 'LineWidth', 1.0);

            hFold = plot(ax1, f(isFold), amp(isFold), 'ms', 'MarkerSize', 6, 'LineWidth', 1.2);
            hPD   = plot(ax1, f(isPD),   amp(isPD),   'b^', 'MarkerSize', 6, 'LineWidth', 1.2);
            hNS   = plot(ax1, f(isNS),   amp(isNS),   'gd', 'MarkerSize', 6, 'LineWidth', 1.2);

            xlabel(ax1,'Frequency (Hz)');
            ylabel(ax1, ['|Y^{(',num2str(hh),')}_{',num2str(ch),'}|']);
            title(ax1,'HBSS\_Cont – NFRC');

            % Dynamic legend
            hLeg = hNFRC;
            legTxt = {'NFRC'};

            if any(isFold), hLeg(end+1) = hFold; legTxt{end+1} = 'Fold'; end
            if any(isPD),   hLeg(end+1) = hPD;   legTxt{end+1} = 'PD';   end
            if any(isNS),   hLeg(end+1) = hNS;   legTxt{end+1} = 'NS';   end
            if any(isUn),   hLeg(end+1) = hUn;   legTxt{end+1} = 'Unstable'; end

            legend(ax1, hLeg, legTxt, 'Location','best');

            % Phase
            ax2 = subplot(2,1,2); hold(ax2,'on'); grid(ax2,'on');

            plot(ax2, f, ph, ...
                'LineStyle', lineSpec, ...
                'Color', lineColor, ...
                'LineWidth', lineWidth);

            % unstable points on phase: lighter shade of lineColor
            plot(ax2, f(isUn), ph(isUn), ...
                '.', 'Color', cUn, 'MarkerSize', unstableMSize, 'LineWidth', 1.0);

            % Same bifurcation markers on phase
            plot(ax2, f(isFold), ph(isFold), 'ms', 'MarkerSize', 6, 'LineWidth', 1.2);
            plot(ax2, f(isPD),   ph(isPD),   'b^', 'MarkerSize', 6, 'LineWidth', 1.2);
            plot(ax2, f(isNS),   ph(isNS),   'gd', 'MarkerSize', 6, 'LineWidth', 1.2);

            xlabel(ax2,'Frequency (Hz)');
            ylabel(ax2,'Phase (rad)');

        else
            % ----- Multiple harmonics: amplitude only, overlaid -----
            nCurves = numel(hVec);
            cols = turbo(nCurves);

            figure('Color','w');
            ax1 = axes(); hold(ax1,'on'); grid(ax1,'on');

            % legend handles/text
            hLeg = gobjects(0);
            legTxt = {};

            % We'll add "Fold/PD/NS/Unstable" only once in legend, but plot for each curve
            hFoldLeg = gobjects(1); hPDLeg = gobjects(1); hNSLeg = gobjects(1); hUnLeg = gobjects(1);
            hasFold = false; hasPD = false; hasNS = false; hasUn = false;

            for i = 1:nCurves
                hh = hVec(i);
                Yh  = squeeze(out.Yh(hh+1, ch, :));
                amp = abs(Yh);

                cMain = cols(i,:);
                cUn   = lightenColor(cMain, unstableAlpha);

                % main curve
                hMain = plot(ax1, f, amp, ...
                    'LineStyle', lineSpec, ...
                    'Color', cMain, ...
                    'LineWidth', lineWidth);

                hLeg(end+1) = hMain; %#ok<AGROW>
                legTxt{end+1} = sprintf('h=%d', hh); %#ok<AGROW>

                % unstable points for this curve (lighter shade, marker '.')
                if any(isUn)
                    hUnTmp = plot(ax1, f(isUn), amp(isUn), ...
                        '.', 'Color', cUn, 'MarkerSize', unstableMSize, 'LineWidth', 1.0);
                    if ~hasUn
                        hUnLeg = hUnTmp; hasUn = true;
                    end
                end

                % bifurcation markers for this curve
                if any(isFold)
                    hTmp = plot(ax1, f(isFold), amp(isFold), ...
                        'ms', 'MarkerSize', 6, 'LineWidth', 1.2);
                    if ~hasFold
                        hFoldLeg = hTmp; hasFold = true;
                    end
                end
                if any(isPD)
                    hTmp = plot(ax1, f(isPD), amp(isPD), ...
                        'b^', 'MarkerSize', 6, 'LineWidth', 1.2);
                    if ~hasPD
                        hPDLeg = hTmp; hasPD = true;
                    end
                end
                if any(isNS)
                    hTmp = plot(ax1, f(isNS), amp(isNS), ...
                        'gd', 'MarkerSize', 6, 'LineWidth', 1.2);
                    if ~hasNS
                        hNSLeg = hTmp; hasNS = true;
                    end
                end
            end

            xlabel(ax1,'Frequency (Hz)');
            ylabel(ax1, sprintf('|Y^{(h)}_{%d}|', ch));
            title(ax1,'HBSS\_Cont – NFRC');

            % Add bif/unstable legend entries only once (if present)
            if hasFold, hLeg(end+1) = hFoldLeg; legTxt{end+1} = 'Fold'; end
            if hasPD,   hLeg(end+1) = hPDLeg;   legTxt{end+1} = 'PD';   end
            if hasNS,   hLeg(end+1) = hNSLeg;   legTxt{end+1} = 'NS';   end
            if hasUn,   hLeg(end+1) = hUnLeg;   legTxt{end+1} = 'Unstable'; end

            legend(ax1, hLeg, legTxt, 'Location','best');
        end
    end
end
end
