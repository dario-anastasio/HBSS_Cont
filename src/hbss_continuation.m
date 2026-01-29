function out = hbss_continuation(HB)
% HBSS_CONTINUATION  Nonlinear frequency response via Harmonic Balance
%                    with Pseudo-Arc-Length Continuation
%
% This function computes nonlinear frequency response curves (NFRCs)
% using a Harmonic Balance formulation in extended state space, combined
% with a pseudo-arc-length continuation algorithm.
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% The output struct 'out' contains:
%   - out.freqHz, out.Omega
%   - out.Yh      : complex harmonic coefficients (one-sided)
%   - out.Ymax/Ymin : time-domain extrema for all output channels
%   - out.isStable: 1/0 stability information
%   - out.bif : bifurcation classification
%
% -------------------------------------------------------------------------
% Author: Dario Anastasio (Politecnico di Torino)
% Part of: HBSS_Cont toolbox
% -------------------------------------------------------------------------

%% -------------------------------------------------
% Validate HB
% --------------------------------------------------
HB = hbss_validate_model(HB);

% If model is discrete -> convert it to continuous
if HB.ssType=='d'
    sysd = struct('A',HB.A,'Be',HB.Be,'C',HB.C,'De',HB.De);
    sysc = hbss_d2c_zoh(sysd, HB.Ts);

    HB.A  = sysc.A;
    HB.Be = sysc.Be;
    HB.C  = sysc.C;
    HB.De = sysc.De;

    HB.ssType = 'c';
    HB.Ts = 0;
end

Harmonics = HB.Harmonics(:).';
nH = numel(Harmonics)-1;
ny = HB.ny;

% Prepare options for bifurcation detection
bifOpts = struct();
bifOpts.tolCross = 1e-3;   % tolerance for crossing detection (abs(sigma)-1)
% bifOpts.tolIm = 1e-3;          % tolerance to treat multiplier as real
bifOpts.tolRad = 0.1;          % tolerance for phase calculation

%% Frequency limits
wMin = 2*pi*min(HB.fLim);
if wMin==0 
    wMin=1e-3; % It must be >0
end
wMax = 2*pi*max(HB.fLim);

HF = 1; % fundamental for plotting & forcing channel indexing
ihF = find(Harmonics==HF,1);
assert(~isempty(ihF),'Harmonics must include HF=1.');

%%  Auto scaleY

% User can set scaleY as numeric or 'auto'
isAutoScaleY = false;
if isfield(HB,'scaleY')
    if (ischar(HB.scaleY) && strcmpi(HB.scaleY,'auto')) || (isstring(HB.scaleY) && strcmpi(HB.scaleY,"auto"))
        isAutoScaleY = true;
    elseif isnumeric(HB.scaleY) && ~isfinite(HB.scaleY)
        isAutoScaleY = true; % NaN/Inf -> auto
    end
else
    isAutoScaleY = true; % default auto if not provided
end

if isAutoScaleY
    % Frequency grid for linear FRF (rad/s)
    nGrid = 500;
    OmegaGrid = linspace(wMin, wMax, nGrid);

    % Metric over all output channels (RMS across channels at each frequency)
    ampLin = zeros(1,nGrid);

    for jj = 1:nGrid
        Om = OmegaGrid(jj);

        % Evaluate force amplitude
        Fh = hbss_eval_forcing(HB, Om);
        Fext = [Fh(:,ihF); zeros(HB.nNL,1)];

        % FRF matrix and outputs
        Tm = 1i*HF*Om*eye(HB.nx) - HB.A;
        [Lt,Ut] = lu(Tm);
        He = HB.De + ((HB.C/Ut)*(Lt\HB.Be));
        ylin = He*Fext; % ny×1 complex

        % Global metric over all outputs
        ampLin(jj) = sqrt(mean(abs(ylin).^2));   % RMS across channels
    end

    % Robust reference amplitude over frequency
    ampLin = ampLin(isfinite(ampLin) & ampLin>0);
    if isempty(ampLin)
        scaleY_auto = 1e3; % fallback
        % Aref = NaN;
    else
        Aref = prctile(ampLin,95);
        scaleY_auto = 10 / max(Aref, 1e-12);

        % Round to power of 10
        scaleY_auto = 10^round(log10(scaleY_auto));
    end

    HB.scaleY = scaleY_auto;
    out.scaleY = HB.scaleY;
    out.scaleY_isAuto = isAutoScaleY;
    disp(['--- Scaling set to auto (based on linear FRFs). Value = ',num2str(HB.scaleY),'.']);
end

%% Initial guess from linear response at wMin

Omega0 = wMin;
Y0_phys = zeros(ny*(2*nH+1),1);

% Target at HF with zero NL
Fh0 = hbss_eval_forcing(HB, Omega0);
Fext0 = [Fh0(:,ihF); zeros(HB.nNL,1)];
Tm = 1i*HF*Omega0*eye(HB.nx) - HB.A;
[Lt,Ut] = lu(Tm);
He = HB.De + ((HB.C/Ut)*(Lt\HB.Be));
ylin = He*Fext0; 

baseHF = ny + (HF-1)*2*ny;
Y0_phys(baseHF + (1:ny))       = real(ylin);
Y0_phys(baseHF + ny + (1:ny))  = -imag(ylin);

Y0_scaled = HB.scaleY * Y0_phys;

%% Output struct

out = struct();
out.Omega  = [];
out.freqHz = [];
out.Y      = [];
out.Yh     = [];   % complex harmonic coefficients
out.Ymax   = [];
out.Ymin   = [];
out.nIter  = [];
out.sigma  = {};
out.stability = [];
out.isStable  = [];       
out.bif = struct('Fold',[],'PD',[],'NS',[],'label',[]);

%% Plot setup
if HB.plot.enable
    
    figHBSS=figure; clf;
    set(figHBSS,'Color','w','Units','normalized','Position',[0.15 0.15 0.7 0.65]);

    ax1 = subplot(2,1,1); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
    hNFRC = plot(ax1,NaN,NaN,'k-','LineWidth',1.5);
    set(ax1,'xlim',HB.fLim);

    % instability and bifurcations
    hUnst = plot(ax1, NaN, NaN, '.', 'color', 0.6*[1 1 1], 'MarkerSize', 6, 'LineWidth', 1.0); % generic unstable
    hFold = plot(ax1, NaN, NaN, 'ms', 'MarkerSize', 6, 'LineWidth', 1.2); % Fold
    hPD   = plot(ax1, NaN, NaN, 'b^', 'MarkerSize', 6, 'LineWidth', 1.2); % PD
    hNS   = plot(ax1, NaN, NaN, 'gd', 'MarkerSize', 6, 'LineWidth', 1.2); % NS

    xlabel(ax1,'Frequency (Hz)');
    ylabel(ax1,sprintf('|Y^{(1)}_%d|', HB.plot.chOut));
    % legend(ax1, {'NFRC','Unstable','Fold','PD','NS'}, 'Location','northeast');
    
    ax2 = subplot(2,2,3); grid(ax2,'on'); hold(ax2,'on'); box(ax2,'on');
    bar(ax2, Harmonics, 0.*Harmonics);
    xlabel(ax2,'Harmonic');
    ylabel(ax2,'Amplitude');
    hold(ax2,'off');

    ax3 = subplot(2,2,4); hold(ax3,'on'); axis(ax3,'equal'); grid(ax3,'on'); box(ax3,'on');
    th = linspace(0,2*pi,300);
    plot(ax3,cos(th),sin(th),'k');
    xline(0,'k'); yline(0,'k');
    xlabel(ax3,'Re(\sigma)');
    ylabel(ax3,'Im(\sigma)');
    set(ax3,'xtick',[-1 0 1], 'ytick',[-1 0 1]);

    stopBtn = uicontrol('Style','togglebutton','String','Stop',...
        'Units','normalized','Position',[0.01 0.01 0.06 0.05]);
else
    stopBtn = [];
end

%% hmax vector for continuation

nyHB  = ny*(2*nH+1);
nyTot = nyHB + 1;

if isscalar(HB.cont.hmax)
    hmaxVec = HB.cont.hmax * ones(nyTot,1);
else
    hmaxVec = HB.cont.hmax(:);
    assert(numel(hmaxVec)==nyTot, 'HB.cont.hmax must be scalar.');
end
hmin = HB.cont.hmin;

%% First point: solve HB at Omega0
nTry=1;

switch HB.solver.name
    case 'fsolve'
        Ysol = fsolve(@(Ys) hbss_residual(Ys,Omega0,HB), Y0_scaled, HB.solver.opts);
     case 'newton'
        [Ysol, info] = hbss_newton_solver(@(Ys) hbss_residual(Ys,Omega0,HB), ...
                                          Y0_scaled, HB.solver.opts);
        if ~info.converged
            error('First point did not converge with Newton solver at Omega0=%.6g (||R||=%.3e).', ...
                  Omega0, info.resNorm);
        end
    otherwise
        error('Unknown solver "%s".', HB.solver.name);
end

contPoint = [Ysol; Omega0];           % current solution
contTangent = zeros(size(contPoint)); contTangent(end)=1; % initial tangent

% initial step
h = max(hmin, max(hmaxVec));

%% Continuation loop
disp('--- Starting continuation ---');

while contPoint(end) <= wMax

    % Store solution
    out.freqHz(end+1,1)   = contPoint(end)/(2*pi);
    out.Omega(end+1,1)    = contPoint(end);

    k = numel(out.freqHz);
    Y_scaled = contPoint(1:end-1);
    Y_phys   = Y_scaled / HB.scaleY;
    out.Y(:,k)   = Y_phys;
    out.nIter(end+1,1) = nTry;

    % df estimate
    if numel(out.freqHz)>1
        df = out.freqHz(end) - out.freqHz(end-1);
    else
        df = NaN;
    end

    % Build complex coefficients and save
    Yh_cur = hbss_Y_to_Yh(Y_phys, ny, Harmonics);   % (nH+1)×ny
    out.Yh(:,:,k) = Yh_cur;

    % Time-domain max/min for all outputs
    [ymax, ymin] = hbss_maxmin_from_Yh(Yh_cur, Harmonics, contPoint(end), HB.nSamples);
    out.Ymax(:,k) = ymax;
    out.Ymin(:,k) = ymin;

    % Stability and bifurcation detection
    stab = hbss_stability_monodromy(contPoint, HB);
    sigma = stab.sigma;
    out.sigma{end+1} = sigma;
        
    if numel(out.sigma) >= 2
        sigma_prev = out.sigma{end-1};
        sigma_cur  = sigma;
        bif = hbss_classify_bifurcation_crossing(sigma_prev, sigma_cur, bifOpts);
    else
        % no previous point available
        bif = hbss_classify_bifurcation_crossing([], sigma, bifOpts);
    end
    
    % Save bifurcation/stability info in out
    out.isStable(end+1,1)  = bif.stable;
    out.bif.Fold(end+1,1)  = bif.Fold;
    out.bif.PD(end+1,1)    = bif.PD;
    out.bif.NS(end+1,1)    = bif.NS;
    out.bif.label{end+1,1} = char(bif.label);
    
    % Print bifurcation info
    if  ~strcmp(char(bif.label),'none')
        disp([' * ',char(bif.label),' bifurcation detected at ',...
            num2str(out.freqHz(end)),'Hz']);
    end

    % Plot update
    if HB.plot.enable && (mod(k, HB.plot.updateEvery) ==0 || k == 1)
        
        ch = HB.plot.chOut;
        HF = 1;   

        % harmonic coefficients at current point
        Acur   = abs(Yh_cur(:,ch));    % magnitudes
        
        % first plot
        YhFund = squeeze(out.Yh(HF+1, ch, :));   % complex vector
        ampFund = abs(YhFund);
        set(hNFRC,'XData', out.freqHz, 'YData', ampFund);
        title(ax1, sprintf( ...
            'f = %.4f Hz   |   df = %.3e Hz   |   iterations = %d', ...
            out.freqHz(end), df, out.nIter(end)));

        % second plot
        cla(ax2);
        bar(ax2, Harmonics, Acur);
        xlabel(ax2,'Harmonic');
        ylabel(ax2,'Amplitude');
        ylim(ax2,[0 max(out.Ymax(ch,end))]); grid(ax2,'on');

        % third plot
        cla(ax3);
        plot(ax3,cos(th),sin(th),'k'); hold(ax3,'on');
        xline(0,'k'); yline(0,'k');
        if bif.stable
            plot(ax3,real(sigma),imag(sigma),'k.','MarkerSize',12);
        else
            plot(ax3,real(sigma),imag(sigma),'.','color',0.6*[1 1 1],'MarkerSize',12);
        end

        % update NFRC markers using bifLabel
        f = out.freqHz(:);
        a = abs(squeeze(out.Yh(HF+1, ch, :)));
        lab = string(out.bif.label(:));
        set(hFold, 'XData', f(lab=="Fold"),     'YData', a(lab=="Fold"));
        set(hPD,   'XData', f(lab=="PD"),       'YData', a(lab=="PD"));
        set(hNS,   'XData', f(lab=="NS"),       'YData', a(lab=="NS"));
        set(hUnst, 'XData', f(out.isStable==0), 'YData', a(out.isStable==0));

        % drawnow;
        drawnow limitrate;
        
        if ~isempty(stopBtn) && get(stopBtn,'Value')==1
            break;
        end
    end

    % Continuation step limiting
    Dy = h * contTangent;
    coeff = 1;
    for i = 1:numel(Dy)
        if Dy(i) > hmaxVec(i)
            coeff = min(coeff, hmaxVec(i) / Dy(i));
        end
    end
    h = max(hmin, coeff * h);

    % Corrector retries
    flagCont = 0; 
    nTry = 0;
    while flagCont ~= 1
        nTry = nTry + 1;

        % Predictor
        yP = contPoint + h*contTangent;

        % Corrector
        switch HB.solver.name
            case 'fsolve'
                [contPointNew,~,flagSolver] = fsolve(@(yy) hbss_continuation_residual(yy, yP, contTangent, HB), ...
                                                     yP, HB.solver.opts);
        
            case 'newton'
                [contPointNew, info] = hbss_newton_solver(@(yy) hbss_continuation_residual(yy, yP, contTangent, HB), ...
                                                          yP, HB.solver.opts);
                flagSolver = double(info.converged);
                if flagSolver == 0
                    flagSolver = -1; % mimic fsolve failure
                end
        
            otherwise
                error('Unknown solver "%s".', HB.solver.name);
        end

        if flagSolver >= 1
            % Compute new secant direction (approximate tangent)
            contTangentNew = (contPointNew - contPoint);
            contTangentNew = contTangentNew / norm(contTangentNew);

            % Tangent consistency check
            if (contTangent.'*contTangentNew <= 0)
                flagCont = 0;            % reject continuation point
                h = max(h/2, hmin);      % reduce step size
            else
                % Step length control
                deltaY = (contPointNew - contPoint);
                if all(deltaY < (hmaxVec + 1e-14))
                    flagCont = 1;        % accept continuation point
                else
                    flagCont = 0;        % reject and retry
                    h = max(h/2, hmin);  % reduce step size
                end
            end
        else
            % Newton did not converge
            h = max(h/2, hmin);          % reduce and retry
        end

        if nTry >= HB.cont.maxIterCont
            % Too many iterations
            YmaxAll = out.Ymax(:);
            YmaxAll = YmaxAll(isfinite(YmaxAll) & YmaxAll > 0);
        
            if isempty(YmaxAll)
                suggestScale = HB.scaleY;
            else
                ratio = mean(HB.fLim) / mean(YmaxAll);
                suggestScale = 10^floor(log10(abs(ratio)));
            end
            
            disp([
                '--- Continuation stopped at f = ', num2str(out.freqHz(end,1),'%.4f'), ' Hz.']);
            disp(['Suggestions:', ...
                ' (1) Try a different scaleY factor (suggested order: ', num2str(suggestScale), ');', ...
                ' (2) Reduce the continuation step cont.hmax;', ...
                ' (3) Increase cont.maxIterCont;', ...
                ' (4) Reduce the number of harmonics.']);
        
            if HB.plot.enable
                title(ax1, sprintf( ...
                    'f = %.4f Hz   |   df = %.3e Hz   |   iterations = %d', ...
                    out.freqHz(end), df, nTry));                
            end              
            return;      
        end

        % Accept
        contPoint = contPointNew;
        contTangent = contTangentNew;
    
        % Increase step after success
        h = min(2*h, max(hmaxVec));
    end
end

disp('--- End ---');

if HB.plot.enable
    close(figHBSS);
end