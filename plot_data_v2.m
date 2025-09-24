%──────────────────────────────────────────────────────────────
%  plot_metrics_style2.m ― 선택 플롯 Style-2 시각화
%──────────────────────────────────────────────────────────────
%  대상 플롯
%   • |Y| vs Time
%   • e_y vs Time
%   • e_\psi vs Time
%   • e_\gamma vs Time          % (Yaw rate error를 γ로 표기)
%   • \psi vs \gamma (Phase Plane)
%   • \beta vs \gamma (Phase Plane)
%   • y vs v_y (Phase Plane)
%   • e_y vs \dot{e}_y (Phase Plane, 보간·스무딩)
%   • Differential Wheel Torque (Front & Rear) vs Time
%──────────────────────────────────────────────────────────────
clear; close all; clc;

%% ═════ USER CONFIG ═════
baseDir        = 'C:\CM_Projects\Integrated_Adaptive_Control_final_v3_SCIE\src_cm4sl\Result';
scenarios      = {'Scenario 1','Scenario 2','Scenario 3','Scenario 4'};
reference_path = 'C:\CM_Projects\Integrated_Adaptive_Control_final_v3_SCIE\src_cm4sl\Reference_Path_data';
outDir  = fullfile(baseDir,'Summary_Figures_Style2');   % 결과 저장 폴더
dpi     = 300;                                          % 해상도
figPos  = [2.54 2.54 41.51 22.72];                      % cm
timeLim = 55;                                           % 데이터 사용 한계 [s]

% (NEW) |Y| 오프셋
Y_OFFSET = 5.25;

% Phase Plane 시간 범위 [s]

phaseTimeRange      = [17.9, 20.8];     % \psi–\gamma, \beta–\gamma 전용

phaseTimeRange      = [4.5, 7.2];     % \psi–\gamma, \beta–\gamma 전용
phaseTimeRangeYVy   = [5, 9];       % y–v_y 전용
phaseTimeRangeEyEyd = [10, 13.3];   % e_y–\dot{e}_y 전용

% 선택적 y축 한계: []면 자동
yLimCfg = struct( ...
  'Y',              [-5 5], ...
  'e_y',            [-2.25 2.25], ...
  'e_psi',          [], ...
  'e_gamma',        [], ...
  'deltaSW',        [], ...
  'phasePsiGamma',  [], ...
  'phaseBetaGamma', [], ...
  'phaseYVy',       [], ...
  'phaseEyEyd',     [], ...
  'torque',         [-400, 800], ...
  'diffTorque',     [-200 200] ...
);

% 컨트롤러별 팔레트/스타일
ctrlCmap  = [1 0 0; 0 0 1; 1 0 1; 0 1 0];   % R, B, M, G
ctrlStyle = {'-','-','-','-'};              % 모두 실선

if ~exist(outDir,'dir'); mkdir(outDir); end
refT = readtable(fullfile(reference_path,'Reference.csv'));
refT = refT(refT.Time<=timeLim,:);

%% ═════ FIGURE TEMPLATE (시나리오 전체 공통) ═════
figY   = initFig(201,'|Y| vs Time');                 % tex
figEy  = initFig(202,'e_y vs Time');                 % tex
figEps = initFig(203,'e_\psi vs Time');              % tex
figEr  = initFig(204,'e_\gamma vs Time');            % tex
figSte = initFig(205,'\delta_{f} vs Time');          % tex

% Phase Plane (Yaw rate=γ 표기)
figPhasePsiGam  = initFig(206,'Phase: \psi vs \gamma');                   % tex
figPhaseBetaGam = initFig(207,'Phase: \beta vs \gamma');                  % tex
figPhaseYVy     = initFig(208,'Phase: y vs v_y');                         % tex
figPhaseEyEyd   = initFig(210,'Phase: {\boldmath$e_y$} vs {\boldmath$\dot{e}_y$}','latex'); % ★ LaTeX

%% ═════ MAIN LOOP (시나리오) ═════
allTags = {};
for s = 1:numel(scenarios)
    rootDir = fullfile(baseDir, scenarios{s});
    cases   = collectCases(rootDir);
    allTags = [allTags, {cases.tag}];

    % 시나리오별 타일 축
    axY   = getTile(figY,  s, scenarios{s});
    axEy  = getTile(figEy, s, scenarios{s});
    axEps = getTile(figEps,s, scenarios{s});
    axEr  = getTile(figEr, s, scenarios{s});
    axSte = getTile(figSte,s, scenarios{s});

    axPhasePsiGam  = getTile(figPhasePsiGam,  s, scenarios{s});
    axPhaseBetaGam = getTile(figPhaseBetaGam, s, scenarios{s});
    axPhaseYVy     = getTile(figPhaseYVy,     s, scenarios{s});
    axPhaseEyEyd   = getTile(figPhaseEyEyd,   s, scenarios{s});   % (LaTeX 라벨 예정)

    % 참조 궤적 (|Y|) ― (NEW) +Y_OFFSET
    plot(axY, refT.Time, refT.Y + Y_OFFSET, '--k', 'LineWidth',2.0);

    %% ── 차동토크 Figure: (1×2) 서브플롯, TV 케이스만 표시 ──
    figDTq = figure(300+s); clf(figDTq);
    set(figDTq,'Units','centimeters','Position',figPos);
    tlDTq  = tiledlayout(figDTq,1,2,'TileSpacing','compact');
    sgtitle(tlDTq, ['Differential Wheel Torque ― ' scenarios{s}], ...
            'FontWeight','bold','Interpreter','tex');

    % 좌/우 축 (왼: 전륜 ΔT_F, 오른: 후륜 ΔT_R)
    axF = nexttile(tlDTq,1); hold(axF,'on'); grid(axF,'on');
    title(axF, '\DeltaT_F = T_{FL} - T_{FR}','Interpreter','tex','FontWeight','bold');
    xlabel(axF,'Time [s]','Interpreter','tex','FontWeight','bold');
    ylabel(axF,'\DeltaT_F [Nm]','Interpreter','tex','FontWeight','bold');
    xlim(axF,[0 timeLim]); applyYLim(axF, yLimCfg.diffTorque); yline(axF,0,'--','LineWidth',1.0,'HandleVisibility','off');

    axR = nexttile(tlDTq,2); hold(axR,'on'); grid(axR,'on');
    title(axR, '\DeltaT_R = T_{RL} - T_{RR}','Interpreter','tex','FontWeight','bold');
    xlabel(axR,'Time [s]','Interpreter','tex','FontWeight','bold');
    ylabel(axR,'\DeltaT_R [Nm]','Interpreter','tex','FontWeight','bold');
    xlim(axR,[0 timeLim]); applyYLim(axR, yLimCfg.diffTorque); yline(axR,0,'--','LineWidth',1.0,'HandleVisibility','off');

    % 범례용 핸들/이름
    hLeg = []; legNames = {};

    %% 케이스 루프 (차동토크: TV 케이스만)
    for k = 1:numel(cases)
        if ~ismember(cases(k).tag, {'LQR_TV_SWA','MRAC_TV_SWA'})
            continue;  % TV 아닌 케이스 스킵
        end

        % 데이터 읽기 & 시간 제한
        absT = readtable(cases(k).abs);
        absT = absT(absT.Time<=timeLim,:);

        % 스타일/색상
        ls = ctrlStyle{mod(k-1,numel(ctrlStyle))+1};
        cc = ctrlCmap (mod(k-1,size(ctrlCmap,1))+1,:);

        % ΔT 계산 + 0~4 s 구간 0 마스킹
        reqCols = {'T_FL','T_FR','T_RL','T_RR'};
        if ~all(ismember(reqCols, absT.Properties.VariableNames))
            warning('Wheel torque columns missing for case %s', cases(k).tag);
            continue;
        end
        dTF = absT.T_FL - absT.T_FR;
        dTR = absT.T_RL - absT.T_RR;
        idxZero = (absT.Time >= 0) & (absT.Time <= 4);
        dTF(idxZero) = 0; dTR(idxZero) = 0;

        % 플로팅
        dispName = getDisplayName(cases(k).tag);
        hF = plot(axF, absT.Time, dTF, 'LineStyle',ls,'Color',cc,'LineWidth',2.0, 'DisplayName', dispName);
        plot(axR, absT.Time, dTR, 'LineStyle',ls,'Color',cc,'LineWidth',2.0, 'DisplayName', dispName);

        hLeg(end+1)     = hF; %#ok<AGROW>
        legNames{end+1} = dispName; %#ok<AGROW>
    end

    % 범례(좌측 축)
    if ~isempty(hLeg)
        lgd = legend(axF, hLeg, legNames, 'Interpreter','none','FontWeight','bold');
        lgd.Orientation = 'horizontal';
        try, lgd.Layout.Tile = 'north'; catch, end
    end

    %% ── 공통 플롯들 ──
    for k = 1:numel(cases)
        absT = readtable(cases(k).abs);
        errT = readtable(cases(k).err);
        absT = absT(absT.Time<=timeLim,:); errT = errT(errT.Time<=timeLim,:);

        ls = ctrlStyle{mod(k-1,numel(ctrlStyle))+1};
        cc = ctrlCmap (mod(k-1,size(ctrlCmap,1))+1,:);

        % |Y| vs Time (+오프셋)
        plot(axY, absT.Time, absT.Y + Y_OFFSET, 'LineStyle',ls,'Color',cc,'LineWidth',2.0);

        % e_y, e_\psi, e_\gamma
        plot(axEy,  errT.Time, errT.e_y,   'LineStyle',ls,'Color',cc,'LineWidth',2.0);
        plot(axEps, errT.Time, errT.e_yaw, 'LineStyle',ls,'Color',cc,'LineWidth',2.0);
        if ismember('e_yawrate', errT.Properties.VariableNames)
            plot(axEr, errT.Time, errT.e_yawrate,'LineStyle',ls,'Color',cc,'LineWidth',2.0);
        end

        % δ_f
        if ismember('Steer', absT.Properties.VariableNames)
            plot(axSte, absT.Time, absT.Steer, 'LineStyle',ls,'Color',cc,'LineWidth',2.0);
        else
            warning('Column ''Steer'' not found for case %s', cases(k).tag);
        end

        % Phase: ψ–γ, β–γ
        idxYaw = absT.Time >= phaseTimeRange(1) & absT.Time <= phaseTimeRange(2);
        if all(ismember({'Yaw','YawRate'}, absT.Properties.VariableNames))
            plot(axPhasePsiGam,  absT.Yaw(idxYaw),           absT.YawRate(idxYaw), 'LineWidth',2.0,'LineStyle',ls,'Color',cc);
        end
        if ismember('SideSlipAngle', absT.Properties.VariableNames)
            plot(axPhaseBetaGam, absT.SideSlipAngle(idxYaw), absT.YawRate(idxYaw), 'LineWidth',2.0,'LineStyle',ls,'Color',cc);
        end

        % Phase: y–v_y
        idxYVy = absT.Time >= phaseTimeRangeYVy(1) & absT.Time <= phaseTimeRangeYVy(2);
        if ismember('Y',absT.Properties.VariableNames) && ismember('Vy',absT.Properties.VariableNames)
            plot(axPhaseYVy, absT.Y(idxYVy), absT.Vy(idxYVy), 'LineWidth',2.0,'LineStyle',ls,'Color',cc);
        end

        % Phase: e_y – \dot{e}_y (보간·스무딩 후 미분)
        if ismember('Time', errT.Properties.VariableNames) && ismember('e_y', errT.Properties.VariableNames)
            t0 = phaseTimeRangeEyEyd(1); t1 = phaseTimeRangeEyEyd(2);
            dtNative = median(diff(errT.Time(~isnan(errT.Time))));
            if isempty(dtNative) || ~isfinite(dtNative) || dtNative<=0, dt = 0.01; else, dt = min(0.01, dtNative/2); end
            tq  = (t0:dt:t1).';
            eyq = interp1(errT.Time, errT.e_y, tq, 'linear', 'extrap');
            win = max(5, 2*floor(0.25/dt)+1);     % 약 0.25 s 창 (홀수)
            eySmooth = smoothdata(eyq, 'movmean', win);
            eydot    = gradient(eySmooth, dt);
            plot(axPhaseEyEyd, eySmooth, eydot, 'LineWidth',2.0,'LineStyle',ls,'Color',cc);
        end
    end

    %% 라벨 & y축 한계 적용
    formatXAxis(axY,  timeLim, 'Y [m]');              applyYLim(axY,  yLimCfg.Y);
    formatXAxis(axEy, timeLim, 'e_y [m]');            applyYLim(axEy, yLimCfg.e_y);
    formatXAxis(axEps,timeLim, 'e_\psi [rad]');       applyYLim(axEps,yLimCfg.e_psi);
    formatXAxis(axEr, timeLim, 'e_{\gamma} [rad/s]'); applyYLim(axEr, yLimCfg.e_gamma);
    formatXAxis(axSte,timeLim, '\delta_{f} [rad]');   applyYLim(axSte,yLimCfg.deltaSW);

    xlabel(axPhasePsiGam,'\psi [rad]','Interpreter','tex','FontWeight','bold');
    ylabel(axPhasePsiGam,'\gamma [rad/s]','Interpreter','tex','FontWeight','bold');  applyYLim(axPhasePsiGam,  yLimCfg.phasePsiGamma);

    xlabel(axPhaseBetaGam,'\beta [rad]','Interpreter','tex','FontWeight','bold');
    ylabel(axPhaseBetaGam,'\gamma [rad/s]','Interpreter','tex','FontWeight','bold'); applyYLim(axPhaseBetaGam, yLimCfg.phaseBetaGamma);

    xlabel(axPhaseYVy,'y [m]','Interpreter','tex','FontWeight','bold');
    ylabel(axPhaseYVy,'v_y [m/s]','Interpreter','tex','FontWeight','bold');         applyYLim(axPhaseYVy,     yLimCfg.phaseYVy);

    % ★ \dot{e}_y 라벨은 LaTeX + boldmath (경고 제거용)
    xlabel(axPhaseEyEyd,'{\boldmath$e_y$} [m]','Interpreter','latex','FontWeight','bold');
    ylabel(axPhaseEyEyd,'{\boldmath$\dot{e}_y$} [m/s]','Interpreter','latex','FontWeight','bold');
    applyYLim(axPhaseEyEyd,   yLimCfg.phaseEyEyd);

    %% 차동토크 Figure 저장
    saveFigures( {figDTq, ['DiffTorque_FR_' scenarios{s}]}, outDir, dpi );
end

%% ═════ LEGENDS (공통 Figure) ═════
uniqueTags = unique(allTags,'stable');
legCtrl    = cellfun(@getDisplayName, uniqueTags,'UniformOutput',false);
legRef     = [{'Reference'}, legCtrl{:}];

addLegend(figY,  legRef);
addLegend(figEy, legCtrl);
addLegend(figEps,legCtrl);
addLegend(figEr, legCtrl);
addLegend(figSte,legCtrl);

addLegend(figPhasePsiGam,  legCtrl);
addLegend(figPhaseBetaGam, legCtrl);
addLegend(figPhaseYVy,     legCtrl);
addLegend(figPhaseEyEyd,   legCtrl);

%% ═════ 공통 Figure 저장 ═════
saveFigures( ...
    {figY,   'Y_vs_Time_Style2'; ...
     figEy,  'e_y_vs_Time_Style2'; ...
     figEps, 'e_psi_vs_Time_Style2'; ...
     figEr,  'e_r_vs_Time_Style_2'; ...
     figSte, 'Steer_vs_Time_Style2'; ...
     figPhasePsiGam,  'Phase_Psi_vs_gamma_Style2'; ...
     figPhaseBetaGam, 'Phase_Beta_vs_gamma_Style2'; ...
     figPhaseYVy,     'Phase_y_vs_vy_Style2'; ...
     figPhaseEyEyd,   'Phase_ey_vs_eyd_Style2' ...
     }, outDir, dpi);

fprintf('\n모든 Figure 저장 완료 → %s\n', outDir);

%% ═════ HEL퍼 FUNCTIONS ═════
function fh = initFig(fNo, sgTitle, intp)
    % intp: 'tex'|'latex' (기본 'tex')
    if nargin < 3 || isempty(intp), intp = 'tex'; end
    figPos = evalin('base','figPos');
    fh = figure(fNo); clf(fh);
    set(fh,'Units','centimeters','Position',figPos);
    tiledlayout(fh,2,2,'TileSpacing','compact');
    sgtitle(sgTitle,'FontWeight','bold','Interpreter',intp);
end

function ax = getTile(figHandle, idx, scenarioName)
    tl = get(figHandle,'Children');
    ax = nexttile(tl, idx); hold(ax,'on'); grid(ax,'on');
    title(ax, scenarioName,'Interpreter','none','FontWeight','bold');
end

function formatXAxis(ax, xMax, yLabelStr)
    xlabel(ax,'Time [s]','Interpreter','tex','FontWeight','bold');
    ylabel(ax,yLabelStr,'Interpreter','tex','FontWeight','bold');
    xlim(ax,[0 xMax]);
end

function applyYLim(ax, lim)
    if ~isempty(lim) && isnumeric(lim) && numel(lim)==2 && all(isfinite(lim))
        ylim(ax, lim);
    end
end

function addLegend(figHandle, entries)
    axesAll = findobj(figHandle,'Type','Axes','-not','Tag','legend');
    if isempty(axesAll), return; end
    lgd = legend(axesAll(end), entries, 'Interpreter','none','FontWeight','bold');
    lgd.Orientation = 'horizontal';
    try, lgd.Layout.Tile = 'north'; catch, end
end

function cases = collectCases(rootDir)
    cases  = struct('tag',{},'abs',{},'err',{},'tire',{});
    eFiles = dir(fullfile(rootDir,'*_Error.csv'));
    for i = 1:numel(eFiles)
        base  = erase(eFiles(i).name,'_Error.csv');
        absP  = fullfile(rootDir,[base '.csv']);
        tireP = fullfile(rootDir,[base '_Tire.csv']);
        if exist(absP,'file') && exist(tireP,'file')
            cases(end+1) = struct('tag',base,'abs',absP,'err',fullfile(rootDir,eFiles(i).name),'tire',tireP);
        end
    end
end

function name = getDisplayName(tag)
    switch tag
        case 'LQR_SWA',     name='LQR (Steering)';
        case 'LQR_TV_SWA',  name='LQR (Steer+TV)';
        case 'MRAC_SWA',    name='MRAC (Steering)';
        case 'MRAC_TV_SWA', name='MRAC (Steer+TV)';
        otherwise,          name=tag;
    end
end

% ===== PNG만 저장 =====
function saveFigures(figCell, outDir, dpi)
    for i = 1:size(figCell,1)
        fh   = figCell{i,1};
        base = figCell{i,2};
        fPng = fullfile(outDir,[base '.png']);
        set(fh,'PaperUnits','centimeters');
        pos = get(fh,'Position');
        set(fh,'PaperPosition',[0 0 pos(3) pos(4)]);
        print(fh, fPng, '-dpng', ['-r' num2str(dpi)]);
        fprintf('✓ %s (PNG)\n', base);
    end
end
