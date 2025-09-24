%% Tyre Friction Circles — R2023b (shared labels centred, PNG only)
clear; clc;close all;

% ── 기본 설정 ─────────────────────────────────────────────
rootDir    = "C:\CM_Projects\Integrated_Adaptive_Control_final_v3_SCIE\src_cm4sl\Result";
timeWin    = [10 20];   sampleStep = 1;
saveDir    = fullfile(rootDir,"Figures"); if ~exist(saveDir,'dir'); mkdir(saveDir); end

scNames    = ["Scenario 1","Scenario 2","Scenario 3","Scenario 4"];
ctrlFiles  = ["LQR_SWA","LQR_TV_SWA","MRAC_SWA","MRAC_TV_SWA"];
ctrlNames  = ["LQR (Steering)","LQR (Steering + Torque Vectoring)", ...
              "MRAC (Steering)","MRAC (Steering + Torque Vectoring)"];
wheelNames = ["FR","FL","RR","RL"];
unitTheta  = linspace(0,2*pi,361);
col        = lines(numel(wheelNames));

% ── Figure & Layout ──────────────────────────────────────
fig = figure('Units','centimeters','Position',[2 2 50 35]);
tlo = tiledlayout(fig,4,4,'TileSpacing','compact','Padding','loose');
title(tlo,'Tire Friction Circles', ...
      'FontWeight','bold','FontSize',14);

hSc = gobjects(1,numel(wheelNames));   % 범례용 핸들

for si = 1:4
    for ci = 1:4
        ax = nexttile(tlo,(si-1)*4+ci); hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');
        plot(ax,cos(unitTheta),sin(unitTheta),'k:','LineWidth',1);     % μFz = 1

        fp = fullfile(rootDir,"Scenario "+si,ctrlFiles(ci)+"_Tire.csv");
        if ~isfile(fp), continue, end
        T = readtable(fp);
        T = T(T.Time>=timeWin(1) & T.Time<=timeWin(2),:);

        for wi = 1:numel(wheelNames)
            w   = wheelNames(wi); idx = 1:sampleStep:height(T);
            Fx  = T{idx,w+"_x"};  Fy  = T{idx,w+"_y"};
            Fz  = abs(T{idx,w+"_z"});  mu = T{idx,"Road_"+w+"_mu"};
            h   = scatter(ax, Fx./(mu.*Fz), Fy./(mu.*Fz), ...
                          14, col(wi,:), 'filled','MarkerFaceAlpha',0.65);
            if si==1 && ci==1, hSc(wi)=h; end
        end
        xlim(ax,[-1.2 1.2]); ylim(ax,[-1.2 1.2]);
        set(ax,'TickDir','out','FontSize',9);
        if si==1, title(ax,ctrlNames(ci),'FontWeight','bold'); end
        if ci==1, ylabel(ax,scNames(si),'FontWeight','bold','FontSize',11);  end
    end
end

% ── 범례 (타이틀 바로 아래) ──────────────────────────────
lgd = legend(hSc(isgraphics(hSc)), wheelNames, 'Orientation','horizontal');
lgd.Layout.Tile = 'north';
if isprop(lgd,'NumColumns'), lgd.NumColumns = numel(wheelNames); end

% ── 공통 레이블을 annotation 으로 배치 (0-1 범위) ───────
annotation(fig,'textbox',[0.1 0.425 0.14 0.14], ...  % x=0.02, y=0.5
           'String','Lateral Friction [-]','Rotation',90, ...
           'LineStyle','none','FontWeight','bold','FontSize',15, ...
           'HorizontalAlignment','center','VerticalAlignment','middle');

annotation(fig,'textbox',[0.43 0.03 0.14 0.05], ...   % x=0.5, y=0.03
           'String','Longitudinal Friction [-]', ...
           'LineStyle','none','FontWeight','bold','FontSize',15, ...
           'HorizontalAlignment','center','VerticalAlignment','middle');

drawnow limitrate;  % flush UI events

% ── PNG 저장 ─────────────────────────────────────────────
exportgraphics(fig,fullfile(saveDir,'FrictionCircle_Tiled.png'),'Resolution',300);
disp("✓ PNG saved → "+saveDir);
