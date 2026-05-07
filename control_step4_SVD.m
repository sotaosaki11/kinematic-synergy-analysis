close all; % すべてのグラフウィンドウを閉じる
clear;     % ワークスペース上の変数をすべて消去する

% 機器設定ファイルの読み込み
make_status_bertec; 
make_status_mocap_29; 

%% =========================================================================
% 【プログラム概要: SVDによるキネマティックシナジー解析】
% 歩行の角度軌道データに対して特異値分解（SVD）を行い、多自由度な運動を
% 少数の低次元な協調構造（キネマティックシナジー）に縮約する。
% 1. 全歩行周期（Global）の解析
% 2. 4つの歩行フェーズ（DS1, SS1, DS2, SS2）ごとの局所解析
% を行い、それぞれについて空間基底（各体節の重み付け）、時間基底（シナジーの
% 活性化パターン）、および寄与率を算出・可視化する。
% =========================================================================

%% ===================================================================
% 1. 初期設定（入出力ファイルとパラメータ設定）
% ===================================================================
tTextFilePath   = './'; 
tTextFileName   = 'IBA_12'; 

% 入出力ファイル名の設定
tTextLoadDataName = [tTextFileName, '_step3_planerlaw']; % 前段(Planar Law)の処理データ
tTextSaveMatName  = [tTextFileName, '_step4_SVD'];       % 本プログラムの保存先

% --- モード反転設定 ---
% SVDの特異ベクトルは数学的性質上、符号が反転して出力されることがある。
% 物理的な解釈（屈曲・伸展など）を統一するため、ここで手動で符号を反転させる。
% (1: 反転する, 0: そのまま)
FLIP_MODES_GLOBAL = [1, 0, 1]; 
FLIP_MODES_DS1    = [1, 1]; 
FLIP_MODES_SS1    = [0, 0]; 
FLIP_MODES_DS2    = [1, 1]; 
FLIP_MODES_SS2    = [0, 1]; 

% --- グラフ保存名設定 ---
tTextGraphName1    = [tTextFileName, '_SVD_Fig1_Avg_Contrib'];
tTextGraphName2    = [tTextFileName, '_SVD_Fig2_Avg_Spatial'];
tTextGraphName3    = [tTextFileName, '_SVD_Fig3_Avg_Temporal']; 
tTextGraphName4    = [tTextFileName, '_SVD_Fig4_Avg_3D_Traj']; 
tTextGraphName4_2D = [tTextFileName, '_SVD_Fig4_2D_Proj']; 

tTextGraphName5    = [tTextFileName, '_SVD_Fig5_Phase_DS1'];
tTextGraphName6    = [tTextFileName, '_SVD_Fig6_Phase_SS1'];
tTextGraphName7    = [tTextFileName, '_SVD_Fig7_Phase_DS2'];
tTextGraphName8    = [tTextFileName, '_SVD_Fig8_Phase_SS2'];

flg_graphSave = 1; % グラフ保存フラグ

% --- フェーズごとの描画カラー設定 ---
color_SS1 = [1, 0, 0];                 % 単脚支持1 (赤)
color_SS2 = [0, 90/255, 1];            % 単脚支持2 (青)
color_DS1 = [77/255, 196/255, 1];      % 両脚支持1 (シアン)
color_DS2 = [3/255, 175/255, 122/255]; % 両脚支持2 (緑)

%% ===================================================================
% 2. データの読込と成形
% ===================================================================
fprintf('Loading data from %s.mat...\n', tTextLoadDataName);
load(tTextLoadDataName); 

if ~exist('averaged_stride', 'var') 
    error('averaged_stride が見つかりません。');
end

% データ列の並び替え
% 下半身から体幹へ、左脚から右脚へと論理的な順序でプロットするためのインデックス
angleNames = {'L Foot', 'L Shank', 'L Thigh', 'Trunk', 'R Thigh', 'R Shank', 'R Foot'};
newOrder = [7, 5, 3, 1, 2, 4, 6]; 

%% ===================================================================
% 3. 代表平均データに対する SVD 解析 (Global & Phase)
% ===================================================================
% グラフのY軸範囲を統一するため、描画前にすべての解析を完了させる。

% --- 3a. Global SVD解析 (全歩行周期: 200pts) ---
% SVD: X = U * S * V'
% U: 時間的な活性化パターン (Temporal Modes)
% S: 特異値 (各モードの振幅・重要度)
% V: 各体節の空間的な重み付け (Spatial Modes)
originalData = averaged_stride; 
X_global = originalData(:, newOrder); 
[T_svd, N_svd] = size(X_global);  

fprintf('Performing Global SVD on AVERAGE stride (200pts)...\n');
[U_global, S_global, V_global] = svd(X_global, 'econ');

% プロット用にスケーリングした変数を定義
U_with_lambda = U_global * S_global; % 時間基底 (特異値を掛けて実際の振幅にする)
V_spatial = V_global;                % 空間基底

% Globalモードの符号反転処理
if length(FLIP_MODES_GLOBAL) < N_svd
    FLIP_MODES_GLOBAL(end+1:N_svd) = 0; 
end
for i = 1:N_svd
    if FLIP_MODES_GLOBAL(i) == 1
        V_spatial(:, i) = -V_spatial(:, i);
        U_with_lambda(:, i) = -U_with_lambda(:, i);
    end
end
W_proj = V_spatial(:, 1:3); % 3D軌道投影用の空間基底 (上位3成分)

% --- 3b. Phase SVD解析 (各フェーズごと: 100pts) ---
X_DS1 = averaged_DS1(:, newOrder); X_SS1 = averaged_SS1(:, newOrder);
X_DS2 = averaged_DS2(:, newOrder); X_SS2 = averaged_SS2(:, newOrder);

[V_DS1, S_DS1, U_DS1, normal_DS1, d_DS1, p_DS1_3D, P_DS1_Proj] = analyze_phase_svd(X_DS1, W_proj, FLIP_MODES_DS1);
[V_SS1, S_SS1, U_SS1, normal_SS1, d_SS1, p_SS1_3D, P_SS1_Proj] = analyze_phase_svd(X_SS1, W_proj, FLIP_MODES_SS1);
[V_DS2, S_DS2, U_DS2, normal_DS2, d_DS2, p_DS2_3D, P_DS2_Proj] = analyze_phase_svd(X_DS2, W_proj, FLIP_MODES_DS2);
[V_SS2, S_SS2, U_SS2, normal_SS2, d_SS2, p_SS2_3D, P_SS2_Proj] = analyze_phase_svd(X_SS2, W_proj, FLIP_MODES_SS2);

% --- 3c. グラフの共通軸範囲・目盛りの決定 ---
% 空間基底 (Spatial: V) の共通Y軸範囲
numModesToPlot_main = min(3, N_svd);
max_V_global = max(abs(V_spatial(:, 1:numModesToPlot_main)), [], 'all');
max_V_DS1 = max(abs(V_DS1(:, 1:min(2, size(V_DS1,2)))), [], 'all'); 
max_V_SS1 = max(abs(V_SS1(:, 1:min(2, size(V_SS1,2)))), [], 'all'); 
max_V_DS2 = max(abs(V_DS2(:, 1:min(2, size(V_DS2,2)))), [], 'all'); 
max_V_SS2 = max(abs(V_SS2(:, 1:min(2, size(V_SS2,2)))), [], 'all');

maxAbs_V_common = max([max_V_global, max_V_DS1, max_V_SS1, max_V_DS2, max_V_SS2]); 
if isempty(maxAbs_V_common)||maxAbs_V_common==0, maxAbs_V_common=1; end; 

yLimit_V_common = [-1.1, 1.1]; 
yTicks_V_common = [-1 0 1];

% 時間基底 (Temporal: U*S) の共通Y軸範囲
US_DS1=U_DS1(:,1:min(2,size(S_DS1,1)))*S_DS1(1:min(2,size(S_DS1,1)),1:min(2,size(S_DS1,2))); 
US_SS1=U_SS1(:,1:min(2,size(S_SS1,1)))*S_SS1(1:min(2,size(S_SS1,1)),1:min(2,size(S_SS1,2))); 
US_DS2=U_DS2(:,1:min(2,size(S_DS2,1)))*S_DS2(1:min(2,size(S_DS2,1)),1:min(2,size(S_DS2,2))); 
US_SS2=U_SS2(:,1:min(2,size(S_SS2,1)))*S_SS2(1:min(2,size(S_SS2,1)),1:min(2,size(S_SS2,2)));
max_U_DS1=max(abs(US_DS1),[],'all'); max_U_SS1=max(abs(US_SS1),[],'all'); 
max_U_DS2=max(abs(US_DS2),[],'all'); max_U_SS2=max(abs(US_SS2),[],'all');
maxAbs_U_phase = max([max_U_DS1, max_U_SS1, max_U_DS2, max_U_SS2]); 
if isempty(maxAbs_U_phase)||maxAbs_U_phase==0, maxAbs_U_phase=1; end; 

yLimit_U_phase = [-maxAbs_U_phase*1.1, maxAbs_U_phase*1.1];


%% ===================================================================
% 4. グラフ描画セクション (Global)
% ===================================================================

% グラフ描画用のレイアウトマージン設定
LAYOUT_MARGIN_B = 0.15; % 下余白
LAYOUT_MARGIN_T = 0.05; % 上余白
LAYOUT_MARGIN_L = 0.22; % 左余白
LAYOUT_MARGIN_R = 0.05; % 右余白
LAYOUT_GAP      = 0.05; % パネル間の隙間

% --- [Fig 1] 寄与率 (Contribution Ratio) ---
singularValues = diag(S_global);
varianceExplained = singularValues.^2 / sum(singularValues.^2) * 100;
cumulativeVariance = cumsum(varianceExplained);

fig1 = figure(1); clf(fig1);
axes('Position', [LAYOUT_MARGIN_L, LAYOUT_MARGIN_B, 1-LAYOUT_MARGIN_L-LAYOUT_MARGIN_R, 1-LAYOUT_MARGIN_B-LAYOUT_MARGIN_T]);
hold on;
b = bar(varianceExplained, 'FaceColor', [0.2 0.4 0.6]);
p_line = plot(cumulativeVariance, 'r-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'MarkerSize', 4);

for i = 1:N_svd
    if i > 1 && varianceExplained(i) > 0.1 
        text(i, varianceExplained(i) + 4, sprintf('%.3g%%', varianceExplained(i)), ...
            'Horiz', 'center', 'FontSize', 8, 'FontWeight', 'bold');
    end
    text(i, cumulativeVariance(i) + 6, sprintf('%.3g%%', cumulativeVariance(i)), ...
        'Horiz', 'center', 'Color', 'r', 'FontSize', 8);
end

xlabel('Mode number', 'FontSize', 11); 
ylabel('Ratio (%)', 'FontSize', 11); 
grid on; ax = gca; ax.GridAlpha = 0.3; 
xticks(1:N_svd); xlim([0.5, N_svd+0.5]);
ylim([0, 110]); yticks(0:50:100); 
legend([b, p_line], {'C_k', 'CC_k'}, 'Location', 'southeast', 'FontSize', 10);
hold off;

% --- [Fig 2] 空間基底 (Spatial Modes) ---
fig2 = figure(2); clf(fig2);
h_panel_global = (1 - LAYOUT_MARGIN_B - LAYOUT_MARGIN_T - (numModesToPlot_main-1)*LAYOUT_GAP) / numModesToPlot_main;

for i = 1:numModesToPlot_main
    pos_bottom = LAYOUT_MARGIN_B + (numModesToPlot_main - i) * (h_panel_global + LAYOUT_GAP);
    axes('Position', [LAYOUT_MARGIN_L, pos_bottom, 1-LAYOUT_MARGIN_L-LAYOUT_MARGIN_R, h_panel_global]);
    
    bar(V_spatial(:, i), 'FaceColor', [0.4 0.6 0.8]); 
    ylabel(['{\bf z}_{' num2str(i) '}'], 'Interpreter', 'tex', 'FontSize', 11);
    
    grid on; ylim(yLimit_V_common); yticks(yTicks_V_common);
    ax = gca; ax.XTick = 1:N_svd; ax.GridAlpha = 0.3;
    
    if i < numModesToPlot_main
        set(ax, 'XTickLabel', []);
    else
        if N_svd <= length(angleNames)
            ax.XTickLabel = angleNames; 
            ax.XTickLabelRotation = 45; 
        end
    end
end

% --- [Fig 3] 時間基底 (Temporal Modes) ---
maxAbs_U_global = max(abs(U_with_lambda(:, 1:numModesToPlot_main)), [], 'all');
yLimit_U_global = [-maxAbs_U_global * 1.1, maxAbs_U_global * 1.1];
numPoints = size(X_global, 1);
time_axis = linspace(0, 100, numPoints); % 0から100%の歩行周期時間軸

fig3 = figure(3); clf(fig3);

for i = 1:numModesToPlot_main
    pos_bottom = LAYOUT_MARGIN_B + (numModesToPlot_main - i) * (h_panel_global + LAYOUT_GAP);
    axes('Position', [LAYOUT_MARGIN_L, pos_bottom, 1-LAYOUT_MARGIN_L-LAYOUT_MARGIN_R, h_panel_global]);
    
    plot(time_axis, U_with_lambda(:, i), 'k-', 'LineWidth', 1.8);
    ylabel(['\lambda_{' num2str(i) '}{\bf v}_{' num2str(i) '}'], 'Interpreter', 'tex', 'FontSize', 11);
    
    grid on; xlim([0, 100]); ylim(yLimit_U_global);
    xticks([0, 50, 100]);
    ax = gca; ax.GridAlpha = 0.3;
    
    if i < numModesToPlot_main
        set(ax, 'XTickLabel', []);
        xlabel(''); 
    else
        xlabel('Time [%]', 'FontSize', 10);
    end
end

%% ===================================================================
% 5. 3D軌道の空間投影プロット (Fig 4 & 4_2D)
% ===================================================================
% 各フェーズの軌道を上位3つの空間基底が張る空間（シナジー空間）に投影する
fig4 = figure(4); clf(fig4);
hold on;
hT_DS1 = plot3(P_DS1_Proj(:,1), P_DS1_Proj(:,2), P_DS1_Proj(:,3), '.', 'Color', color_DS1, 'MarkerSize', 8);
hT_SS1 = plot3(P_SS1_Proj(:,1), P_SS1_Proj(:,2), P_SS1_Proj(:,3), '.', 'Color', color_SS1, 'MarkerSize', 8);
hT_DS2 = plot3(P_DS2_Proj(:,1), P_DS2_Proj(:,2), P_DS2_Proj(:,3), '.', 'Color', color_DS2, 'MarkerSize', 8);
hT_SS2 = plot3(P_SS2_Proj(:,1), P_SS2_Proj(:,2), P_SS2_Proj(:,3), '.', 'Color', color_SS2, 'MarkerSize', 8);

all_m1 = [P_DS1_Proj(:,1); P_SS1_Proj(:,1); P_DS2_Proj(:,1); P_SS2_Proj(:,1)];
all_m2 = [P_DS1_Proj(:,2); P_SS1_Proj(:,2); P_DS2_Proj(:,2); P_SS2_Proj(:,2)];
all_m3 = [P_DS1_Proj(:,3); P_SS1_Proj(:,3); P_DS2_Proj(:,3); P_SS2_Proj(:,3)];
mode1_lim = [min(all_m1), max(all_m1)]; mode2_lim = [min(all_m2), max(all_m2)]; mode3_lim = [min(all_m3), max(all_m3)]; 
padding = 0.1;
final_x_lim = [mode1_lim(1)-(mode1_lim(2)-mode1_lim(1))*padding, mode1_lim(2)+(mode1_lim(2)-mode1_lim(1))*padding];
final_y_lim = [mode2_lim(1)-(mode2_lim(2)-mode2_lim(1))*padding, mode2_lim(2)+(mode2_lim(2)-mode2_lim(1))*padding];
final_z_lim = [mode3_lim(1)-(mode3_lim(2)-mode3_lim(1))*padding, mode3_lim(2)+(mode3_lim(2)-mode3_lim(1))*padding];

[grid_x, grid_y] = meshgrid(linspace(final_x_lim(1), final_x_lim(2), 10), linspace(final_y_lim(1), final_y_lim(2), 10));

hP_DS1 = surf(grid_x, grid_y, (d_DS1 - normal_DS1(1)*grid_x - normal_DS1(2)*grid_y)/normal_DS1(3), 'FaceColor', color_DS1, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hP_SS1 = surf(grid_x, grid_y, (d_SS1 - normal_SS1(1)*grid_x - normal_SS1(2)*grid_y)/normal_SS1(3), 'FaceColor', color_SS1, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hP_DS2 = surf(grid_x, grid_y, (d_DS2 - normal_DS2(1)*grid_x - normal_DS2(2)*grid_y)/normal_DS2(3), 'FaceColor', color_DS2, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hP_SS2 = surf(grid_x, grid_y, (d_SS2 - normal_SS2(1)*grid_x - normal_SS2(2)*grid_y)/normal_DS2(3), 'FaceColor', color_SS2, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

hStart = plot3(P_DS1_Proj(1,1), P_DS1_Proj(1,2), P_DS1_Proj(1,3), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
hEnd   = plot3(P_SS2_Proj(end,1), P_SS2_Proj(end,2), P_SS2_Proj(end,3), 'x', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'LineWidth', 2);

hold off; 
xlabel('Mode 1', 'FontSize', 12); ylabel('Mode 2', 'FontSize', 12); zlabel('Mode 3', 'FontSize', 12);
title('3D Trajectory in Synergy Space', 'FontSize', 16); 
legend([hT_DS1, hT_SS1, hT_DS2, hT_SS2, hP_DS1, hP_SS1, hP_DS2, hP_SS2, hStart, hEnd], ...
       'DS1 Traj.', 'SS1 Traj.', 'DS2 Traj.', 'SS2 Traj.', 'DS1 Plane', 'SS1 Plane', 'DS2 Plane', 'SS2 Plane', 'Start', 'End', ...
       'Location', 'bestoutside'); 
grid on; xlim(final_x_lim); ylim(final_y_lim); zlim(final_z_lim); axis vis3d; view(30, 20);

% --- 2D Projections (XY, XZ, YZ平面への投影) ---
fig4_2D = figure(42); clf(fig4_2D); fig4_2D.Position = [100, 100, 600, 1000];
subplot(3, 1, 1); hold on; 
plot(P_DS1_Proj(:,1), P_DS1_Proj(:,2), '.', 'Color', color_DS1, 'MarkerSize', 8); 
plot(P_SS1_Proj(:,1), P_SS1_Proj(:,2), '.', 'Color', color_SS1, 'MarkerSize', 8); 
plot(P_DS2_Proj(:,1), P_DS2_Proj(:,2), '.', 'Color', color_DS2, 'MarkerSize', 8); 
plot(P_SS2_Proj(:,1), P_SS2_Proj(:,2), '.', 'Color', color_SS2, 'MarkerSize', 8); 
hold off; grid on; axis equal; xlabel('Mode 1'); ylabel('Mode 2'); title('Mode 1 vs Mode 2');

subplot(3, 1, 2); hold on; 
plot(P_DS1_Proj(:,1), P_DS1_Proj(:,3), '.', 'Color', color_DS1, 'MarkerSize', 8); 
plot(P_SS1_Proj(:,1), P_SS1_Proj(:,3), '.', 'Color', color_SS1, 'MarkerSize', 8); 
plot(P_DS2_Proj(:,1), P_DS2_Proj(:,3), '.', 'Color', color_DS2, 'MarkerSize', 8); 
plot(P_SS2_Proj(:,1), P_SS2_Proj(:,3), '.', 'Color', color_SS2, 'MarkerSize', 8); 
hold off; grid on; axis equal; xlabel('Mode 1'); ylabel('Mode 3'); title('Mode 1 vs Mode 3');

subplot(3, 1, 3); hold on; 
plot(P_DS1_Proj(:,2), P_DS1_Proj(:,3), '.', 'Color', color_DS1, 'MarkerSize', 8); 
plot(P_SS1_Proj(:,2), P_SS1_Proj(:,3), '.', 'Color', color_SS1, 'MarkerSize', 8); 
plot(P_DS2_Proj(:,2), P_DS2_Proj(:,3), '.', 'Color', color_DS2, 'MarkerSize', 8); 
plot(P_SS2_Proj(:,2), P_SS2_Proj(:,3), '.', 'Color', color_SS2, 'MarkerSize', 8); 
hold off; grid on; axis equal; xlabel('Mode 2'); ylabel('Mode 3'); title('Mode 2 vs Mode 3');

%% ===================================================================
% 6. 各フェーズの SVD 解析結果のグラフ描画 (Fig 5-8)
% ===================================================================
% サブ関数を使用して、各フェーズごとの寄与率・空間基底・時間基底を描画する
[f5c, f5s, f5t] = plot_phase_svd_separated_manual(50, S_DS1, V_DS1, U_DS1, angleNames, yLimit_V_common, yLimit_U_phase, yTicks_V_common);
[f6c, f6s, f6t] = plot_phase_svd_separated_manual(60, S_SS1, V_SS1, U_SS1, angleNames, yLimit_V_common, yLimit_U_phase, yTicks_V_common);
[f7c, f7s, f7t] = plot_phase_svd_separated_manual(70, S_DS2, V_DS2, U_DS2, angleNames, yLimit_V_common, yLimit_U_phase, yTicks_V_common);
[f8c, f8s, f8t] = plot_phase_svd_separated_manual(80, S_SS2, V_SS2, U_SS2, angleNames, yLimit_V_common, yLimit_U_phase, yTicks_V_common);


%% ===================================================================
% 7. 全ストライドに対する個別の SVD 解析 (統計処理用)
% ===================================================================
% 平均軌道だけでなく、すべての試行（ストライド）に対して個別にSVDを適用し、
% 分散や標準偏差などの統計的指標を算出する。
fprintf('Performing SVD on EACH stride (Global & Phase)...\n');
n_strides = size(all_normalized_strides, 3);
stride_svd = struct(); 

% Global用保存配列
rec_lambda = zeros(n_strides, 3); % 特異値 (1st-3rd)
rec_cum    = zeros(n_strides, 3); % 累積寄与率 (1st-3rd)

% Phase用保存配列
stats_DS1_lam = zeros(n_strides, 2); stats_DS1_cum = zeros(n_strides, 2);
stats_SS1_lam = zeros(n_strides, 2); stats_SS1_cum = zeros(n_strides, 2);
stats_DS2_lam = zeros(n_strides, 2); stats_DS2_cum = zeros(n_strides, 2);
stats_SS2_lam = zeros(n_strides, 2); stats_SS2_cum = zeros(n_strides, 2);

for k = 1:n_strides
    % --- (A) Global Analysis (全歩行周期) ---
    X_glob_k = all_normalized_strides(:,:,k);
    X_glob_k = X_glob_k(:, newOrder);
    [U_k, S_k, V_k] = svd(X_glob_k, 'econ');
    
    % 個別ストライドの特異ベクトルが平均ベクトルと同じ向きを向くように補正
    for m = 1:size(V_k, 2)
        if dot(V_k(:,m), V_spatial(:,m)) < 0
            V_k(:,m) = -V_k(:,m);
            U_k(:,m) = -U_k(:,m);
        end
    end
    stride_svd(k).global.S = S_k;
    stride_svd(k).global.explained = diag(S_k).^2 / sum(diag(S_k).^2) * 100;
    
    s_vals = diag(S_k);
    cum_vals = cumsum(s_vals.^2) / sum(s_vals.^2);
    rec_lambda(k, :) = s_vals(1:3)';   
    rec_cum(k, :)    = cum_vals(1:3)'; 

    % --- (B) Phase Analysis (各支持期) ---
    if exist('all_norm_DS1','var')
        X_p = all_norm_DS1(:,:,k); X_p = X_p(:, newOrder);
        [v, s, u, n, d] = analyze_phase_svd_single(X_p, W_proj, V_DS1); 
        s_vec = diag(s); total_var = sum(s_vec.^2); cum_vec = cumsum(s_vec.^2)/total_var;
        stats_DS1_lam(k, :) = s_vec(1:2)'; stats_DS1_cum(k, :) = cum_vec(1:2)';
    end
    
    if exist('all_norm_SS1','var')
        X_p = all_norm_SS1(:,:,k); X_p = X_p(:, newOrder);
        [v, s, u, n, d] = analyze_phase_svd_single(X_p, W_proj, V_SS1); 
        s_vec = diag(s); total_var = sum(s_vec.^2); cum_vec = cumsum(s_vec.^2)/total_var;
        stats_SS1_lam(k, :) = s_vec(1:2)'; stats_SS1_cum(k, :) = cum_vec(1:2)';
    end

    if exist('all_norm_DS2','var')
        X_p = all_norm_DS2(:,:,k); X_p = X_p(:, newOrder);
        [v, s, u, n, d] = analyze_phase_svd_single(X_p, W_proj, V_DS2); 
        s_vec = diag(s); total_var = sum(s_vec.^2); cum_vec = cumsum(s_vec.^2)/total_var;
        stats_DS2_lam(k, :) = s_vec(1:2)'; stats_DS2_cum(k, :) = cum_vec(1:2)';
    end

    if exist('all_norm_SS2','var')
        X_p = all_norm_SS2(:,:,k); X_p = X_p(:, newOrder);
        [v, s, u, n, d] = analyze_phase_svd_single(X_p, W_proj, V_SS2); 
        s_vec = diag(s); total_var = sum(s_vec.^2); cum_vec = cumsum(s_vec.^2)/total_var;
        stats_SS2_lam(k, :) = s_vec(1:2)'; stats_SS2_cum(k, :) = cum_vec(1:2)';
    end
end
fprintf('Stride-by-stride SVD complete.\n');

%% ===================================================================
% 8. 論文表作成用データのコンソール出力
% ===================================================================
fmt_val = @(m, s) sprintf('%.2f ± %.2f', m, s);
fmt_lam = @(m, s) sprintf('%.1f ± %.1f', m, s); 

fprintf('\n==========================================================================\n');
fprintf('   Table Data Output: 1. Full Stride (Global) & 2. Stance Phases   \n');
fprintf('==========================================================================\n');
fprintf('Subject ID: %s\n', tTextFileName);

% --- 1. Global (全歩行周期) の結果出力 ---
mean_lambda = mean(rec_lambda, 1); std_lambda  = std(rec_lambda, 0, 1);
mean_cum    = mean(rec_cum, 1);    std_cum     = std(rec_cum, 0, 1);

fprintf('\n[1] Full Stride (Global SVD)\n');
fprintf('----------------------------------------------------------------\n');
fprintf('Mode\t| Singular Value (Mean ± SD)\t| Cumulative Prop (Mean ± SD)\n');
fprintf('----------------------------------------------------------------\n');
for i = 1:3
    fprintf('%dst/rd\t| %6.1f ± %4.1f\t\t\t| %.2f ± %.2f\n', ...
        i, mean_lambda(i), std_lambda(i), mean_cum(i), std_cum(i));
end
fprintf('----------------------------------------------------------------\n');

% --- 2. Phase (支持期詳細) の結果出力 ---
m_DS1_L = mean(stats_DS1_lam); s_DS1_L = std(stats_DS1_lam); m_DS1_C = mean(stats_DS1_cum); s_DS1_C = std(stats_DS1_cum);
m_SS1_L = mean(stats_SS1_lam); s_SS1_L = std(stats_SS1_lam); m_SS1_C = mean(stats_SS1_cum); s_SS1_C = std(stats_SS1_cum);
m_DS2_L = mean(stats_DS2_lam); s_DS2_L = std(stats_DS2_lam); m_DS2_C = mean(stats_DS2_cum); s_DS2_C = std(stats_DS2_cum);
m_SS2_L = mean(stats_SS2_lam); s_SS2_L = std(stats_SS2_lam); m_SS2_C = mean(stats_SS2_cum); s_SS2_C = std(stats_SS2_cum);

fprintf('\n[2] Stance Phases (DS1, SS1, DS2, SS2)\n');
fprintf('-----------------------------------------------------------------------------------------------------\n');
fprintf(' Phase | Mode | Lambda (Mean±SD) | CumProp (Mean±SD) || Phase | Mode | Lambda (Mean±SD) | CumProp (Mean±SD)\n');
fprintf('-----------------------------------------------------------------------------------------------------\n');
fprintf(' DS1   | 1st  | %s      | %s        || SS1   | 1st  | %s      | %s\n', ...
    fmt_lam(m_DS1_L(1), s_DS1_L(1)), fmt_val(m_DS1_C(1), s_DS1_C(1)), fmt_lam(m_SS1_L(1), s_SS1_L(1)), fmt_val(m_SS1_C(1), s_SS1_C(1)));
fprintf('       | 2nd  | %s      | %s        ||       | 2nd  | %s      | %s\n', ...
    fmt_lam(m_DS1_L(2), s_DS1_L(2)), fmt_val(m_DS1_C(2), s_DS1_C(2)), fmt_lam(m_SS1_L(2), s_SS1_L(2)), fmt_val(m_SS1_C(2), s_SS1_C(2)));
fprintf('-----------------------------------------------------------------------------------------------------\n');
fprintf(' DS2   | 1st  | %s      | %s        || SS2   | 1st  | %s      | %s\n', ...
    fmt_lam(m_DS2_L(1), s_DS2_L(1)), fmt_val(m_DS2_C(1), s_DS2_C(1)), fmt_lam(m_SS2_L(1), s_SS2_L(1)), fmt_val(m_SS2_C(1), s_SS2_C(1)));
fprintf('       | 2nd  | %s      | %s        ||       | 2nd  | %s      | %s\n', ...
    fmt_lam(m_DS2_L(2), s_DS2_L(2)), fmt_val(m_DS2_C(2), s_DS2_C(2)), fmt_lam(m_SS2_L(2), s_SS2_L(2)), fmt_val(m_SS2_C(2), s_SS2_C(2)));
fprintf('-----------------------------------------------------------------------------------------------------\n');
fprintf('※ LaTeX コピー用 (Values for Table):\n');
fprintf('DS1 1st: & %s & %s & %s & %s \n', fmt_lam(m_DS1_L(1), s_DS1_L(1)), fmt_lam(m_DS1_L(2), s_DS1_L(2)), fmt_val(m_DS1_C(1), s_DS1_C(1)), fmt_val(m_DS1_C(2), s_DS1_C(2)));
fprintf('SS1 1st: & %s & %s & %s & %s \n', fmt_lam(m_SS1_L(1), s_SS1_L(1)), fmt_lam(m_SS1_L(2), s_SS1_L(2)), fmt_val(m_SS1_C(1), s_SS1_C(1)), fmt_val(m_SS1_C(2), s_SS1_C(2)));
fprintf('DS2 1st: & %s & %s & %s & %s \n', fmt_lam(m_DS2_L(1), s_DS2_L(1)), fmt_lam(m_DS2_L(2), s_DS2_L(2)), fmt_val(m_DS2_C(1), s_DS2_C(1)), fmt_val(m_DS2_C(2), s_DS2_C(2)));
fprintf('SS2 1st: & %s & %s & %s & %s \n', fmt_lam(m_SS2_L(1), s_SS2_L(1)), fmt_lam(m_SS2_L(2), s_SS2_L(2)), fmt_val(m_SS2_C(1), s_SS2_C(1)), fmt_val(m_SS2_C(2), s_SS2_C(2)));
fprintf('\n');

%% ===================================================================
% 9. 計算結果のデータ保存
% ===================================================================
fprintf('Saving results to %s.mat...\n', tTextSaveMatName);
save(tTextSaveMatName, ...
     'data', 'theta', 'centeredtheta', 'meanposture', ...
     'averaged_stride', 'std_stride', 'all_normalized_strides', ...
     'averaged_DS1', 'averaged_SS1', 'averaged_DS2', 'averaged_SS2', ...
     'all_norm_DS1', 'all_norm_SS1', 'all_norm_DS2', 'all_norm_SS2', ...
     'stride_svd', ... 
     'U_with_lambda', 'V_spatial', 'S_global', 'U_global', 'V_global', ...
     'W_proj', ...
     'V_DS1', 'S_DS1', 'U_DS1', 'normal_DS1', 'd_DS1', 'p_DS1_3D', 'P_DS1_Proj', ...
     'V_SS1', 'S_SS1', 'U_SS1', 'normal_SS1', 'd_SS1', 'p_SS1_3D', 'P_SS1_Proj', ...
     'V_DS2', 'S_DS2', 'U_DS2', 'normal_DS2', 'd_DS2', 'p_DS2_3D', 'P_DS2_Proj', ...
     'V_SS2', 'S_SS2', 'U_SS2', 'normal_SS2', 'd_SS2', 'p_SS2_3D', 'P_SS2_Proj', ...
     'FLIP_MODES_GLOBAL', 'FLIP_MODES_DS1', 'FLIP_MODES_SS1', 'FLIP_MODES_DS2', 'FLIP_MODES_SS2', ...
     'mean_orbit_R', 'std_orbit_R', 'all_normalized_R_orbits', ... 
     'mean_orbit_L', 'std_orbit_L', 'all_normalized_L_orbits', ...
     'mean_orbit_R_centered', 'mean_orbit_L_centered', ...
     'mean_orbit_R_RAW', 'std_orbit_R_RAW', ...
     'mean_orbit_L_RAW', 'std_orbit_L_RAW', ...
     'axis_labels_3D', 'cols_R_Leg', 'cols_L_Leg', ...
     'coeff_R', 'explained_R', 'score_R', ...
     'coeff_L', 'explained_L', 'score_L', ...
     'ratio_R', 'ratio_L', 'area_R', 'area_L', ...
     'mean_normalized_LO_R', 'mean_LO_index_R', ...
     'mean_normalized_LO_L', 'mean_LO_index_L', ...
     'normal_dot_product', ...
     'res_Ratio_R', 'res_Area_R', 'res_Ratio_L', 'res_Area_L', ...
     'res_Planarity_R', 'res_Planarity_L', ...
     'res_DotProd', 'res_Angle');

%% ===================================================================
% 10. グラフ画像の一括保存 (サイズ調整込み)
% ===================================================================
if flg_graphSave
    fprintf('Finalizing figures and saving...\n');
    func_graphSave2(fig1, tTextGraphName1, flg_graphSave, 0, [400, 300]); 
    func_graphSave2(fig2, tTextGraphName2, flg_graphSave, 0, [180, 300]); 
    func_graphSave2(fig3, tTextGraphName3, flg_graphSave, 0, [180, 300]); 
    func_graphSave2(fig4, tTextGraphName4, flg_graphSave, 0, [600, 400]);
    func_graphSave2(fig4_2D, tTextGraphName4_2D, flg_graphSave, 0, [600, 1000]); 
    
    size_common = [270, 250]; 
    
    func_graphSave2(f5c, [tTextGraphName5 '_Contrib'], flg_graphSave, 0, size_common);
    func_graphSave2(f5s, [tTextGraphName5 '_Spatial'], flg_graphSave, 0, size_common);
    func_graphSave2(f5t, [tTextGraphName5 '_Temporal'], flg_graphSave, 0, size_common);
    
    func_graphSave2(f6c, [tTextGraphName6 '_Contrib'], flg_graphSave, 0, size_common);
    func_graphSave2(f6s, [tTextGraphName6 '_Spatial'], flg_graphSave, 0, size_common);
    func_graphSave2(f6t, [tTextGraphName6 '_Temporal'], flg_graphSave, 0, size_common);
    
    func_graphSave2(f7c, [tTextGraphName7 '_Contrib'], flg_graphSave, 0, size_common);
    func_graphSave2(f7s, [tTextGraphName7 '_Spatial'], flg_graphSave, 0, size_common);
    func_graphSave2(f7t, [tTextGraphName7 '_Temporal'], flg_graphSave, 0, size_common);
    
    func_graphSave2(f8c, [tTextGraphName8 '_Contrib'], flg_graphSave, 0, size_common);
    func_graphSave2(f8s, [tTextGraphName8 '_Spatial'], flg_graphSave, 0, size_common);
    func_graphSave2(f8t, [tTextGraphName8 '_Temporal'], flg_graphSave, 0, size_common);
end

fprintf('All analysis and saving processes are complete.\n');

%% ===================================================================
% 11. サブ関数定義
% ===================================================================

% 平均軌道用 Phase SVD解析関数
function [V_sub, S_sub, U_sub, normal, d, p_3D, P_proj] = analyze_phase_svd(X_phase_data, W_global, flip_m)
    X_mean = mean(X_phase_data, 1);
    X_var = X_phase_data - X_mean;
    [U_sub, S_sub, V_sub] = svd(X_var, 'econ');
    for i = 1:min(length(flip_m), size(V_sub, 2))
        if flip_m(i)==1, V_sub(:,i) = -V_sub(:,i); U_sub(:,i) = -U_sub(:,i); end
    end
    v1_3D = (V_sub(:,1)' * W_global);
    v2_3D = (V_sub(:,2)' * W_global);
    normal = cross(v1_3D, v2_3D); normal = normal / norm(normal);
    p_3D = X_mean * W_global;
    d = dot(normal, p_3D);
    P_proj = X_phase_data * W_global;
end

% 個別軌道用 Phase SVD解析関数
function [V_sub, S_sub, U_sub, normal, d] = analyze_phase_svd_single(X_phase_data, W_global, V_ref)
    X_mean = mean(X_phase_data, 1);
    X_var = X_phase_data - X_mean;
    [U_sub, S_sub, V_sub] = svd(X_var, 'econ');
    if nargin >= 3
        for m = 1:size(V_sub, 2)
            if dot(V_sub(:,m), V_ref(:,m)) < 0
                V_sub(:,m) = -V_sub(:,m);
                U_sub(:,m) = -U_sub(:,m);
            end
        end
    end
    v1_3D = (V_sub(:,1)' * W_global);
    v2_3D = (V_sub(:,2)' * W_global);
    normal = cross(v1_3D, v2_3D); normal = normal / norm(normal);
    p_3D = X_mean * W_global;
    d = dot(normal, p_3D);
end

% 分割版グラフ描画関数 (寄与率・空間・時間を個別のFigureに生成)
function [figC, figS, figT] = plot_phase_svd_separated_manual(baseFigIdx, S_matrix, V_matrix, U_matrix, angleNames, yLimit_V, yLimit_U, yTicks_V)
    margin_B = 0.15; 
    margin_T = 0.05; 
    margin_L = 0.22; 
    margin_R = 0.05; 
    gap      = 0.05; 
    
    numPoints = size(U_matrix, 1);
    numModesToPlot_phase = min(2, size(V_matrix, 2));
    N_phase_dim = size(V_matrix, 1);
    s_val = diag(S_matrix);
    var_exp = s_val.^2 / sum(s_val.^2) * 100;
    cum_var = cumsum(var_exp);
    
    % --- 1. 寄与率 (Contribution Ratio) ---
    figC = figure(baseFigIdx + 1); clf(figC);
    axes('Position', [margin_L, margin_B, 1-margin_L-margin_R, 1-margin_B-margin_T]);
    hold on;
    b = bar(var_exp, 'FaceColor', [0.2 0.4 0.6]);
    p_line = plot(cum_var, 'r-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    
    for i = 1:length(s_val)
        if i > 1 && var_exp(i) > 0.1 
            text(i, var_exp(i) + 4, sprintf('%.3g%%', var_exp(i)), ...
                'Horiz', 'center', 'FontSize', 8, 'FontWeight', 'bold');
        end
        text(i, cum_var(i) + 6, sprintf('%.3g%%', cum_var(i)), ...
            'Horiz', 'center', 'Color', 'r', 'FontSize', 8);
    end
    
    xlabel('Mode number', 'FontSize', 11); 
    ylabel('Ratio (%)', 'FontSize', 11);
    
    grid on; ax = gca; ax.GridAlpha = 0.3; xticks(1:length(s_val)); xlim([0.5, length(s_val)+0.5]); 
    ylim([0, 110]); yticks(0:50:100); 
    legend([b, p_line], {'C_k', 'CC_k'}, 'Location', 'southeast', 'FontSize', 8); 
    hold off;

    % --- 2. 空間パターン (Spatial Pattern) ---
    figS = figure(baseFigIdx + 2); clf(figS);
    h_panel = (1 - margin_B - margin_T - (numModesToPlot_phase-1)*gap) / numModesToPlot_phase;
    
    for i = 1:numModesToPlot_phase
        pos_bottom = margin_B + (numModesToPlot_phase - i) * (h_panel + gap);
        axes('Position', [margin_L, pos_bottom, 1-margin_L-margin_R, h_panel]);
        
        bar(V_matrix(:, i), 'FaceColor', [0.4 0.6 0.8]);
        ylabel(['{\bf z}_{' num2str(i) '}'], 'Interpreter', 'tex', 'FontSize', 11);
        
        grid on; ylim(yLimit_V); yticks(yTicks_V);
        ax = gca; ax.XTick = 1:N_phase_dim; ax.GridAlpha = 0.3;
        
        if i < numModesToPlot_phase
            set(ax, 'XTickLabel', []); 
        else
            if N_phase_dim <= length(angleNames)
                ax.XTickLabel = angleNames; 
                ax.XTickLabelRotation = 45; 
            end
        end
    end
    
    % --- 3. 時間パターン (Temporal Pattern) ---
    figT = figure(baseFigIdx + 3); clf(figT);
    
    for i = 1:numModesToPlot_phase
        pos_bottom = margin_B + (numModesToPlot_phase - i) * (h_panel + gap);
        axes('Position', [margin_L, pos_bottom, 1-margin_L-margin_R, h_panel]);
        
        plot(1:numPoints, U_matrix(:,i)*S_matrix(i,i), 'k-', 'LineWidth', 1.8);
        ylabel(['\lambda_{' num2str(i) '}{\bf v}_{' num2str(i) '}'], 'Interpreter', 'tex', 'FontSize', 11);
        
        grid on; xlim([1, numPoints]); ylim(yLimit_U);
        xticks([1, round(numPoints/2), numPoints]);
        
        ax = gca; ax.GridAlpha = 0.3; 
        if i < numModesToPlot_phase
            set(ax, 'XTickLabel', []);
        else
            xlabel('Time [%]', 'FontSize', 10);
        end
    end
end