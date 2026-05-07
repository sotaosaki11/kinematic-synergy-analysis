% =========================================================================
% 【プログラム概要: 片麻痺歩行用 平面則 (Planar Law) 解析 (Step 3)】
% 前段（Step 2）でセグメンテーションしたストライドデータに対し、
% 主成分分析 (PCA) および 特異値分解 (SVD) を用いてキネマティックシナジーを解析する。
% 
% 1. 大腿・下腿・足部の3体節が形成する「共変動面 (Covariation Plane)」の算出。
% 2. 健側（Right）と患側（Left）の軌道・平面の3D可視化および比較。
% 3. 平面性の度合い (Planarity)、軌道の形状 (Ratio, Area)、および
%    左右の対称性 (法線ベクトルの内積) を個別のストライドごとに統計解析する。
% =========================================================================
close all; 
clear;
clc;

%% ===================================================================
% 1. 初期設定 (入出力ファイルと解析パラメータ)
% ===================================================================
% --- 解析対象の被験者ID ---
SubjectName = 'KM_ID3';

% --- 入出力ファイル名 (バケツリレー方式) ---
tTextLoadDataName = [SubjectName, '_step2_segmented'];  % 読み込みファイル
outputMatName     = [SubjectName, '_step3_planerlaw'];  % 保存ファイル名

% --- グラフ保存名 ---
tTextGraphName_Fig1_L        = [SubjectName, '_Step3_Fig1_2D_Angles_Left']; 
tTextGraphName_Fig2_R        = [SubjectName, '_Step3_Fig2_2D_Angles_Right']; 
tTextGraphName_Fig3_3D_Split = [SubjectName, '_Step3_Fig3_3D_Orbits_Split']; 
tTextGraphName_Fig4_PCA_Split= [SubjectName, '_Step3_Fig4_PCA_Split'];
tTextGraphName_Fig5_3D_Comb  = [SubjectName, '_Step3_Fig5_3D_Orbits_Combined']; 
tTextGraphName_Fig6_PCA_Comb = [SubjectName, '_Step3_Fig6_PCA_Combined'];       
tTextGraphName_Fig7_Stats    = [SubjectName, '_Step3_Fig7_Planar_Stats']; 

flg_graphSave = 1; % グラフ保存フラグ

%% ===================================================================
% 2. 描画・解析パラメータの設定
% ===================================================================
N_RESAMPLE_POINTS = 101; % 1歩行周期を101点 (0~100%) にリサンプリング

param_ViewAngle = [160, 10]; % 3Dグラフの初期カメラ視点 (Azimuth, Elevation)

% 3D軌道描画に使用する関節インデックス
cols_R_Leg = [2, 4, 6]; % 右脚 (R-Thigh, R-Shank, R-Foot)
cols_L_Leg = [3, 5, 7]; % 左脚 (L-Thigh, L-Shank, L-Foot)
axis_labels_3D = {'Thigh Angle (deg)', 'Shank Angle (deg)', 'Foot Angle (deg)'};

% --- 基本カラー設定 ---
color_Plane_R = [1, 0.7, 0.7]; 
color_Plane_L = [0.7, 0.7, 1]; 
color_Shank = [0, 0.6, 0];   color_Thigh = [0, 0, 0.8];   color_Foot  = [0.8, 0, 0];
color_Shank_Shade = [0.7, 1, 0.7]; color_Thigh_Shade = [0.7, 0.7, 1]; color_Foot_Shade  = [1, 0.7, 0.7]; 
color_Stance = [0, 0, 0.8]; color_Swing  = [0, 0.7, 0]; 

% --- 左右比較用カラー設定 (健側: 赤, 患側: 青) ---
col_Traj_R = [0.8, 0, 0]; 
col_Traj_L = [0, 0, 0.8]; 
col_Plane_R_Face = col_Traj_R; 
col_Plane_L_Face = col_Traj_L; 

param_Fig1_AspectRatio = [1.6, 1, 1]; % 2D時系列グラフの縦横比

%% ===================================================================
% 3. セグメンテーション済みデータの読み込み
% ===================================================================
fprintf('Loading data from "%s.mat" ...\n', tTextLoadDataName);
if exist([tTextLoadDataName, '.mat'], 'file')
    load([tTextLoadDataName, '.mat']); 
else
    error('ファイルが見つからない: %s.mat\nstep2を実行すること。', tTextLoadDataName);
end

if ~exist('valid_strides_R', 'var') || ~exist('valid_strides_L', 'var')
    error('個別ストライドデータ (valid_strides_R/L) が見つからない。');
end

num_strides_R = length(valid_strides_R);
num_strides_L = length(valid_strides_L);
fprintf('Using Pre-Segmented Strides:\n  Right: %d strides\n  Left : %d strides\n', num_strides_R, num_strides_L);

%% ===================================================================
% 4. 3D軌道の正規化と主成分分析 (PCA)
% ===================================================================
% 平面則 (Planar Law) の抽出には、平均姿勢を引いた「角度変動成分」を用いる。
fprintf('Calculating CENTERED orbits for 3D/PCA from Valid Strides...\n');
resampled_time_vector = linspace(0, 1, N_RESAMPLE_POINTS)';

% --- 健側 (Right Leg) の処理 ---
all_normalized_R_orbits = zeros(N_RESAMPLE_POINTS, 3, num_strides_R);
LO_pct_list_R = zeros(num_strides_R, 1);

% 統計用構造体 (ストライドごとのPCA結果を格納)
indiv_R = struct('coeff', zeros(3,3,num_strides_R), 'score', zeros(N_RESAMPLE_POINTS,3,num_strides_R), ...
                 'explained', zeros(3,num_strides_R), 'ratio', zeros(num_strides_R,1), 'area', zeros(num_strides_R,1), ...
                 'planarity', zeros(num_strides_R,1));

for i = 1:num_strides_R
    st = valid_strides_R(i);
    n_len = size(st.raw_stride, 1);
    
    all_normalized_R_orbits(:,:,i) = interp1(linspace(0,1,n_len), st.raw_stride(:, cols_R_Leg), resampled_time_vector, 'pchip');
    LO_pct_list_R(i) = st.pct_LO;

    % ストライド単位のPCA実行
    traj_cent = all_normalized_R_orbits(:,:,i);
    [coeff, score, ~, ~, explained] = pca(traj_cent - mean(traj_cent,1));
    indiv_R.coeff(:,:,i)   = coeff;
    indiv_R.score(:,:,i)   = score;
    indiv_R.explained(:,i) = explained(1:3);
    indiv_R.ratio(i)       = explained(1) / explained(2);
    indiv_R.area(i)        = polyarea(score(:,1), score(:,2));
    indiv_R.planarity(i)   = explained(1) + explained(2);
end

% 平均軌道に対するPCA
mean_orbit_R = mean(all_normalized_R_orbits, 3); 
std_orbit_R  = std(all_normalized_R_orbits, 0, 3); 
mean_orbit_R_centered = mean_orbit_R - mean(mean_orbit_R, 1); 
[coeff_R, score_R, ~, ~, explained_R] = pca(mean_orbit_R_centered);

mean_normalized_LO_R = mean(LO_pct_list_R); 
mean_LO_index_R = round(mean_normalized_LO_R * (N_RESAMPLE_POINTS - 1)) + 1; 

% --- 患側 (Left Leg) の処理 ---
all_normalized_L_orbits = zeros(N_RESAMPLE_POINTS, 3, num_strides_L);
LO_pct_list_L = zeros(num_strides_L, 1);

indiv_L = struct('coeff', zeros(3,3,num_strides_L), 'score', zeros(N_RESAMPLE_POINTS,3,num_strides_L), ...
                 'explained', zeros(3,num_strides_L), 'ratio', zeros(num_strides_L,1), 'area', zeros(num_strides_L,1), ...
                 'planarity', zeros(num_strides_L,1));

for i = 1:num_strides_L
    st = valid_strides_L(i);
    n_len = size(st.raw_stride, 1);
    
    all_normalized_L_orbits(:,:,i) = interp1(linspace(0,1,n_len), st.raw_stride(:, cols_L_Leg), resampled_time_vector, 'pchip');
    LO_pct_list_L(i) = st.pct_LO;

    traj_cent = all_normalized_L_orbits(:,:,i);
    [coeff, score, ~, ~, explained] = pca(traj_cent - mean(traj_cent,1));
    indiv_L.coeff(:,:,i)   = coeff;
    indiv_L.score(:,:,i)   = score;
    indiv_L.explained(:,i) = explained(1:3);
    indiv_L.ratio(i)       = explained(1) / explained(2);
    indiv_L.area(i)        = polyarea(score(:,1), score(:,2));
    indiv_L.planarity(i)   = explained(1) + explained(2);
end

mean_orbit_L = mean(all_normalized_L_orbits, 3); 
std_orbit_L  = std(all_normalized_L_orbits, 0, 3); 
mean_orbit_L_centered = mean_orbit_L - mean(mean_orbit_L, 1); 
[coeff_L, score_L, ~, ~, explained_L] = pca(mean_orbit_L_centered);

mean_normalized_LO_L = mean(LO_pct_list_L); 
mean_LO_index_L = round(mean_normalized_LO_L * (N_RESAMPLE_POINTS - 1)) + 1; 

%% ===================================================================
% 5. 2Dグラフ用 生角度(Absolute Theta)軌道の計算
% ===================================================================
% PCAには用いないが、時系列の変動を直感的に確認するために実際の角度を計算する。
fprintf('Calculating RAW (Absolute) orbits for Fig 1 & 2...\n');

% 健側 (Right)
all_normalized_R_orbits_RAW = zeros(N_RESAMPLE_POINTS, 3, num_strides_R);
for i = 1:num_strides_R
    st = valid_strides_R(i);
    n_len = size(st.raw_stride_theta, 1);
    all_normalized_R_orbits_RAW(:,:,i) = interp1(linspace(0,1,n_len), st.raw_stride_theta(:, cols_R_Leg), resampled_time_vector, 'pchip');
end
mean_orbit_R_RAW = mean(all_normalized_R_orbits_RAW, 3);
std_orbit_R_RAW  = std(all_normalized_R_orbits_RAW, 0, 3);

% 患側 (Left)
all_normalized_L_orbits_RAW = zeros(N_RESAMPLE_POINTS, 3, num_strides_L);
for i = 1:num_strides_L
    st = valid_strides_L(i);
    n_len = size(st.raw_stride_theta, 1);
    all_normalized_L_orbits_RAW(:,:,i) = interp1(linspace(0,1,n_len), st.raw_stride_theta(:, cols_L_Leg), resampled_time_vector, 'pchip');
end
mean_orbit_L_RAW = mean(all_normalized_L_orbits_RAW, 3);
std_orbit_L_RAW  = std(all_normalized_L_orbits_RAW, 0, 3);

% 後方互換性(エイリアス)の定義
mean_orbit_R_cent = mean_orbit_R_centered;
mean_orbit_L_cent = mean_orbit_L_centered;
mean_orbit_R_raw  = mean_orbit_R_RAW;
mean_orbit_L_raw  = mean_orbit_L_RAW;
std_orbit_R_raw   = std_orbit_R_RAW;
std_orbit_L_raw   = std_orbit_L_RAW;

% 左右のグラフでスケールを統一するための共通Y軸範囲の算出
all_vals_L = [mean_orbit_L_RAW + std_orbit_L_RAW; mean_orbit_L_RAW - std_orbit_L_RAW];
all_vals_R = [mean_orbit_R_RAW + std_orbit_R_RAW; mean_orbit_R_RAW - std_orbit_R_RAW];
global_max_Y = max([all_vals_L(:); all_vals_R(:)]);
global_min_Y = min([all_vals_L(:); all_vals_R(:)]);
margin_Y = (global_max_Y - global_min_Y) * 0.1;
common_ylim = [global_min_Y - margin_Y, global_max_Y + margin_Y];

%% ===================================================================
% 6. グラフ生成: [Fig 1 & Fig 2] 2D 時系列仰角
% ===================================================================
time_axis = linspace(0, 100, N_RESAMPLE_POINTS)'; 
plot_x_axis_fill = [time_axis; flipud(time_axis)]; 

% --- [Fig 1] 患側 (Left Leg) ---
fprintf('Plotting Fig 1: 2D Time Series (Left)...\n');
fig_1 = figure(1); clf(fig_1);
set(fig_1, 'Name', '2D Angles Left', 'NumberTitle', 'off');
movegui(fig_1, 'west'); 

hold on; grid on;
fill(plot_x_axis_fill, [mean_orbit_L_RAW(:,1)+std_orbit_L_RAW(:,1); flipud(mean_orbit_L_RAW(:,1)-std_orbit_L_RAW(:,1))], color_Thigh_Shade, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
fill(plot_x_axis_fill, [mean_orbit_L_RAW(:,2)+std_orbit_L_RAW(:,2); flipud(mean_orbit_L_RAW(:,2)-std_orbit_L_RAW(:,2))], color_Shank_Shade, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
fill(plot_x_axis_fill, [mean_orbit_L_RAW(:,3)+std_orbit_L_RAW(:,3); flipud(mean_orbit_L_RAW(:,3)-std_orbit_L_RAW(:,3))], color_Foot_Shade, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hT_L = plot(time_axis, mean_orbit_L_RAW(:,1), 'Color', color_Thigh, 'LineWidth', 2);
hS_L = plot(time_axis, mean_orbit_L_RAW(:,2), 'Color', color_Shank, 'LineWidth', 2);
hF_L = plot(time_axis, mean_orbit_L_RAW(:,3), 'Color', color_Foot, 'LineWidth', 2);
ylabel('Angle (deg)'); xlabel('Gait Cycle (%)'); xlim([0 100]); 
ylim(common_ylim); 
pbaspect(param_Fig1_AspectRatio); 
legend([hT_L, hS_L, hF_L], axis_labels_3D, 'Location', 'best'); 

% --- [Fig 2] 健側 (Right Leg) ---
fprintf('Plotting Fig 2: 2D Time Series (Right)...\n');
fig_2 = figure(2); clf(fig_2);
set(fig_2, 'Name', '2D Angles Right', 'NumberTitle', 'off');
movegui(fig_2, 'east'); 

hold on; grid on;
fill(plot_x_axis_fill, [mean_orbit_R_RAW(:,1)+std_orbit_R_RAW(:,1); flipud(mean_orbit_R_RAW(:,1)-std_orbit_R_RAW(:,1))], color_Thigh_Shade, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
fill(plot_x_axis_fill, [mean_orbit_R_RAW(:,2)+std_orbit_R_RAW(:,2); flipud(mean_orbit_R_RAW(:,2)-std_orbit_R_RAW(:,2))], color_Shank_Shade, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
fill(plot_x_axis_fill, [mean_orbit_R_RAW(:,3)+std_orbit_R_RAW(:,3); flipud(mean_orbit_R_RAW(:,3)-std_orbit_R_RAW(:,3))], color_Foot_Shade, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hT_R = plot(time_axis, mean_orbit_R_RAW(:,1), 'Color', color_Thigh, 'LineWidth', 2);
hS_R = plot(time_axis, mean_orbit_R_RAW(:,2), 'Color', color_Shank, 'LineWidth', 2);
hF_R = plot(time_axis, mean_orbit_R_RAW(:,3), 'Color', color_Foot, 'LineWidth', 2);
ylabel('Angle (deg)'); xlabel('Gait Cycle (%)'); xlim([0 100]); 
ylim(common_ylim); 
pbaspect(param_Fig1_AspectRatio); 
legend([hT_R, hS_R, hF_R], axis_labels_3D, 'Location', 'best'); 

%% ===================================================================
% 7. グラフ生成: [Fig 3 & Fig 4] 左右個別の3D空間とPCA平面への射影
% ===================================================================
% --- [Fig 3] 3D Orbit Split ---
fprintf('Plotting Fig 3: 3D Orbit Split...\n');
fig_3 = figure(3); clf(fig_3);
set(fig_3, 'Name', '3D Split', 'NumberTitle', 'off');
movegui(fig_3, 'center'); 

pad_ratio = 1.1; 
pos_L = [0.07, 0.15, 0.38, 0.75]; pos_R = [0.50, 0.15, 0.38, 0.75]; 

% 患側 (Left Leg)
axL = subplot('Position', pos_L); hold on; grid on; axis equal; view(param_ViewAngle);
orbit_L = mean_orbit_L_centered; nL = coeff_L(:, 3); % 第3主成分(法線)
min_L = min(orbit_L); max_L = max(orbit_L); center_L = (min_L + max_L) / 2;
max_span = max(max_L - min_L) * pad_ratio; half_span = max_span / 2;
lims_x = [center_L(1)-half_span, center_L(1)+half_span];
lims_y = [center_L(2)-half_span, center_L(2)+half_span];
lims_z = [center_L(3)-half_span, center_L(3)+half_span];
[gx, gy] = meshgrid(lims_x, lims_y); gz = (-nL(1)*gx - nL(2)*gy)/nL(3);
surf(gx, gy, gz, 'FaceColor', color_Plane_L, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
plot3(orbit_L(1:mean_LO_index_L,1), orbit_L(1:mean_LO_index_L,2), orbit_L(1:mean_LO_index_L,3), '-', 'Color', color_Stance, 'LineWidth', 2);
plot3(orbit_L(mean_LO_index_L:end,1), orbit_L(mean_LO_index_L:end,2), orbit_L(mean_LO_index_L:end,3), '-', 'Color', color_Swing, 'LineWidth', 2);
plot3(orbit_L(1,1), orbit_L(1,2), orbit_L(1,3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
plot3(orbit_L(mean_LO_index_L,1), orbit_L(mean_LO_index_L,2), orbit_L(mean_LO_index_L,3), 'ks', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
xlim(lims_x); ylim(lims_y); zlim(lims_z);
xlabel(axis_labels_3D{1}); ylabel(axis_labels_3D{2}); zlabel(axis_labels_3D{3});
set(gca, 'XDir', 'reverse', 'YDir', 'reverse'); 

% 健側 (Right Leg)
axR = subplot('Position', pos_R); hold on; grid on; axis equal; view(param_ViewAngle);
orbit_R = mean_orbit_R_centered; nR = coeff_R(:, 3);
min_R = min(orbit_R); max_R = max(orbit_R); center_R = (min_R + max_R) / 2;
max_span = max(max_R - min_R) * pad_ratio; half_span = max_span / 2;
lims_x = [center_R(1)-half_span, center_R(1)+half_span];
lims_y = [center_R(2)-half_span, center_R(2)+half_span];
lims_z = [center_R(3)-half_span, center_R(3)+half_span];
[gx, gy] = meshgrid(lims_x, lims_y); gz = (-nR(1)*gx - nR(2)*gy)/nR(3);
surf(gx, gy, gz, 'FaceColor', color_Plane_R, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
hSt = plot3(orbit_R(1:mean_LO_index_R,1), orbit_R(1:mean_LO_index_R,2), orbit_R(1:mean_LO_index_R,3), '-', 'Color', color_Stance, 'LineWidth', 2);
hSw = plot3(orbit_R(mean_LO_index_R:end,1), orbit_R(mean_LO_index_R:end,2), orbit_R(mean_LO_index_R:end,3), '-', 'Color', color_Swing, 'LineWidth', 2);
hTD = plot3(orbit_R(1,1), orbit_R(1,2), orbit_R(1,3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
hLO = plot3(orbit_R(mean_LO_index_R,1), orbit_R(mean_LO_index_R,2), orbit_R(mean_LO_index_R,3), 'ks', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
lgd = legend([hSt, hSw, hTD, hLO], {'Stance Phase', 'Swing Phase', 'TD', 'LO'}, 'Location', 'best');
set(axR, 'Position', pos_R); xlim(lims_x); ylim(lims_y); zlim(lims_z);
xlabel(axis_labels_3D{1}); ylabel(axis_labels_3D{2}); zlabel(axis_labels_3D{3});
set(gca, 'XDir', 'reverse', 'YDir', 'reverse'); 

% --- [Fig 4] PCA Projection Split ---
% 抽出された平面（第1主成分と第2主成分）に投影されたループ形状を描画する。
fprintf('Plotting Fig 4: PCA Projection Split...\n');
fig_4 = figure(4); clf(fig_4);
set(fig_4, 'Name', 'PCA Split', 'NumberTitle', 'off');
movegui(fig_4, 'center'); 

pos_L = [0.15, 0.55, 0.70, 0.40]; pos_R = [0.15, 0.08, 0.70, 0.40]; 

% 患側
axL = subplot('Position', pos_L); hold on; grid on; axis equal;
plot(score_L(:,1), score_L(:,2), 'k-');
plot(score_L(1:mean_LO_index_L,1), score_L(1:mean_LO_index_L,2), '-', 'Color', color_Stance, 'LineWidth', 2);
plot(score_L(mean_LO_index_L:end,1), score_L(mean_LO_index_L:end,2), '-', 'Color', color_Swing, 'LineWidth', 2);
plot(score_L(1,1), score_L(1,2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
plot(score_L(mean_LO_index_L,1), score_L(mean_LO_index_L,2), 'ks', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
xlabel(sprintf('PC1 (%.1f %%)', explained_L(1))); ylabel(sprintf('PC2 (%.1f %%)', explained_L(2))); 

% 健側
axR = subplot('Position', pos_R); hold on; grid on; axis equal;
plot(score_R(:,1), score_R(:,2), 'k-');
hSt = plot(score_R(1:mean_LO_index_R,1), score_R(1:mean_LO_index_R,2), '-', 'Color', color_Stance, 'LineWidth', 2);
hSw = plot(score_R(mean_LO_index_R:end,1), score_R(mean_LO_index_R:end,2), '-', 'Color', color_Swing, 'LineWidth', 2);
hTD = plot(score_R(1,1), score_R(1,2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
hLO = plot(score_R(mean_LO_index_R,1), score_R(mean_LO_index_R,2), 'ks', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
legend([hSt, hSw, hTD, hLO], {'Stance', 'Swing', 'TD', 'LO'}, 'Location', 'best');
set(axR, 'Position', pos_R); 
xlabel(sprintf('PC1 (%.1f %%)', explained_R(1))); ylabel(sprintf('PC2 (%.1f %%)', explained_R(2))); 

%% ===================================================================
% 8. グラフ生成: [Fig 5 & Fig 6] 左右同軸での比較 (Combined)
% ===================================================================
% 健側と患側の平面の「ずれ」を視覚的に評価するため、同一の3D空間にプロットする。
% 平面の境界が自然な四角形として描画されるよう、ボックスの交差計算（最適化）を行う。
fprintf('Plotting Fig 5: Combined 3D Orbits...\n');
fig_5 = figure(5); clf(fig_5);
set(fig_5, 'Name', '3D Combined', 'NumberTitle', 'off');
movegui(fig_5, 'center'); 
hold on; grid on; axis equal; view(param_ViewAngle);

orbit_R = mean_orbit_R_centered; orbit_L = mean_orbit_L_centered; 

% データの包含範囲から基準スパンを決定
all_pts = [orbit_R; orbit_L];
min_pt = min(all_pts); max_pt = max(all_pts);
center_pt = (min_pt + max_pt) / 2;
base_span = (max_pt - min_pt) * 1.1; 
max_axis_span = max(base_span);
current_span_x = max_axis_span; current_span_y = max_axis_span; current_span_z = max_axis_span;

% [幾何学的最適化] 平面が途切れずに四角形として表示されるまで、仮想的な描画ボックスを拡張する
fprintf('  Optimizing box size to ensure quadrangular plane frames...\n');
max_iter = 50; expand_step = 1.05; 
legs_data = struct('coeff', {coeff_R, coeff_L}, 'col', {col_Traj_R, col_Traj_L});
final_intersect_R = []; final_intersect_L = [];

for iter = 1:max_iter
    half_x = current_span_x / 2; half_y = current_span_y / 2; half_z = current_span_z / 2;
    box_min = center_pt - [half_x, half_y, half_z]; box_max = center_pt + [half_x, half_y, half_z];
    
    corners = [
        box_min(1), box_min(2), box_min(3); box_max(1), box_min(2), box_min(3);
        box_min(1), box_max(2), box_min(3); box_max(1), box_max(2), box_min(3);
        box_min(1), box_min(2), box_max(3); box_max(1), box_min(2), box_max(3);
        box_min(1), box_max(2), box_max(3); box_max(1), box_max(2), box_max(3);
    ];
    edges = [1,2; 3,4; 5,6; 7,8; 1,3; 2,4; 5,7; 6,8; 1,5; 2,6; 3,7; 4,8];
    
    is_hex_R = false; is_hex_L = false;
    current_intersects = cell(1, 2);
    for lg = 1:2
        n = legs_data(lg).coeff(:, 3); 
        pts = [];
        for i = 1:size(edges, 1)
            p1 = corners(edges(i,1), :); p2 = corners(edges(i,2), :);
            vec = p2 - p1; denom = dot(n, vec);
            if abs(denom) > 1e-6
                t = -dot(n, p1) / denom;
                if t >= 0 && t <= 1, pts = [pts; p1 + t * vec]; end
            end
        end
        if ~isempty(pts)
            pts = unique(round(pts, 5), 'rows');
            if size(pts, 1) > 2
                pc1 = legs_data(lg).coeff(:, 1); pc2 = legs_data(lg).coeff(:, 2);
                projs_x = pts * pc1; projs_y = pts * pc2;
                angles = atan2(projs_y, projs_x);
                [~, sort_idx] = sort(angles);
                pts = pts(sort_idx, :); pts(end+1, :) = pts(1, :);
            end
        end
        current_intersects{lg} = pts;
        if size(pts, 1) - 1 > 4
            if lg == 1, is_hex_R = true; else, is_hex_L = true; end
        end
    end
    if ~is_hex_R && ~is_hex_L
        final_intersect_R = current_intersects{1}; final_intersect_L = current_intersects{2};
        fprintf('    Converged at iteration %d. Frame is quadrangular.\n', iter);
        break;
    else
        current_span_y = current_span_y * expand_step;
    end
end

% 描画スケールの設定
lims_x = [center_pt(1)-current_span_x/2, center_pt(1)+current_span_x/2];
lims_y = [center_pt(2)-current_span_y/2, center_pt(2)+current_span_y/2];
lims_z = [center_pt(3)-current_span_z/2, center_pt(3)+current_span_z/2];
xlim(lims_x); ylim(lims_y); zlim(lims_z);

% 目盛りを20度刻みで綺麗に整える
xticks(ceil(lims_x(1)/20)*20 : 20 : floor(lims_x(2)/20)*20);
yticks(ceil(lims_y(1)/20)*20 : 20 : floor(lims_y(2)/20)*20);
zticks(ceil(lims_z(1)/20)*20 : 20 : floor(lims_z(2)/20)*20);

xlabel(axis_labels_3D{1}); ylabel(axis_labels_3D{2}); zlabel(axis_labels_3D{3});
set(gca, 'XDir', 'reverse'); set(gca, 'YDir', 'reverse');

% 平面枠の描画
if ~isempty(final_intersect_R), hP_R = plot3(final_intersect_R(:,1), final_intersect_R(:,2), final_intersect_R(:,3), '-', 'Color', col_Traj_R, 'LineWidth', 0.5); end
if ~isempty(final_intersect_L), hP_L = plot3(final_intersect_L(:,1), final_intersect_L(:,2), final_intersect_L(:,3), '-', 'Color', col_Traj_L, 'LineWidth', 0.5); end

% 軌道の描画
hT_R = plot3(orbit_R(:,1), orbit_R(:,2), orbit_R(:,3), '-', 'Color', col_Traj_R, 'LineWidth', 1.0);
plot3(orbit_R(1,1), orbit_R(1,2), orbit_R(1,3), 'o', 'MarkerEdgeColor', col_Traj_R, 'MarkerFaceColor', col_Traj_R, 'MarkerSize', 8);
plot3(orbit_R(mean_LO_index_R,1), orbit_R(mean_LO_index_R,2), orbit_R(mean_LO_index_R,3), 's', 'MarkerEdgeColor', col_Traj_R, 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.0);

hT_L = plot3(orbit_L(:,1), orbit_L(:,2), orbit_L(:,3), '-', 'Color', col_Traj_L, 'LineWidth', 1.0);
plot3(orbit_L(1,1), orbit_L(1,2), orbit_L(1,3), 'o', 'MarkerEdgeColor', col_Traj_L, 'MarkerFaceColor', col_Traj_L, 'MarkerSize', 8);
plot3(orbit_L(mean_LO_index_L,1), orbit_L(mean_LO_index_L,2), orbit_L(mean_LO_index_L,3), 's', 'MarkerEdgeColor', col_Traj_L, 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.0);

h_dummy_PlaneR = plot3(nan, nan, nan, '-', 'Color', col_Traj_R, 'LineWidth', 0.5);
h_dummy_PlaneL = plot3(nan, nan, nan, '-', 'Color', col_Traj_L, 'LineWidth', 0.5);
h_dummy_TD = plot3(nan, nan, nan, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
h_dummy_LO = plot3(nan, nan, nan, 'ks', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.0);

legend([hT_R, h_dummy_PlaneR, hT_L, h_dummy_PlaneL, h_dummy_TD, h_dummy_LO], ...
       {'Right Orbit (Non-paretic)', 'Right Plane Frame', 'Left Orbit (Paretic)', 'Left Plane Frame', 'TD (Touch Down)', 'LO (Lift Off)'}, ...
       'Location', 'bestoutside');

% --- [Fig 6] PCA Combined ---
fprintf('Plotting Fig 6: Combined PCA Projection...\n');
fig_6 = figure(6); clf(fig_6);
set(fig_6, 'Name', 'PCA Combined', 'NumberTitle', 'off');
movegui(fig_6, 'center'); 
hold on; grid on; axis equal;

% 健側 (Red)
hT_R = plot(score_R(:,1), score_R(:,2), '-', 'Color', col_Traj_R, 'LineWidth', 1.0);
plot(score_R(1,1), score_R(1,2), 'o', 'MarkerEdgeColor', col_Traj_R, 'MarkerFaceColor', col_Traj_R, 'MarkerSize', 8);
plot(score_R(mean_LO_index_R,1), score_R(mean_LO_index_R,2), 's', 'MarkerEdgeColor', col_Traj_R, 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.0);

% 患側 (Blue)
hT_L = plot(score_L(:,1), score_L(:,2), '-', 'Color', col_Traj_L, 'LineWidth', 1.0);
plot(score_L(1,1), score_L(1,2), 'o', 'MarkerEdgeColor', col_Traj_L, 'MarkerFaceColor', col_Traj_L, 'MarkerSize', 8);
plot(score_L(mean_LO_index_L,1), score_L(mean_LO_index_L,2), 's', 'MarkerEdgeColor', col_Traj_L, 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.0);

xlabel(sprintf('PC1 (R: %.1f%%, L: %.1f%%)', explained_R(1), explained_L(1)));
ylabel(sprintf('PC2 (R: %.1f%%, L: %.1f%%)', explained_R(2), explained_L(2)));

% 軸の余白と目盛りの調整
all_scores = [score_R(:,1:2); score_L(:,1:2)];
min_xy = min(all_scores); max_xy = max(all_scores);
margin_ratio = 0.1; span = max_xy - min_xy;
lims_x = [min_xy(1) - span(1)*margin_ratio, max_xy(1) + span(1)*margin_ratio];
lims_y = [min_xy(2) - span(2)*margin_ratio, max_xy(2) + span(2)*margin_ratio];
xlim(lims_x); ylim(lims_y);
xticks(floor(lims_x(1)/20)*20 : 20 : ceil(lims_x(2)/20)*20);
yticks(floor(lims_y(1)/20)*20 : 20 : ceil(lims_y(2)/20)*20);

%% ===================================================================
% 9. [Fig 7] 個別ストライドの統計プロット
% ===================================================================
% 全ストライドに対する個別のPCA結果を展開し、分散や平面の類似度を評価する。
n_strides_valid = min(num_strides_R, num_strides_L); 

res_Ratio_R = indiv_R.ratio(1:n_strides_valid);     
res_Ratio_L = indiv_L.ratio(1:n_strides_valid);
res_Area_R  = indiv_R.area(1:n_strides_valid);
res_Area_L  = indiv_L.area(1:n_strides_valid);
res_Planarity_R = indiv_R.planarity(1:n_strides_valid);
res_Planarity_L = indiv_L.planarity(1:n_strides_valid);

res_DotProd = zeros(n_strides_valid, 1);
res_Angle   = zeros(n_strides_valid, 1);
for k = 1:n_strides_valid
    n_vec_R = indiv_R.coeff(:, 3, k);
    n_vec_L = indiv_L.coeff(:, 3, k);
    dot_val = abs(dot(n_vec_R, n_vec_L));
    res_DotProd(k) = dot_val;
    res_Angle(k) = acos(min(dot_val, 1)) * 180/pi;
end

fprintf('Plotting Fig 7: Stride-by-Stride Statistics...\n');
fig_7 = figure(7); clf(fig_7);
set(fig_7, 'Name', 'Planar Metrics Stats', 'NumberTitle', 'off');
movegui(fig_7, 'center');

col_BoxLine = 'k'; col_Left = [0, 0, 0.8]; col_Right = [0.8, 0, 0]; col_Single = [0, 0, 0];
jit_width = 0.2;
calc_jitter = @(n) (rand(n,1)-0.5)*jit_width;
center_jitter = @(x) x - mean(x); 
x_L = 1 + center_jitter(calc_jitter(n_strides_valid));
x_R = 2 + center_jitter(calc_jitter(n_strides_valid));
x_S = 1 + center_jitter(calc_jitter(n_strides_valid));
calc_ylim = @(d) [min(d) - (max(d)-min(d))*0.2, max(d) + (max(d)-min(d))*0.2];

% 1. Planarity
subplot(2, 2, 1); hold on; grid on;
boxplot([res_Planarity_L, res_Planarity_R], 'Positions', [1, 2], 'Labels', {'Left', 'Right'}, 'Symbol', '', 'Colors', col_BoxLine);
scatter(x_L, res_Planarity_L, 40, col_Left, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(x_R, res_Planarity_R, 40, col_Right, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('Planarity (%)'); title('Planarity (PC1+PC2)');
xlim([0.5, 2.5]);
if max([res_Planarity_L; res_Planarity_R]) > min([res_Planarity_L; res_Planarity_R]), ylim(calc_ylim([res_Planarity_L; res_Planarity_R])); else, ylim([90, 100]); end

% 2. Ratio
subplot(2, 2, 2); hold on; grid on;
boxplot([res_Ratio_L, res_Ratio_R], 'Positions', [1, 2], 'Labels', {'Left', 'Right'}, 'Symbol', '', 'Colors', col_BoxLine);
scatter(x_L, res_Ratio_L, 40, col_Left, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(x_R, res_Ratio_R, 40, col_Right, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('Ratio (PV1 / PV2)'); title('Planarity Ratio');
xlim([0.5, 2.5]); 
if max([res_Ratio_L; res_Ratio_R]) > min([res_Ratio_L; res_Ratio_R]), ylim(calc_ylim([res_Ratio_L; res_Ratio_R])); end

% 3. Area
subplot(2, 2, 3); hold on; grid on;
boxplot([res_Area_L, res_Area_R], 'Positions', [1, 2], 'Labels', {'Left', 'Right'}, 'Symbol', '', 'Colors', col_BoxLine);
scatter(x_L, res_Area_L, 40, col_Left, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(x_R, res_Area_R, 40, col_Right, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('Area (deg^2)'); title('Geometric Area');
xlim([0.5, 2.5]);
if max([res_Area_L; res_Area_R]) > min([res_Area_L; res_Area_R]), ylim(calc_ylim([res_Area_L; res_Area_R])); end

% 4. Dot Product
subplot(2, 2, 4); hold on; grid on;
boxplot(res_DotProd, 'Positions', 1, 'Labels', {'R-L Similarity'}, 'Symbol', '', 'Colors', col_BoxLine);
scatter(x_S, res_DotProd, 40, col_Single, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('Dot Product'); title('Normal Vector Dot Product');
xlim([0.5, 1.5]);
if max(res_DotProd) > min(res_DotProd), ylim(calc_ylim(res_DotProd)); else, ylim([0.95, 1.05]); end

sgtitle('Stride-by-Stride Planar Law Metrics', 'FontSize', 16, 'FontWeight', 'bold');

%% ===================================================================
% 10. コンソールへの統計結果出力
% ===================================================================
val_Planarity_MeanOrbit_L = explained_L(1) + explained_L(2);
val_Planarity_MeanOrbit_R = explained_R(1) + explained_R(2);

if ~exist('ratio_L', 'var'), ratio_L = explained_L(1) / explained_L(2); end
if ~exist('ratio_R', 'var'), ratio_R = explained_R(1) / explained_R(2); end
if ~exist('area_L', 'var'),  area_L  = polyarea(score_L(:,1), score_L(:,2)); end
if ~exist('area_R', 'var'),  area_R  = polyarea(score_R(:,1), score_R(:,2)); end
if ~exist('normal_dot_product', 'var'), normal_dot_product = dot(coeff_R(:,3), coeff_L(:,3)); end

fprintf('\n');
fprintf('========================================================================================\n');
fprintf('                               Planar Law Analysis Report                               \n');
fprintf('========================================================================================\n');
fprintf('%-28s | %-16s | %-25s\n', 'Metric', 'Mean Orbit', 'Individual Strides');
fprintf('%-28s | %-16s | %-25s\n', '', '(Step 3 Result)', '(Mean +/- Std Dev)');
fprintf('----------------------------------------------------------------------------------------\n');

fprintf('Planarity (Left) [%%]        : %10.2f       |  %6.2f  +/- %5.2f\n', ...
    val_Planarity_MeanOrbit_L, mean(res_Planarity_L), std(res_Planarity_L));
fprintf('Planarity (Right) [%%]       : %10.2f       |  %6.2f  +/- %5.2f\n', ...
    val_Planarity_MeanOrbit_R, mean(res_Planarity_R), std(res_Planarity_R));
fprintf('----------------------------------------------------------------------------------------\n');

fprintf('Eigenvalue Ratio (Left)     : %10.2f       |  %6.2f  +/- %5.2f\n', ...
    ratio_L, mean(res_Ratio_L), std(res_Ratio_L));
fprintf('Eigenvalue Ratio (Right)    : %10.2f       |  %6.2f  +/- %5.2f\n', ...
    ratio_R, mean(res_Ratio_R), std(res_Ratio_R));
fprintf('----------------------------------------------------------------------------------------\n');

fprintf('Orbital Area (Left) [deg^2]  : %10.1f       |  %6.1f  +/- %5.1f\n', ...
    area_L, mean(res_Area_L), std(res_Area_L));
fprintf('Orbital Area (Right) [deg^2] : %10.1f       |  %6.1f  +/- %5.1f\n', ...
    area_R, mean(res_Area_R), std(res_Area_R));
fprintf('----------------------------------------------------------------------------------------\n');

fprintf('Symmetry Index              : %10.4f       |  %6.4f  +/- %5.4f\n', ...
    normal_dot_product, mean(res_DotProd), std(res_DotProd));
fprintf('========================================================================================\n');

%% ===================================================================
% 11. 論文用テーブルデータの出力 (Individual Stats)
% ===================================================================
% 平均軌道に対してではなく、個別のストライドごとにSVDを適用した際の特異値と累積寄与率を出力する。
vals_SV_R = zeros(num_strides_R, 3); vals_CP_R = zeros(num_strides_R, 3); 
for i = 1:num_strides_R
    orb_R_centered = all_normalized_R_orbits(:,:,i) - mean(all_normalized_R_orbits(:,:,i), 1);
    s_vals = svd(orb_R_centered, 'econ'); 
    vals_SV_R(i, :) = s_vals';
    s_sq = s_vals.^2; vals_CP_R(i, :) = (cumsum(s_sq) / sum(s_sq))';
end
mean_SV_R = mean(vals_SV_R, 1); std_SV_R = std(vals_SV_R, 0, 1);
mean_CP_R = mean(vals_CP_R, 1); std_CP_R = std(vals_CP_R, 0, 1);

vals_SV_L = zeros(num_strides_L, 3); vals_CP_L = zeros(num_strides_L, 3);
for i = 1:num_strides_L
    orb_L_centered = all_normalized_L_orbits(:,:,i) - mean(all_normalized_L_orbits(:,:,i), 1);
    s_vals = svd(orb_L_centered, 'econ');
    vals_SV_L(i, :) = s_vals';
    s_sq = s_vals.^2; vals_CP_L(i, :) = (cumsum(s_sq) / sum(s_sq))';
end
mean_SV_L = mean(vals_SV_L, 1); std_SV_L = std(vals_SV_L, 0, 1);
mean_CP_L = mean(vals_CP_L, 1); std_CP_L = std(vals_CP_L, 0, 1);

fprintf('\n======================================================================================================\n');
fprintf('      Additional Table Data: Singular Values & Cumulative Proportion (Individual Statistics)          \n');
fprintf('======================================================================================================\n');
fprintf('Subject: %s\n', SubjectName);
fprintf('Calculation: SVD performed on EACH stride independently.\n');
fprintf('Values shown as: Mean +/- S.D. over all valid strides.\n');
fprintf('------------------------------------------------------------------------------------------------------\n');
fprintf('Side   | Comp |      Singular Value (lambda)      |   Cumulative Proportion (Lambda)   \n');
fprintf('       |      |      (Mean +/- SD)                |      (Mean +/- SD)                 \n');
fprintf('------------------------------------------------------------------------------------------------------\n');
for k = 1:3
    fprintf('Right  |   %d  | %10.4f +/- %8.4f           | %10.4f +/- %8.4f\n', ...
        k, mean_SV_R(k), std_SV_R(k), mean_CP_R(k), std_CP_R(k));
end
fprintf('------------------------------------------------------------------------------------------------------\n');
for k = 1:3
    fprintf('Left   |   %d  | %10.4f +/- %8.4f           | %10.4f +/- %8.4f\n', ...
        k, mean_SV_L(k), std_SV_L(k), mean_CP_L(k), std_CP_L(k));
end
fprintf('======================================================================================================\n\n');

%% ===================================================================
% 12. グラフ画像の一括保存
% ===================================================================
if flg_graphSave
    fprintf('Saving figures using func_graphSave2...\n');
    size_Fig1 = [400, 250]; size_Fig2 = [400, 250]; size_Fig3 = [1100, 500];  
    size_Fig4 = [500, 800]; size_Fig5 = [700, 300]; size_Fig6 = [400, 300];   
    size_Fig7 = [1200, 800];  
    
    func_graphSave2(fig_1, tTextGraphName_Fig1_L,   flg_graphSave, 0, size_Fig1);
    func_graphSave2(fig_2, tTextGraphName_Fig2_R,   flg_graphSave, 0, size_Fig2);
    func_graphSave2(fig_3, tTextGraphName_Fig3_3D_Split, flg_graphSave, 0, size_Fig3);
    func_graphSave2(fig_4, tTextGraphName_Fig4_PCA_Split, flg_graphSave, 0, size_Fig4);
    func_graphSave2(fig_5, tTextGraphName_Fig5_3D_Comb, flg_graphSave, 0, size_Fig5);
    func_graphSave2(fig_6, tTextGraphName_Fig6_PCA_Comb, flg_graphSave, 0, size_Fig6);
    func_graphSave2(fig_7, tTextGraphName_Fig7_Stats, flg_graphSave, 0, size_Fig7);
end

%% ===================================================================
% 13. 解析データの保存
% ===================================================================
fprintf('Saving results to %s.mat...\n', outputMatName);

save(outputMatName, ...
    'SubjectName', 'samplingFrequency', 'data', 'rawData', ...
    'theta', 'centeredtheta', 'mean_posture', 'contact_flags', ...
    'theta_centered_interp', 'theta_interp', 'contact_flags_interp', ...
    'valid_strides_R', 'valid_strides_L', 'num_valid_R', ...
    'all_normalized_strides', 'all_normalized_strides_RAW', ...
    'averaged_stride', 'std_stride', ...
    'averaged_DS1', 'averaged_SS1', 'averaged_DS2', 'averaged_SS2', ...
    'all_norm_DS1', 'all_norm_SS1', 'all_norm_DS2', 'all_norm_SS2', ...
    'mean_phase_transition_pcts', ...
    'mean_orbit_R_cent', 'mean_orbit_L_cent', ...
    'mean_orbit_R_raw', 'mean_orbit_L_raw', ...
    'std_orbit_R_raw', 'std_orbit_L_raw', ...
    'coeff_R', 'score_R', 'explained_R', ...
    'coeff_L', 'score_L', 'explained_L', ...
    'mean_LO_index_R', 'mean_LO_index_L', ...
    'indiv_R', 'indiv_L', ... 
    'ratio_R', 'ratio_L', 'area_R', 'area_L', 'normal_dot_product', ...
    'res_Ratio_R', 'res_Ratio_L', 'res_Area_R', 'res_Area_L', ...
    'res_Planarity_R', 'res_Planarity_L', ... 
    'res_DotProd', 'res_Angle', ...
    'axis_labels_3D', 'cols_R_Leg', 'cols_L_Leg');

fprintf('=== ALL PROCESS COMPLETED ===\n');