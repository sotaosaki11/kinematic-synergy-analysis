close all; 
clear;     

% 機器設定ファイルの読み込み
if exist('make_status_bertec', 'file'), make_status_bertec; end
if exist('make_status_mocap_29', 'file'), make_status_mocap_29; end

%% =========================================================================
% 【プログラム概要: 平面則 (Planar Law) 解析】
% 大腿・下腿・足部の3体節の仰角データに対して主成分分析(PCA)および特異値分解(SVD)
% を行い、歩行中の姿勢が形成する「共変動面（Covariation Plane）」を算出・評価する。
% 平面性の度合い（Planarity）、形状（Ratio, Area）、および左右の対称性（内積）
% などの指標を計算し、2D/3Dグラフによる視覚化と統計データの出力を行う。
% =========================================================================

%% ===================================================================
% 1. ファイル名とパスの設定
% ===================================================================
tTextFilePath  = './'; 

% 解析対象の被験者名
tTextFileName   = 'IBA'; 

% 読み込むファイル（segment_control.m の出力結果）
tTextLoadDataName = [tTextFileName, '_12_step2_segmented'];

% 保存するデータファイル名
tTextSaveMatName  = [tTextFileName, '_12_step3_planerlaw']; 

% グラフ画像ファイルの保存名
tTextGraphName_Fig1_L         = [tTextFileName, '_Fig1_2D_Angles_Left']; 
tTextGraphName_Fig2_R         = [tTextFileName, '_Fig2_2D_Angles_Right']; 
tTextGraphName_Fig3_3D_Split  = [tTextFileName, '_Fig3_3D_Orbits_Split']; 
tTextGraphName_Fig4_PCA_Split = [tTextFileName, '_Fig4_PCA_Split'];
tTextGraphName_Fig5_3D_Comb   = [tTextFileName, '_Fig5_3D_Orbits_Combined']; 
tTextGraphName_Fig6_PCA_Comb  = [tTextFileName, '_Fig6_PCA_Combined'];       
tTextGraphName_Fig7_Stats     = [tTextFileName, '_Fig7_Planar_Stats']; 

flg_graphSave = 1; % グラフ保存フラグ

%% ===================================================================
% 2. 解析・描画パラメータの設定
% ===================================================================
N_RESAMPLE_POINTS = 101;  % 1歩行周期を101点(0~100%)にリサンプリング

% 3Dグラフのカメラ視点 (Azimuth, Elevation)
param_ViewAngle = [160, 10]; 

% 解析対象の体節インデックス（preprocess.mで定義したthetaの列番号）
% [大腿(Thigh), 下腿(Shank), 足部(Foot)]
cols_R_Leg = [2, 4, 6]; % 右脚
cols_L_Leg = [3, 5, 7]; % 左脚
axis_labels_3D = {'Thigh Angle (deg)', 'Shank Angle (deg)', 'Foot Angle (deg)'};

% グラフ描画用のカラー設定
color_Plane_R = [1, 0.7, 0.7]; 
color_Plane_L = [0.7, 0.7, 1]; 
color_Shank = [0, 0.6, 0];   color_Thigh = [0, 0, 0.8];   color_Foot  = [0.8, 0, 0];
color_Shank_Shade = [0.7, 1, 0.7]; color_Thigh_Shade = [0.7, 0.7, 1]; color_Foot_Shade  = [1, 0.7, 0.7]; 
color_Stance = [0, 0, 0.8]; color_Swing  = [0, 0.7, 0]; 

% 左右比較グラフ（Fig 5, 6）用のカラー設定
col_Traj_R = [0.8, 0, 0];      % 右脚軌道（赤）
col_Traj_L = [0, 0, 0.8];      % 左脚軌道（青）
col_Plane_R_Face = col_Traj_R; 
col_Plane_L_Face = col_Traj_L; 

% Fig 1, 2の縦横比（Aspect Ratio）
param_Fig1_AspectRatio = [1.6, 1, 1]; 

%% ===================================================================
% 3. データの読込とバリデーション
% ===================================================================
fprintf('Loading data from %s.mat...\n', tTextLoadDataName);
load(tTextLoadDataName); 

% segment_control.m で選別された有効なストライド（valid_strides）が存在するか確認
if ~exist('valid_strides_R', 'var') || ~exist('valid_strides_L', 'var')
    error('valid_strides_R / L がMATファイルに見つかりません。segment_control.mを更新して実行してください。');
end

num_strides_R = length(valid_strides_R);
num_strides_L = length(valid_strides_L);
fprintf('Using Pre-Segmented Strides:\n  Right: %d strides\n  Left : %d strides\n', num_strides_R, num_strides_L);

%% ===================================================================
% 4. PCA用の中心化軌道（Centered Orbits）の算出
% ===================================================================
% 平面則の抽出には、時間平均を引いた「角度の変動成分（Centered Theta）」を用いる。
fprintf('Calculating CENTERED orbits for 3D/PCA from Valid Strides...\n');
resampled_time_vector = linspace(0, 1, N_RESAMPLE_POINTS)';

% --- 右脚 (Right Leg) の処理 ---
all_normalized_R_orbits = zeros(N_RESAMPLE_POINTS, 3, num_strides_R);
LO_pct_list_R = zeros(num_strides_R, 1); 

for i = 1:num_strides_R
    st = valid_strides_R(i);
    n_len = size(st.raw_stride, 1);
    all_normalized_R_orbits(:,:,i) = interp1(linspace(0,1,n_len), st.raw_stride(:, cols_R_Leg), resampled_time_vector, 'pchip');
    LO_pct_list_R(i) = st.pct_LO;
end

% 全ストライドの平均軌道を算出し、さらにその平均姿勢で中心化する
mean_orbit_R = mean(all_normalized_R_orbits, 3); 
std_orbit_R  = std(all_normalized_R_orbits, 0, 3); 
mean_orbit_R_centered = mean_orbit_R - mean(mean_orbit_R, 1);

% 平均軌道に対するPCA（主成分分析）の実行
% coeff_R: 主成分ベクトル（3列目が平面の法線ベクトルとなる）
% score_R: 平面上に射影された軌道の座標
% explained_R: 各主成分の寄与率（分散の割合）
[coeff_R, score_R, ~, ~, explained_R] = pca(mean_orbit_R_centered);

mean_normalized_LO_R = mean(LO_pct_list_R); 
mean_LO_index_R = round(mean_normalized_LO_R * (N_RESAMPLE_POINTS - 1)) + 1; 

% --- 左脚 (Left Leg) の処理 ---
all_normalized_L_orbits = zeros(N_RESAMPLE_POINTS, 3, num_strides_L);
LO_pct_list_L = zeros(num_strides_L, 1);

for i = 1:num_strides_L
    st = valid_strides_L(i);
    n_len = size(st.raw_stride, 1);
    all_normalized_L_orbits(:,:,i) = interp1(linspace(0,1,n_len), st.raw_stride(:, cols_L_Leg), resampled_time_vector, 'pchip');
    LO_pct_list_L(i) = st.pct_LO;
end

mean_orbit_L = mean(all_normalized_L_orbits, 3); 
std_orbit_L  = std(all_normalized_L_orbits, 0, 3); 
mean_orbit_L_centered = mean_orbit_L - mean(mean_orbit_L, 1);

[coeff_L, score_L, ~, ~, explained_L] = pca(mean_orbit_L_centered);

mean_normalized_LO_L = mean(LO_pct_list_L); 
mean_LO_index_L = round(mean_normalized_LO_L * (N_RESAMPLE_POINTS - 1)) + 1; 

%% ===================================================================
% 5. 2Dプロット用の絶対角度軌道（Absolute Orbits）の算出
% ===================================================================
% 平面分析には使わないが、Fig 1, 2の時系列プロット用に実際の角度（中心化なし）を用意する。
fprintf('Calculating RAW (Absolute) orbits for Fig 1 & 2...\n');

% 右脚
all_normalized_R_orbits_RAW = zeros(N_RESAMPLE_POINTS, 3, num_strides_R);
for i = 1:num_strides_R
    st = valid_strides_R(i);
    n_len = size(st.raw_stride_theta, 1);
    all_normalized_R_orbits_RAW(:,:,i) = interp1(linspace(0,1,n_len), st.raw_stride_theta(:, cols_R_Leg), resampled_time_vector, 'pchip');
end
mean_orbit_R_RAW = mean(all_normalized_R_orbits_RAW, 3);
std_orbit_R_RAW  = std(all_normalized_R_orbits_RAW, 0, 3);

% 左脚
all_normalized_L_orbits_RAW = zeros(N_RESAMPLE_POINTS, 3, num_strides_L);
for i = 1:num_strides_L
    st = valid_strides_L(i);
    n_len = size(st.raw_stride_theta, 1);
    all_normalized_L_orbits_RAW(:,:,i) = interp1(linspace(0,1,n_len), st.raw_stride_theta(:, cols_L_Leg), resampled_time_vector, 'pchip');
end
mean_orbit_L_RAW = mean(all_normalized_L_orbits_RAW, 3);
std_orbit_L_RAW  = std(all_normalized_L_orbits_RAW, 0, 3);

%% ===================================================================
% 6. 各種グラフの生成と描画
% ===================================================================

% --- [Fig 1] 2D 仰角の時系列変化 (Left Leg) ---
fprintf('Plotting Fig 1: 2D Time Series (Left)...\n');
fig_1 = figure(1); clf(fig_1);
set(fig_1, 'Name', '2D Angles Left', 'NumberTitle', 'off');
movegui(fig_1, 'west'); 

time_axis = linspace(0, 100, N_RESAMPLE_POINTS)'; 
plot_x_axis_fill = [time_axis; flipud(time_axis)]; 

hold on; grid on;
fill(plot_x_axis_fill, [mean_orbit_L_RAW(:,1)+std_orbit_L_RAW(:,1); flipud(mean_orbit_L_RAW(:,1)-std_orbit_L_RAW(:,1))], color_Thigh_Shade, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
fill(plot_x_axis_fill, [mean_orbit_L_RAW(:,2)+std_orbit_L_RAW(:,2); flipud(mean_orbit_L_RAW(:,2)-std_orbit_L_RAW(:,2))], color_Shank_Shade, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
fill(plot_x_axis_fill, [mean_orbit_L_RAW(:,3)+std_orbit_L_RAW(:,3); flipud(mean_orbit_L_RAW(:,3)-std_orbit_L_RAW(:,3))], color_Foot_Shade, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hT_L = plot(time_axis, mean_orbit_L_RAW(:,1), 'Color', color_Thigh, 'LineWidth', 2);
hS_L = plot(time_axis, mean_orbit_L_RAW(:,2), 'Color', color_Shank, 'LineWidth', 2);
hF_L = plot(time_axis, mean_orbit_L_RAW(:,3), 'Color', color_Foot, 'LineWidth', 2);
ylabel('Angle (deg)'); xlabel('Gait Cycle (%)'); xlim([0 100]); 
pbaspect(param_Fig1_AspectRatio); 
legend([hT_L, hS_L, hF_L], axis_labels_3D, 'Location', 'best'); 

% --- [Fig 2] 2D 仰角の時系列変化 (Right Leg) ---
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
pbaspect(param_Fig1_AspectRatio); 
legend([hT_R, hS_R, hF_R], axis_labels_3D, 'Location', 'best'); 


% --- [Fig 3] 3D 角度軌道と共変動面 (左右分割) ---
% 3体節の角度を3次元空間にプロットし、PCAによって得られた平面（Planar Law）を描画する。
fprintf('Plotting Fig 3: 3D Orbit Split...\n');
fig_3 = figure(3); clf(fig_3);
set(fig_3, 'Name', '3D Split', 'NumberTitle', 'off');
movegui(fig_3, 'center'); 

pad_ratio = 1.1; 
pos_L = [0.07, 0.15, 0.38, 0.75]; 
pos_R = [0.50, 0.15, 0.38, 0.75]; 

% 左脚
axL = subplot('Position', pos_L); hold on; grid on; axis equal; view(param_ViewAngle);
orbit_L = mean_orbit_L_centered; nL = coeff_L(:, 3); % 第3主成分が法線ベクトル

min_L = min(orbit_L); max_L = max(orbit_L); center_L = (min_L + max_L) / 2;
max_span = max(max_L - min_L) * pad_ratio; half_span = max_span / 2;
lims_x = [center_L(1)-half_span, center_L(1)+half_span];
lims_y = [center_L(2)-half_span, center_L(2)+half_span];
lims_z = [center_L(3)-half_span, center_L(3)+half_span];
[gx, gy] = meshgrid(lims_x, lims_y); gz = (-nL(1)*gx - nL(2)*gy)/nL(3); % 平面の方程式
surf(gx, gy, gz, 'FaceColor', color_Plane_L, 'FaceAlpha', 0.4, 'EdgeColor', 'none');

plot3(orbit_L(1:mean_LO_index_L,1), orbit_L(1:mean_LO_index_L,2), orbit_L(1:mean_LO_index_L,3), '-', 'Color', color_Stance, 'LineWidth', 2);
plot3(orbit_L(mean_LO_index_L:end,1), orbit_L(mean_LO_index_L:end,2), orbit_L(mean_LO_index_L:end,3), '-', 'Color', color_Swing, 'LineWidth', 2);
plot3(orbit_L(1,1), orbit_L(1,2), orbit_L(1,3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
plot3(orbit_L(mean_LO_index_L,1), orbit_L(mean_LO_index_L,2), orbit_L(mean_LO_index_L,3), 'ks', 'MarkerFaceColor', 'w', 'MarkerSize', 8);

xlim(lims_x); ylim(lims_y); zlim(lims_z);
xlabel(axis_labels_3D{1}); ylabel(axis_labels_3D{2}); zlabel(axis_labels_3D{3});
set(gca, 'XDir', 'reverse', 'YDir', 'reverse'); 

% 右脚
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
set(axR, 'Position', pos_R);
xlim(lims_x); ylim(lims_y); zlim(lims_z);
xlabel(axis_labels_3D{1}); ylabel(axis_labels_3D{2}); zlabel(axis_labels_3D{3});
set(gca, 'XDir', 'reverse', 'YDir', 'reverse'); 


% --- [Fig 4] PCA 平面上への射影 (左右分割) ---
% 抽出された平面（PC1-PC2平面）を正面から見た際のループ軌道を描画する。
fprintf('Plotting Fig 4: PCA Projection Split...\n');
fig_4 = figure(4); clf(fig_4);
set(fig_4, 'Name', 'PCA Split', 'NumberTitle', 'off');
movegui(fig_4, 'center'); 

pos_L = [0.15, 0.55, 0.70, 0.40]; 
pos_R = [0.15, 0.08, 0.70, 0.40]; 

% 左脚
axL = subplot('Position', pos_L); hold on; grid on; axis equal;
plot(score_L(:,1), score_L(:,2), 'k-');
plot(score_L(1:mean_LO_index_L,1), score_L(1:mean_LO_index_L,2), '-', 'Color', color_Stance, 'LineWidth', 2);
plot(score_L(mean_LO_index_L:end,1), score_L(mean_LO_index_L:end,2), '-', 'Color', color_Swing, 'LineWidth', 2);
plot(score_L(1,1), score_L(1,2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
plot(score_L(mean_LO_index_L,1), score_L(mean_LO_index_L,2), 'ks', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
xlabel(sprintf('PC1 (%.1f %%)', explained_L(1))); ylabel(sprintf('PC2 (%.1f %%)', explained_L(2))); 

% 右脚
axR = subplot('Position', pos_R); hold on; grid on; axis equal;
plot(score_R(:,1), score_R(:,2), 'k-');
hSt = plot(score_R(1:mean_LO_index_R,1), score_R(1:mean_LO_index_R,2), '-', 'Color', color_Stance, 'LineWidth', 2);
hSw = plot(score_R(mean_LO_index_R:end,1), score_R(mean_LO_index_R:end,2), '-', 'Color', color_Swing, 'LineWidth', 2);
hTD = plot(score_R(1,1), score_R(1,2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
hLO = plot(score_R(mean_LO_index_R,1), score_R(mean_LO_index_R,2), 'ks', 'MarkerFaceColor', 'w', 'MarkerSize', 8);

legend([hSt, hSw, hTD, hLO], {'Stance', 'Swing', 'TD', 'LO'}, 'Location', 'best');
set(axR, 'Position', pos_R); 
xlabel(sprintf('PC1 (%.1f %%)', explained_R(1))); ylabel(sprintf('PC2 (%.1f %%)', explained_R(2))); 


% --- [Fig 5] 3D 軌道と共変動面の左右比較 (Combined) ---
% 左右の平面がどの程度ずれているか（対称性）を視覚的に比較するため、同じ空間にプロットする。
% 平面を四角形の枠線として描画するための幾何学的な最適化処理を含む。
fprintf('Plotting Fig 5: Combined 3D Orbits...\n');
fig_5 = figure(5); clf(fig_5);
set(fig_5, 'Name', '3D Combined', 'NumberTitle', 'off');
movegui(fig_5, 'center'); 

hold on; grid on; axis equal; view(param_ViewAngle);

orbit_R = mean_orbit_R_centered; 
orbit_L = mean_orbit_L_centered; 

% 両脚データを包含する基本ボックス範囲を取得
all_pts = [orbit_R; orbit_L];
min_pt = min(all_pts); max_pt = max(all_pts);
center_pt = (min_pt + max_pt) / 2;
base_span = (max_pt - min_pt) * 1.1; 

max_axis_span = max(base_span);
current_span_x = max_axis_span;
current_span_y = max_axis_span;
current_span_z = max_axis_span;

% [幾何最適化] 平面が六角形等で途切れないよう、四角形になるまで表示境界を拡張する
fprintf('  Optimizing box size to ensure quadrangular plane frames...\n');
max_iter = 50; 
expand_step = 1.05; 

legs_data = struct('coeff', {coeff_R, coeff_L}, 'col', {col_Traj_R, col_Traj_L});
final_intersect_R = [];
final_intersect_L = [];

for iter = 1:max_iter
    half_x = current_span_x / 2; half_y = current_span_y / 2; half_z = current_span_z / 2;
    box_min = center_pt - [half_x, half_y, half_z];
    box_max = center_pt + [half_x, half_y, half_z];
    
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
            vec = p2 - p1;
            denom = dot(n, vec);
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
        final_intersect_R = current_intersects{1};
        final_intersect_L = current_intersects{2};
        fprintf('    Converged at iteration %d. Frame is quadrangular.\n', iter);
        break;
    else
        current_span_y = current_span_y * expand_step;
    end
end

% 描画範囲と目盛りの設定
lims_x = [center_pt(1)-current_span_x/2, center_pt(1)+current_span_x/2];
lims_y = [center_pt(2)-current_span_y/2, center_pt(2)+current_span_y/2];
lims_z = [center_pt(3)-current_span_z/2, center_pt(3)+current_span_z/2];
xlim(lims_x); ylim(lims_y); zlim(lims_z);

xticks(ceil(lims_x(1)/20)*20 : 20 : floor(lims_x(2)/20)*20);
yticks(ceil(lims_y(1)/20)*20 : 20 : floor(lims_y(2)/20)*20);
zticks(ceil(lims_z(1)/20)*20 : 20 : floor(lims_z(2)/20)*20);

xlabel(axis_labels_3D{1}); ylabel(axis_labels_3D{2}); zlabel(axis_labels_3D{3});
set(gca, 'XDir', 'reverse'); set(gca, 'YDir', 'reverse');

% 平面枠の描画
if ~isempty(final_intersect_R)
    hP_R = plot3(final_intersect_R(:,1), final_intersect_R(:,2), final_intersect_R(:,3), '-', 'Color', col_Traj_R, 'LineWidth', 0.5);
end
if ~isempty(final_intersect_L)
    hP_L = plot3(final_intersect_L(:,1), final_intersect_L(:,2), final_intersect_L(:,3), '-', 'Color', col_Traj_L, 'LineWidth', 0.5);
end

% 軌道の描画
hT_R = plot3(orbit_R(:,1), orbit_R(:,2), orbit_R(:,3), '-', 'Color', col_Traj_R, 'LineWidth', 1.0);
plot3(orbit_R(1,1), orbit_R(1,2), orbit_R(1,3), 'o', 'MarkerEdgeColor', col_Traj_R, 'MarkerFaceColor', col_Traj_R, 'MarkerSize', 8);
plot3(orbit_R(mean_LO_index_R,1), orbit_R(mean_LO_index_R,2), orbit_R(mean_LO_index_R,3), 's', 'MarkerEdgeColor', col_Traj_R, 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.0);

hT_L = plot3(orbit_L(:,1), orbit_L(:,2), orbit_L(:,3), '-', 'Color', col_Traj_L, 'LineWidth', 1.0);
plot3(orbit_L(1,1), orbit_L(1,2), orbit_L(1,3), 'o', 'MarkerEdgeColor', col_Traj_L, 'MarkerFaceColor', col_Traj_L, 'MarkerSize', 8);
plot3(orbit_L(mean_LO_index_L,1), orbit_L(mean_LO_index_L,2), orbit_L(mean_LO_index_L,3), 's', 'MarkerEdgeColor', col_Traj_L, 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.0);

% 左右平面の法線ベクトル同士の内積（歩行の非対称性の指標となる）
nR = coeff_R(:, 3); nL = coeff_L(:, 3);
normal_dot_product = dot(nR, nL);

h_dummy_PlaneR = plot3(nan, nan, nan, '-', 'Color', col_Traj_R, 'LineWidth', 0.5);
h_dummy_PlaneL = plot3(nan, nan, nan, '-', 'Color', col_Traj_L, 'LineWidth', 0.5);
h_dummy_TD = plot3(nan, nan, nan, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
h_dummy_LO = plot3(nan, nan, nan, 'ks', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.0);

legend([hT_R, h_dummy_PlaneR, hT_L, h_dummy_PlaneL, h_dummy_TD, h_dummy_LO], ...
       {'Right Orbit', 'Right Plane Frame', 'Left Orbit', 'Left Plane Frame', 'TD (Touch Down)', 'LO (Lift Off)'}, ...
       'Location', 'bestoutside');


% --- [Fig 6] PCA 平面上の左右軌道比較 (Combined) ---
fprintf('Plotting Fig 6: Combined PCA Projection...\n');
fig_6 = figure(6); clf(fig_6);
set(fig_6, 'Name', 'PCA Combined', 'NumberTitle', 'off');
movegui(fig_6, 'center'); 
hold on; grid on; axis equal;

hT_R = plot(score_R(:,1), score_R(:,2), '-', 'Color', col_Traj_R, 'LineWidth', 1.0);
plot(score_R(1,1), score_R(1,2), 'o', 'MarkerEdgeColor', col_Traj_R, 'MarkerFaceColor', col_Traj_R, 'MarkerSize', 8);
plot(score_R(mean_LO_index_R,1), score_R(mean_LO_index_R,2), 's', 'MarkerEdgeColor', col_Traj_R, 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.0);

hT_L = plot(score_L(:,1), score_L(:,2), '-', 'Color', col_Traj_L, 'LineWidth', 1.0);
plot(score_L(1,1), score_L(1,2), 'o', 'MarkerEdgeColor', col_Traj_L, 'MarkerFaceColor', col_Traj_L, 'MarkerSize', 8);
plot(score_L(mean_LO_index_L,1), score_L(mean_LO_index_L,2), 's', 'MarkerEdgeColor', col_Traj_L, 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.0);

xlim_curr = xlim; ylim_curr = ylim;
xticks(ceil(xlim_curr(1)/20)*20 : 20 : floor(xlim_curr(2)/20)*20);
yticks(ceil(ylim_curr(1)/20)*20 : 20 : floor(ylim_curr(2)/20)*20);

xlabel(sprintf('PC1 (R: %.1f%%, L: %.1f%%)', explained_R(1), explained_L(1)));
ylabel(sprintf('PC2 (R: %.1f%%, L: %.1f%%)', explained_R(2), explained_L(2)));


%% ===================================================================
% 7. 平均軌道に対する Planar Law 指標の計算
% ===================================================================
% Ratio: 第1主成分と第2主成分の比率（ループの形状・幅を示す）
% Area:  ループが囲む面積（可動域の大きさを示す）
ratio_R = explained_R(1) / explained_R(2);
area_R  = polyarea(score_R(:,1), score_R(:,2));
ratio_L = explained_L(1) / explained_L(2);
area_L  = polyarea(score_L(:,1), score_L(:,2));

fprintf('Average Orbit Metrics:\n');
fprintf('  Right Ratio: %.2f, Area: %.2f\n', ratio_R, area_R);
fprintf('  Left  Ratio: %.2f, Area: %.2f\n', ratio_L, area_L);

%% ===================================================================
% 8. [Fig 7] 全ストライドに対する個別のPlanar Law統計処理
% ===================================================================
% 全試行の平均軌道だけでなく、1ストライドごとに個別にPCAを行い、
% 得られた各種指標（Planarity, Ratio, Area, 内積）のばらつきを評価する。
fprintf('Calculating Planar Law metrics for EACH valid stride...\n');

n_strides_valid = min(num_strides_R, num_strides_L); 

res_Ratio_R = zeros(n_strides_valid, 1); res_Area_R  = zeros(n_strides_valid, 1); res_Planarity_R = zeros(n_strides_valid, 1);
res_Ratio_L = zeros(n_strides_valid, 1); res_Area_L  = zeros(n_strides_valid, 1); res_Planarity_L = zeros(n_strides_valid, 1);
res_DotProd = zeros(n_strides_valid, 1); res_Angle   = zeros(n_strides_valid, 1); 

for i = 1:n_strides_valid
    % 右脚の個別処理
    orbit_R_i = all_normalized_R_orbits(:,:,i); 
    orbit_R_i_centered = orbit_R_i - mean(orbit_R_i, 1);
    [coeff_R_i, score_R_i, ~, ~, expl_R_i] = pca(orbit_R_i_centered);
    
    res_Ratio_R(i) = expl_R_i(1) / expl_R_i(2);
    res_Area_R(i)  = polyarea(score_R_i(:,1), score_R_i(:,2));
    res_Planarity_R(i) = expl_R_i(1) + expl_R_i(2); % 第1+第2主成分の累積寄与率（平面性）
    normal_R_i     = coeff_R_i(:,3);

    % 左脚の個別処理
    orbit_L_i = all_normalized_L_orbits(:,:,i);
    orbit_L_i_centered = orbit_L_i - mean(orbit_L_i, 1);
    [coeff_L_i, score_L_i, ~, ~, expl_L_i] = pca(orbit_L_i_centered);
    
    res_Ratio_L(i) = expl_L_i(1) / expl_L_i(2);
    res_Area_L(i)  = polyarea(score_L_i(:,1), score_L_i(:,2));
    res_Planarity_L(i) = expl_L_i(1) + expl_L_i(2); 
    normal_L_i     = coeff_L_i(:,3);

    % 左右の平面の類似度（内積）
    dot_val = abs(dot(normal_R_i, normal_L_i));
    res_DotProd(i) = dot_val;
    res_Angle(i) = acos(min(dot_val, 1)) * 180/pi;
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

% 1. Planarity (PC1+PC2)
subplot(2, 2, 1); hold on; grid on;
boxplot([res_Planarity_L, res_Planarity_R], 'Positions', [1, 2], 'Labels', {'Left', 'Right'}, 'Symbol', '', 'Colors', col_BoxLine);
scatter(x_L, res_Planarity_L, 40, col_Left, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(x_R, res_Planarity_R, 40, col_Right, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('Planarity (%)'); 
title('Planarity (PC1+PC2)'); 
xlim([0.5, 2.5]);
if max([res_Planarity_L; res_Planarity_R]) > min([res_Planarity_L; res_Planarity_R])
    ylim(calc_ylim([res_Planarity_L; res_Planarity_R])); 
else
    ylim([90, 100]); 
end

% 2. Ratio
subplot(2, 2, 2); hold on; grid on;
boxplot([res_Ratio_L, res_Ratio_R], 'Positions', [1, 2], 'Labels', {'Left', 'Right'}, 'Symbol', '', 'Colors', col_BoxLine);
scatter(x_L, res_Ratio_L, 40, col_Left, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(x_R, res_Ratio_R, 40, col_Right, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('Ratio (PV1 / PV2)'); 
title('Planarity Ratio'); 
xlim([0.5, 2.5]); 
if max([res_Ratio_L; res_Ratio_R]) > min([res_Ratio_L; res_Ratio_R]), ylim(calc_ylim([res_Ratio_L; res_Ratio_R])); end

% 3. Area
subplot(2, 2, 3); hold on; grid on;
boxplot([res_Area_L, res_Area_R], 'Positions', [1, 2], 'Labels', {'Left', 'Right'}, 'Symbol', '', 'Colors', col_BoxLine);
scatter(x_L, res_Area_L, 40, col_Left, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(x_R, res_Area_R, 40, col_Right, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('Area (deg^2)'); 
title('Geometric Area'); 
xlim([0.5, 2.5]);
if max([res_Area_L; res_Area_R]) > min([res_Area_L; res_Area_R]), ylim(calc_ylim([res_Area_L; res_Area_R])); end

% 4. Dot Product
subplot(2, 2, 4); hold on; grid on;
boxplot(res_DotProd, 'Positions', 1, 'Labels', {'R-L Similarity'}, 'Symbol', '', 'Colors', col_BoxLine);
scatter(x_S, res_DotProd, 40, col_Single, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('Dot Product'); 
title('Normal Vector Dot Product'); 
xlim([0.5, 1.5]);
if max(res_DotProd) > min(res_DotProd), ylim(calc_ylim(res_DotProd)); else, ylim([0.95, 1.05]); end

sgtitle('Stride-by-Stride Planar Law Metrics', 'FontSize', 16, 'FontWeight', 'bold'); 

%% ===================================================================
% 9. コンソールへの統計結果出力
% ===================================================================
val_Planarity_MeanOrbit_L = explained_L(1) + explained_L(2);
val_Planarity_MeanOrbit_R = explained_R(1) + explained_R(2);

fprintf('\n');
fprintf('========================================================================================\n');
fprintf('                               Planar Law Analysis Report                               \n');
fprintf('========================================================================================\n');
fprintf('%-28s | %-16s | %-25s\n', 'Metric', 'Mean Orbit', 'Individual Strides');
fprintf('%-28s | %-16s | %-25s\n', '', '(Step 5 Result)', '(Mean +/- Std Dev)');
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
% 10. 論文出力用：個別ストライドごとのSVD統計処理
% ===================================================================
% 論文で提示するための、特異値(Singular Value)と累積寄与率(Cumulative Proportion)を計算する。

% --- 右脚 (Right Leg) ---
vals_SV_R = zeros(num_strides_R, 3); 
vals_CP_R = zeros(num_strides_R, 3); 

for i = 1:num_strides_R
    orb_R = all_normalized_R_orbits(:,:,i);
    orb_R_centered = orb_R - mean(orb_R, 1);
    s_vals = svd(orb_R_centered, 'econ'); 
    vals_SV_R(i, :) = s_vals';
    s_sq = s_vals.^2;
    vals_CP_R(i, :) = (cumsum(s_sq) / sum(s_sq))';
end

mean_SV_R = mean(vals_SV_R, 1); std_SV_R = std(vals_SV_R, 0, 1);
mean_CP_R = mean(vals_CP_R, 1); std_CP_R = std(vals_CP_R, 0, 1);

% --- 左脚 (Left Leg) ---
vals_SV_L = zeros(num_strides_L, 3);
vals_CP_L = zeros(num_strides_L, 3);

for i = 1:num_strides_L
    orb_L = all_normalized_L_orbits(:,:,i);
    orb_L_centered = orb_L - mean(orb_L, 1);
    s_vals = svd(orb_L_centered, 'econ');
    vals_SV_L(i, :) = s_vals';
    s_sq = s_vals.^2;
    vals_CP_L(i, :) = (cumsum(s_sq) / sum(s_sq))';
end

mean_SV_L = mean(vals_SV_L, 1); std_SV_L = std(vals_SV_L, 0, 1);
mean_CP_L = mean(vals_CP_L, 1); std_CP_L = std(vals_CP_L, 0, 1);

% --- 結果出力 ---
fprintf('\n');
fprintf('======================================================================================================\n');
fprintf('      Additional Table Data: Singular Values & Cumulative Proportion (Individual Statistics)          \n');
fprintf('======================================================================================================\n');
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
fprintf('======================================================================================================\n');
fprintf('\n');

%% ===================================================================
% 11. データの保存
% ===================================================================
if ~exist('meanposture', 'var')
    if exist('theta', 'var'), meanposture = mean(theta, 1); else, meanposture = []; end
end
save(tTextSaveMatName, ...
    'data', 'theta', 'centeredtheta', 'meanposture', ...
    'averaged_stride', 'std_stride', 'all_normalized_strides', ...
    'averaged_DS1', 'averaged_SS1', 'averaged_DS2', 'averaged_SS2', ...
    'all_norm_DS1', 'all_norm_SS1', 'all_norm_DS2', 'all_norm_SS2', ...
    'mean_orbit_R', 'std_orbit_R', 'all_normalized_R_orbits', ...
    'mean_orbit_L', 'std_orbit_L', 'all_normalized_L_orbits', ...
    'mean_orbit_R_centered', 'mean_orbit_L_centered', ...
    'mean_orbit_R_RAW', 'std_orbit_R_RAW', ... 
    'mean_orbit_L_RAW', 'std_orbit_L_RAW', ... 
    'axis_labels_3D', 'cols_R_Leg', 'cols_L_Leg', ...
    'coeff_R', 'explained_R', 'score_R', ... 
    'coeff_L', 'explained_L', 'score_L', ... 
    'ratio_R', 'ratio_L', 'area_R', 'area_L', ... 
    'mean_LO_index_R', 'mean_LO_index_L', ...
    'mean_normalized_LO_R', 'mean_normalized_LO_L', ...
    'normal_dot_product', ...
    'res_Ratio_R', 'res_Area_R', 'res_Ratio_L', 'res_Area_L', ...
    'res_Planarity_R', 'res_Planarity_L', ... 
    'res_DotProd', 'res_Angle');


%% ===================================================================
% 12. グラフ画像の一括保存
% ===================================================================
if flg_graphSave
    fprintf('Saving figures using func_graphSave2...\n');
    
    % 保存する画像サイズの指定 [width, height]
    size_Fig1 = [400, 250];   
    size_Fig2 = [400, 250];   
    size_Fig3 = [1100, 500];  
    size_Fig4 = [500, 800];   
    size_Fig5 = [700, 300];   
    size_Fig6 = [400, 300];   
    size_Fig7 = [1200, 800];  
    
    % カスタム関数で一括保存（第4引数の0はサイズ指定を有効にするフラグ）
    func_graphSave2(fig_1, tTextGraphName_Fig1_L,   flg_graphSave, 0, size_Fig1);
    func_graphSave2(fig_2, tTextGraphName_Fig2_R,   flg_graphSave, 0, size_Fig2);
    func_graphSave2(fig_3, tTextGraphName_Fig3_3D_Split, flg_graphSave, 0, size_Fig3);
    func_graphSave2(fig_4, tTextGraphName_Fig4_PCA_Split, flg_graphSave, 0, size_Fig4);
    func_graphSave2(fig_5, tTextGraphName_Fig5_3D_Comb, flg_graphSave, 0, size_Fig5);
    func_graphSave2(fig_6, tTextGraphName_Fig6_PCA_Comb, flg_graphSave, 0, size_Fig6);
    func_graphSave2(fig_7, tTextGraphName_Fig7_Stats, flg_graphSave, 0, size_Fig7);
end