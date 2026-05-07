close all; 
clear;     

% 機器設定ファイルの読み込み
if exist('make_status_bertec', 'file'), make_status_bertec; end
if exist('make_status_mocap_29', 'file'), make_status_mocap_29; end

%% =========================================================================
% 【プログラム概要: 半周期（Half-Cycle）の対称性評価】
% 1歩行周期を左右の半周期（Right Half: DS1+SS1, Left Half: DS2+SS2）に分割し、
% 歩行のキネマティックな非対称性を評価するプログラム。
% 左脚の関節角度配列を右脚の構造に合わせて反転（Flip）させることで、直接比較を可能にする。
% 前段で求めたGlobal SVD空間（シナジー空間）への射影を行い、左右が形成する
% 協調平面の類似度（法線ベクトルの内積）や、軌道形状の比率を統計的に解析する。
% =========================================================================

%% ===================================================================
% 1. 初期設定（入出力ファイルとパラメータ）
% ===================================================================
tTextFilePath  = './'; 
tTextFileName  = 'IBA_12'; 

% 読み込むファイル（SVD解析済みデータ）
tTextLoadDataName   = [tTextFileName, '_step4_SVD'];

% 最終的な保存ファイル名（全データ統合版）
tTextSaveMatName    = [tTextFileName, '_Final_AllData'];

% グラフ保存名設定
tTextGraphName_Check_Angles = [tTextFileName, '_HalfCycles_Fig1_Check_Angles']; 
tTextGraphName_Symmetry_3D  = [tTextFileName, '_HalfCycles_Fig2_Symmetry_3D'];
tTextGraphName_Symmetry_2D  = [tTextFileName, '_HalfCycles_Fig3_Symmetry_2D'];
tTextGraphName_Proj_DS      = [tTextFileName, '_HalfCycles_Fig4_Proj_DS'];
tTextGraphName_Proj_SS      = [tTextFileName, '_HalfCycles_Fig5_Proj_SS'];
tTextGraphName_Stats        = [tTextFileName, '_HalfCycles_Fig6_Stats_Combined']; 

flg_graphSave = 1; % グラフ保存フラグ

% --- 色の設定 ---
color_R_Stance         = [1, 0, 0];       % 赤 (右脚半周期)
color_L_Stance_Flipped = [0, 0, 0.8];     % 青 (左脚半周期・反転済み)
color_Full_Traj        = [0.8, 0.8, 0.8]; % 薄いグレー (全周期軌道)
color_L_Raw            = [0.7, 0.7, 0.7]; % 少し濃いグレー (左脚半周期・生データ)

%% ===================================================================
% 2. SVDデータの読込とバリデーション
% ===================================================================
fprintf('Loading SVD data from %s.mat...\n', tTextLoadDataName);
load(tTextLoadDataName); 

% 必要な変数がワークスペースに存在するか確認する
required_vars = {'averaged_DS1', 'averaged_SS1', 'averaged_DS2', 'averaged_SS2', 'W_proj', 'U_with_lambda', 'normal_DS1', 'normal_SS1'};
for i = 1:length(required_vars)
    if ~exist(required_vars{i}, 'var')
        error('必要な変数 "%s" が見つからない。SVD.m を実行してデータを更新すること。', required_vars{i});
    end
end

% データ列の定義
angleNames_SVD = {'L Foot', 'L Shank', 'L Thigh', 'Trunk', 'R Thigh', 'R Shank', 'R Foot'};

% ★重要: 左脚データを右脚データと直接比較するための反転インデックス
% 左足(1)→右足(7), 左下腿(2)→右下腿(6), 左大腿(3)→右大腿(5), 体幹(4)はそのまま
flip_map_indices = [7, 6, 5, 4, 3, 2, 1];

%% ===================================================================
% 3. 半周期（Half-Cycle）データの構築と反転処理
% ===================================================================
% 右脚基準の半周期（DS1 + SS1）と、左脚基準の半周期（DS2 + SS2）を構築する。
fprintf('Constructing Half-Cycles from Phase Data...\n');

% SVD空間に合わせて列順序を整理
newOrder = [7, 5, 3, 1, 2, 4, 6]; 
data_DS1 = averaged_DS1(:, newOrder);
data_SS1 = averaged_SS1(:, newOrder);
data_DS2 = averaged_DS2(:, newOrder);
data_SS2 = averaged_SS2(:, newOrder);

% 半周期データの結合
Right_Half_Data = [data_DS1; data_SS1];                 % 右脚の立脚期
Left_Half_Data_Raw = [data_DS2; data_SS2];              % 左脚の立脚期（生データ）
Left_Half_Data_Flipped = Left_Half_Data_Raw(:, flip_map_indices); % 左脚を右脚構造に反転

[T_half, N_dim] = size(Right_Half_Data);

%% ===================================================================
% 4. [Fig 1] 各体節角度の左右対称性の確認
% ===================================================================
% 反転処理が正しく行われ、左右の関節角度が重なるか（対称か）を時系列で確認する。
fig1 = figure(1); clf(fig1);
set(fig1, 'Name', 'Check_Angles_Symmetry', 'NumberTitle', 'off');
movegui(fig1, 'west');

plot_x_axis = 1:T_half;
max_val = max(max(abs(Right_Half_Data(:))), max(abs(Left_Half_Data_Flipped(:))));
ylim_common = [-max_val*1.1, max_val*1.1];

for j = 1:N_dim
    subplot(N_dim, 1, j); hold on;
    plot(plot_x_axis, Right_Half_Data(:, j), 'r-', 'LineWidth', 2);
    plot(plot_x_axis, Left_Half_Data_Flipped(:, j), 'b--', 'LineWidth', 2);
    xline(100, 'k:', 'LineWidth', 1); % DSとSSの境界線
    hold off; ylabel(angleNames_SVD{j}); ylim(ylim_common); xlim([1, T_half]); grid on;
    if j == 1
        title('Half-Cycle Symmetry Check (Angles)'); 
        legend('Right', 'Left(Flip)', 'Location', 'best'); 
    end
end
xlabel('Normalized Time (points)'); 

%% ===================================================================
% 5. Global SVD空間（シナジー空間）への射影と拘束平面の計算
% ===================================================================
% 前段で求めた全歩行周期の空間基底（W_proj）を用いて、半周期データを3次元空間に射影する。
Traj_Right_Proj        = Right_Half_Data * W_proj;
Traj_Left_Flipped_Proj = Left_Half_Data_Flipped * W_proj;
Traj_Full_Proj         = U_with_lambda(:, 1:3); 

% 左脚（反転済み）の各フェーズ（DS2, SS2）が形成する平面の法線ベクトルと原点距離を計算する。
% ※右脚側（DS1, SS1）の法線はSVD.mで計算済みの normal_DS1, normal_SS1 を使用。
calc_plane = @(X_data, W) calculate_plane_for_halfcycle(X_data, W);
[n_DS2_F, d_DS2_F] = calc_plane(Left_Half_Data_Flipped(1:100, :), W_proj);
[n_SS2_F, d_SS2_F] = calc_plane(Left_Half_Data_Flipped(101:200, :), W_proj);

%% ===================================================================
% 6. 平均平面の類似度計算 (コンソール出力)
% ===================================================================
% 左右の平面がどの程度平行かを示す指標として、法線ベクトル同士の内積を求める。
% 内積が1に近いほど、左右の運動制御戦略（シナジー）が対称であることを意味する。
dot_DS = abs(dot(normal_DS1, n_DS2_F));
dot_SS = abs(dot(normal_SS1, n_SS2_F));

angle_DS_deg = rad2deg(acos(min(dot_DS, 1)));
angle_SS_deg = rad2deg(acos(min(dot_SS, 1)));

fprintf('\n=== Plane Symmetry Analysis ===\n');
fprintf('DS Phase: Dot=%.4f, Angle=%.2f deg\n', dot_DS, angle_DS_deg);
fprintf('SS Phase: Dot=%.4f, Angle=%.2f deg\n', dot_SS, angle_SS_deg);
fprintf('===============================\n\n');

%% ===================================================================
% 7. [Fig 2] 3D軌道と共変動面の対称性プロット
% ===================================================================
% シナジー空間内における左右の軌道と平面を同一空間に描画し、視覚的に対称性を比較する。
fprintf('Plotting Fig 2: 3D Symmetry...\n');
fig_3D = figure(2); clf(fig_3D);
set(fig_3D, 'Name', 'Symmetry_3D', 'NumberTitle', 'off');
movegui(fig_3D, 'center');

hold on; grid on; pbaspect([1 1 1]); view(50, 10);

Traj_Left_Raw_Proj = Left_Half_Data_Raw * W_proj; 
idx_DS = 1:100; idx_SS = 101:200;

% 軌道の描画（実線:DS, 破線:SS）
h_R_DS = plot3(Traj_Right_Proj(idx_DS,1), Traj_Right_Proj(idx_DS,2), Traj_Right_Proj(idx_DS,3), '-', 'Color', color_R_Stance, 'LineWidth', 1.0);
h_R_SS = plot3(Traj_Right_Proj(idx_SS,1), Traj_Right_Proj(idx_SS,2), Traj_Right_Proj(idx_SS,3), '--', 'Color', color_R_Stance, 'LineWidth', 1.0);
h_L_DS = plot3(Traj_Left_Flipped_Proj(idx_DS,1), Traj_Left_Flipped_Proj(idx_DS,2), Traj_Left_Flipped_Proj(idx_DS,3), '-', 'Color', color_L_Stance_Flipped, 'LineWidth', 1.0);
h_L_SS = plot3(Traj_Left_Flipped_Proj(idx_SS,1), Traj_Left_Flipped_Proj(idx_SS,2), Traj_Left_Flipped_Proj(idx_SS,3), '--', 'Color', color_L_Stance_Flipped, 'LineWidth', 1.0);
h_Orig = plot3(Traj_Left_Raw_Proj(:,1), Traj_Left_Raw_Proj(:,2), Traj_Left_Raw_Proj(:,3), '-', 'Color', color_L_Raw, 'LineWidth', 1.0);

% スタート地点のマーカー
plot3(Traj_Right_Proj(1,1), Traj_Right_Proj(1,2), Traj_Right_Proj(1,3), 'o', 'MarkerFaceColor', color_R_Stance, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
plot3(Traj_Left_Flipped_Proj(1,1), Traj_Left_Flipped_Proj(1,2), Traj_Left_Flipped_Proj(1,3), 's', 'MarkerFaceColor', color_L_Stance_Flipped, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);

% 幾何学的最適化による平面（枠線）の描画
all_pts = [Traj_Right_Proj; Traj_Left_Flipped_Proj; Traj_Left_Raw_Proj];
min_pt = min(all_pts); max_pt = max(all_pts); center_pt = (min_pt + max_pt) / 2;
base_span = max(max_pt - min_pt) * 1.2; 
current_span_x = base_span; current_span_y = base_span; current_span_z = base_span;

planes_def = {
    struct('n', normal_DS1, 'd', d_DS1,   'col', color_R_Stance,         'style', '-'),  
    struct('n', normal_SS1, 'd', d_SS1,   'col', color_R_Stance,         'style', '--'), 
    struct('n', n_DS2_F,    'd', d_DS2_F, 'col', color_L_Stance_Flipped, 'style', '-'),  
    struct('n', n_SS2_F,    'd', d_SS2_F, 'col', color_L_Stance_Flipped, 'style', '--')  
};

final_intersects = cell(1, 4);
for iter = 1:50
    box_min = center_pt - [current_span_x, current_span_y, current_span_z]/2;
    box_max = center_pt + [current_span_x, current_span_y, current_span_z]/2;
    corners = [box_min(1), box_min(2), box_min(3); box_max(1), box_min(2), box_min(3); box_min(1), box_max(2), box_min(3); box_max(1), box_max(2), box_min(3); box_min(1), box_min(2), box_max(3); box_max(1), box_min(2), box_max(3); box_min(1), box_max(2), box_max(3); box_max(1), box_max(2), box_max(3)];
    edges = [1,2; 3,4; 5,6; 7,8; 1,3; 2,4; 5,7; 6,8; 1,5; 2,6; 3,7; 4,8];
    is_hex = false;
    for p = 1:4
        pts = [];
        for e = 1:12
            p1 = corners(edges(e,1),:); p2 = corners(edges(e,2),:);
            vec = p2-p1; denom = dot(planes_def{p}.n, vec);
            if abs(denom) > 1e-6
                t = (planes_def{p}.d - dot(planes_def{p}.n, p1)) / denom;
                if t >= -1e-5 && t <= 1+1e-5, pts = [pts; p1 + t*vec]; end
            end
        end
        if ~isempty(pts)
            pts = unique(round(pts,5), 'rows');
            if size(pts,1)>2
                c = mean(pts,1); [v,~] = pca(pts-c);
                coords = (pts-c)*v(:,1:2); [~,si] = sort(atan2(coords(:,2),coords(:,1)));
                pts = pts(si,:); pts(end+1,:) = pts(1,:);
            end
            if size(pts,1)-1 > 4, is_hex = true; end
        end
        final_intersects{p} = pts;
    end
    if ~is_hex, break; else, current_span_y = current_span_y * 1.05; end
end

plot3(final_intersects{1}(:,1), final_intersects{1}(:,2), final_intersects{1}(:,3), '-',  'Color', color_R_Stance, 'LineWidth', 0.5);
plot3(final_intersects{2}(:,1), final_intersects{2}(:,2), final_intersects{2}(:,3), '--', 'Color', color_R_Stance, 'LineWidth', 0.5);
plot3(final_intersects{3}(:,1), final_intersects{3}(:,2), final_intersects{3}(:,3), '-',  'Color', color_L_Stance_Flipped, 'LineWidth', 0.5);
plot3(final_intersects{4}(:,1), final_intersects{4}(:,2), final_intersects{4}(:,3), '--', 'Color', color_L_Stance_Flipped, 'LineWidth', 0.5);

% 軸の調整
lim_x = [center_pt(1)-current_span_x/2, center_pt(1)+current_span_x/2];
lim_y = [center_pt(2)-current_span_y/2, center_pt(2)+current_span_y/2];
lim_z = [center_pt(3)-current_span_z/2, center_pt(3)+current_span_z/2];
xlim(lim_x); ylim(lim_y); zlim(lim_z);
tick_step = 20;
xticks(ceil(lim_x(1)/tick_step)*tick_step : tick_step : floor(lim_x(2)/tick_step)*tick_step);
yticks(ceil(lim_y(1)/tick_step)*tick_step : tick_step : floor(lim_y(2)/tick_step)*tick_step);
zticks(ceil(lim_z(1)/tick_step)*tick_step : tick_step : floor(lim_z(2)/tick_step)*tick_step);

xlabel('Synergy 1'); ylabel('Synergy 2'); zlabel('Synergy 3');
legend([h_R_DS, h_R_SS, h_L_DS, h_L_SS, h_Orig], {'Right DS', 'Right SS', 'Left DS(Flip)', 'Left SS(Flip)', 'Left Raw'}, 'Location', 'bestoutside', 'FontSize', 8);

%% ===================================================================
% 8. [Fig 3] 2D射影による対称性プロット (各主成分ペア)
% ===================================================================
fprintf('Plotting Fig 3: 2D Symmetry...\n');
fig_2D = figure(3); clf(fig_2D);
set(fig_2D, 'Name', 'Symmetry_2D', 'NumberTitle', 'off');
movegui(fig_2D, 'center');

all_data_for_limit = [Traj_Full_Proj; Traj_Left_Raw_Proj];
max_val = max(abs(all_data_for_limit(:))) * 1.1;
lim_common = [-max_val, max_val];
ticks_val = ceil(lim_common(1)/tick_step)*tick_step : tick_step : floor(lim_common(2)/tick_step)*tick_step;

m1_r = Traj_Right_Proj(:,1); m2_r = Traj_Right_Proj(:,2); m3_r = Traj_Right_Proj(:,3);
m1_l = Traj_Left_Flipped_Proj(:,1); m2_l = Traj_Left_Flipped_Proj(:,2); m3_l = Traj_Left_Flipped_Proj(:,3);
m1_l_raw = Traj_Left_Raw_Proj(:,1); m2_l_raw = Traj_Left_Raw_Proj(:,2); m3_l_raw = Traj_Left_Raw_Proj(:,3);

subplot(1, 3, 1); hold on; grid on;
plot(Traj_Full_Proj(:,1), Traj_Full_Proj(:,2), '.', 'Color', [0.9 0.9 0.9], 'MarkerSize', 4);
plot(m1_l_raw, m2_l_raw, '-', 'Color', color_L_Raw, 'LineWidth', 1.0);
plot(m1_r, m2_r, '-', 'Color', color_R_Stance, 'LineWidth', 1.0);
plot(m1_l, m2_l, '-', 'Color', color_L_Stance_Flipped, 'LineWidth', 1.0);
axis equal; pbaspect([1 1 1]); xlim(lim_common); ylim(lim_common); xticks(ticks_val); yticks(ticks_val);
xlabel('Synergy 1'); ylabel('Synergy 2'); title('Synergy 1 vs 2');

subplot(1, 3, 2); hold on; grid on;
plot(Traj_Full_Proj(:,1), Traj_Full_Proj(:,3), '.', 'Color', [0.9 0.9 0.9], 'MarkerSize', 4);
plot(m1_l_raw, m3_l_raw, '-', 'Color', color_L_Raw, 'LineWidth', 1.0);
plot(m1_r, m3_r, '-', 'Color', color_R_Stance, 'LineWidth', 1.0);
plot(m1_l, m3_l, '-', 'Color', color_L_Stance_Flipped, 'LineWidth', 1.0);
axis equal; pbaspect([1 1 1]); xlim(lim_common); ylim(lim_common); xticks(ticks_val); yticks(ticks_val);
xlabel('Synergy 1'); ylabel('Synergy 3'); title('Synergy 1 vs 3');

subplot(1, 3, 3); hold on; grid on;
plot(Traj_Full_Proj(:,2), Traj_Full_Proj(:,3), '.', 'Color', [0.9 0.9 0.9], 'MarkerSize', 4);
plot(m2_l_raw, m3_l_raw, '-', 'Color', color_L_Raw, 'LineWidth', 1.0);
plot(m2_r, m3_r, '-', 'Color', color_R_Stance, 'LineWidth', 1.0);
plot(m2_l, m3_l, '-', 'Color', color_L_Stance_Flipped, 'LineWidth', 1.0);
axis equal; pbaspect([1 1 1]); xlim(lim_common); ylim(lim_common); xticks(ticks_val); yticks(ticks_val);
xlabel('Synergy 2'); ylabel('Synergy 3'); title('Synergy 2 vs 3');

%% ===================================================================
% 9. [Fig 4 & 5] 局所平面（Local Plane）への射影プロット
% ===================================================================
% 右脚側のDS平面・SS平面を基準とし、そこに左右の軌道を射影することで
% 「平面内でのループ形状」がどう違うか（対称性）を比較する。
fprintf('Plotting Fig 4 & 5: Projection onto Local Phase Planes...\n');

% 基準となる平面（右脚）の基底ベクトルを3D空間で定義する
v1_DS_3D = (V_DS1(:,1)' * W_proj)'; v1_DS_3D = v1_DS_3D / norm(v1_DS_3D);
v2_DS_3D = (V_DS1(:,2)' * W_proj)'; v2_DS_3D = v2_DS_3D / norm(v2_DS_3D);
v1_SS_3D = (V_SS1(:,1)' * W_proj)'; v1_SS_3D = v1_SS_3D / norm(v1_SS_3D);
v2_SS_3D = (V_SS1(:,2)' * W_proj)'; v2_SS_3D = v2_SS_3D / norm(v2_SS_3D);

% データの準備
traj_R_DS_3D = Traj_Right_Proj(1:100, :);
traj_R_SS_3D = Traj_Right_Proj(101:200, :);
traj_L_DS_3D = Traj_Left_Flipped_Proj(1:100, :);
traj_L_SS_3D = Traj_Left_Flipped_Proj(101:200, :);

% 基準平面への射影計算（内積を用いて平面座標系の(X, Y)に変換）
diff_R_DS = traj_R_DS_3D - p_DS1_3D; diff_L_DS = traj_L_DS_3D - p_DS1_3D;
proj_R_DS(:,1) = diff_R_DS * v1_DS_3D; proj_R_DS(:,2) = diff_R_DS * v2_DS_3D;
proj_L_DS(:,1) = diff_L_DS * v1_DS_3D; proj_L_DS(:,2) = diff_L_DS * v2_DS_3D;

diff_R_SS = traj_R_SS_3D - p_SS1_3D; diff_L_SS = traj_L_SS_3D - p_SS1_3D;
proj_R_SS(:,1) = diff_R_SS * v1_SS_3D; proj_R_SS(:,2) = diff_R_SS * v2_SS_3D;
proj_L_SS(:,1) = diff_L_SS * v1_SS_3D; proj_L_SS(:,2) = diff_L_SS * v2_SS_3D;

% プロット共通設定
all_proj_vals = [proj_R_DS; proj_L_DS; proj_R_SS; proj_L_SS];
lim_max = max(abs(all_proj_vals(:))) * 1.1;
lim_plane = [-lim_max, lim_max];
ticks_plane = ceil(lim_plane(1)/20)*20 : 20 : floor(lim_plane(2)/20)*20;

% --- [Fig 4] DS Phase Projection ---
fig4 = figure(4); clf(fig4);
set(fig4, 'Name', 'Proj_DS', 'NumberTitle', 'off');
movegui(fig4, 'center');

hold on; grid on; axis equal;
plot(0, 0, 'k+', 'MarkerSize', 10);
plot(proj_R_DS(:,1), proj_R_DS(:,2), '-', 'Color', color_R_Stance, 'LineWidth', 1.5);
plot(proj_L_DS(:,1), proj_L_DS(:,2), '-', 'Color', color_L_Stance_Flipped, 'LineWidth', 1.5);
plot(proj_R_DS(1,1), proj_R_DS(1,2), 'o', 'MarkerFaceColor', color_R_Stance, 'MarkerEdgeColor', 'k');
plot(proj_L_DS(1,1), proj_L_DS(1,2), 's', 'MarkerFaceColor', color_L_Stance_Flipped, 'MarkerEdgeColor', 'k');
xlim(lim_plane); ylim(lim_plane); xticks(ticks_plane); yticks(ticks_plane);
xlabel('Local Synergy 1'); ylabel('Local Synergy 2');

% --- [Fig 5] SS Phase Projection ---
fig5 = figure(5); clf(fig5);
set(fig5, 'Name', 'Proj_SS', 'NumberTitle', 'off');
movegui(fig5, 'center');

hold on; grid on; axis equal;
plot(0, 0, 'k+', 'MarkerSize', 10);
plot(proj_R_SS(:,1), proj_R_SS(:,2), '--', 'Color', color_R_Stance, 'LineWidth', 1.5);
plot(proj_L_SS(:,1), proj_L_SS(:,2), '--', 'Color', color_L_Stance_Flipped, 'LineWidth', 1.5);
plot(proj_R_SS(1,1), proj_R_SS(1,2), 'o', 'MarkerFaceColor', color_R_Stance, 'MarkerEdgeColor', 'k');
plot(proj_L_SS(1,1), proj_L_SS(1,2), 's', 'MarkerFaceColor', color_L_Stance_Flipped, 'MarkerEdgeColor', 'k');
xlim(lim_plane); ylim(lim_plane); xticks(ticks_plane); yticks(ticks_plane);
xlabel('Local Synergy 1'); ylabel('Local Synergy 2');

%% ===================================================================
% 10. 全ストライドに対する統計解析 (DS/SSフェーズごと)
% ===================================================================
% 個々のストライドごとに平面の内積（類似度）と軌道比率（Ratio）を計算する。
fprintf('Performing Stride-by-Stride Symmetry Analysis (Corrected SVD Space)...\n');

if ~exist('all_norm_DS1', 'var') || ~exist('all_norm_DS2', 'var')
    warning('ストライドデータが見つからない。');
    n_strides_valid = 0;
else
    n_strides_valid = size(all_norm_DS1, 3);
end

res_Sym_Dot_DS = zeros(n_strides_valid, 1); res_Sym_Ang_DS = zeros(n_strides_valid, 1);
res_Sym_Dot_SS = zeros(n_strides_valid, 1); res_Sym_Ang_SS = zeros(n_strides_valid, 1);
res_Ratio_R_DS = zeros(n_strides_valid, 1); res_Ratio_L_DS = zeros(n_strides_valid, 1);
res_Ratio_R_SS = zeros(n_strides_valid, 1); res_Ratio_L_SS = zeros(n_strides_valid, 1);

calc_plane_local = @(X) calculate_plane_for_halfcycle(X, W_proj);

for k = 1:n_strides_valid
    % --- 1. DS Phase Analysis ---
    dat_R_DS = all_norm_DS1(:, newOrder, k);
    dat_L_DS = all_norm_DS2(:, newOrder, k); 
    dat_L_DS_F = dat_L_DS(:, flip_map_indices); % 左脚の反転
    
    [n_R_DS_3D, ~] = calc_plane_local(dat_R_DS);
    [n_L_DS_3D, ~] = calc_plane_local(dat_L_DS_F);
    
    dot_ds = abs(dot(n_R_DS_3D, n_L_DS_3D)); 
    res_Sym_Dot_DS(k) = dot_ds; 
    res_Sym_Ang_DS(k) = rad2deg(acos(min(dot_ds, 1)));

    % Trajectory Ratio の計算 (Global空間射影データの分散比)
    traj_R_DS_proj = (dat_R_DS - mean(dat_R_DS)) * W_proj;
    traj_L_DS_proj = (dat_L_DS_F - mean(dat_L_DS_F)) * W_proj;
    
    [~, ~, lat_R_DS] = pca(traj_R_DS_proj); 
    if length(lat_R_DS) >= 2, res_Ratio_R_DS(k) = lat_R_DS(1)/lat_R_DS(2); else, res_Ratio_R_DS(k) = NaN; end
    
    [~, ~, lat_L_DS] = pca(traj_L_DS_proj);
    if length(lat_L_DS) >= 2, res_Ratio_L_DS(k) = lat_L_DS(1)/lat_L_DS(2); else, res_Ratio_L_DS(k) = NaN; end

    % --- 2. SS Phase Analysis ---
    dat_R_SS = all_norm_SS1(:, newOrder, k);
    dat_L_SS = all_norm_SS2(:, newOrder, k); 
    dat_L_SS_F = dat_L_SS(:, flip_map_indices); 
    
    [n_R_SS_3D, ~] = calc_plane_local(dat_R_SS);
    [n_L_SS_3D, ~] = calc_plane_local(dat_L_SS_F);
    
    dot_ss = abs(dot(n_R_SS_3D, n_L_SS_3D)); 
    res_Sym_Dot_SS(k) = dot_ss; 
    res_Sym_Ang_SS(k) = rad2deg(acos(min(dot_ss, 1)));

    traj_R_SS_proj = (dat_R_SS - mean(dat_R_SS)) * W_proj;
    traj_L_SS_proj = (dat_L_SS_F - mean(dat_L_SS_F)) * W_proj;
    
    [~, ~, lat_R_SS] = pca(traj_R_SS_proj);
    if length(lat_R_SS) >= 2, res_Ratio_R_SS(k) = lat_R_SS(1)/lat_R_SS(2); else, res_Ratio_R_SS(k) = NaN; end
    
    [~, ~, lat_L_SS] = pca(traj_L_SS_proj);
    if length(lat_L_SS) >= 2, res_Ratio_L_SS(k) = lat_L_SS(1)/lat_L_SS(2); else, res_Ratio_L_SS(k) = NaN; end
end

%% ===================================================================
% 11. [Fig 6] ストライドごとの統計結果プロット (Combined 2x2)
% ===================================================================
fprintf('Plotting Fig 6: Combined Statistics (Corrected)...\n');
fig6 = figure(6); clf(fig6);
set(fig6, 'Name', 'Stats_Combined', 'NumberTitle', 'off');
movegui(fig6, 'center');

col_BoxLine = [0.3 0.3 0.3];
col_Sim     = [0.6 0.6 0.6];
col_R       = color_R_Stance;
col_L       = color_L_Stance_Flipped;

% 散布図のプロットが重ならないようにX座標を分散（Jitter処理）させる
rng(1); 
jitter_val = 0.15;
x_1 = 1 + (rand(n_strides_valid, 1) - 0.5) * jitter_val;
x_2 = 2 + (rand(n_strides_valid, 1) - 0.5) * jitter_val;

calc_ylim = @(x) [min(x) - 0.1*range(x), max(x) + 0.1*range(x)];

% 1. Dot Product (平面の類似度)
subplot(2, 2, 1); hold on; grid on;
data_Dot = [res_Sym_Dot_DS; res_Sym_Dot_SS];
boxplot([res_Sym_Dot_DS, res_Sym_Dot_SS], 'Positions', [1, 2], 'Labels', {'DS Phase', 'SS Phase'}, 'Symbol', '', 'Colors', col_BoxLine);
scatter(x_1, res_Sym_Dot_DS, 20, col_Sim, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(x_2, res_Sym_Dot_SS, 20, col_Sim, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('Dot Product'); title('Plane Similarity (Dot Product)'); xlim([0.5, 2.5]);
if max(data_Dot)>min(data_Dot), ylim(calc_ylim(data_Dot)); else, ylim([0.9, 1.05]); end

% 2. Angle (平面のなす角)
subplot(2, 2, 2); hold on; grid on;
data_Ang = [res_Sym_Ang_DS; res_Sym_Ang_SS];
boxplot([res_Sym_Ang_DS, res_Sym_Ang_SS], 'Positions', [1, 2], 'Labels', {'DS Phase', 'SS Phase'}, 'Symbol', '', 'Colors', col_BoxLine);
scatter(x_1, res_Sym_Ang_DS, 20, col_Sim, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(x_2, res_Sym_Ang_SS, 20, col_Sim, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('Angle (deg)'); title('Plane Symmetry (Angle)'); xlim([0.5, 2.5]);
if max(data_Ang)>min(data_Ang), r=calc_ylim(data_Ang); ylim([max(0, r(1)), r(2)]); else, ylim([0, 5]); end

% 3. Ratio (DS Phase: 軌道形状の比率)
subplot(2, 2, 3); hold on; grid on;
data_RDS = [res_Ratio_R_DS; res_Ratio_L_DS];
boxplot([res_Ratio_R_DS, res_Ratio_L_DS], 'Positions', [1, 2], 'Labels', {'Right', 'Left(Flip)'}, 'Symbol', '', 'Colors', col_BoxLine);
scatter(x_1, res_Ratio_R_DS, 20, col_R, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(x_2, res_Ratio_L_DS, 20, col_L, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('Ratio (\lambda_1/\lambda_2)'); title('DS Phase: Trajectory Ratio'); xlim([0.5, 2.5]);
if max(data_RDS)>min(data_RDS), ylim(calc_ylim(data_RDS)); end

% 4. Ratio (SS Phase: 軌道形状の比率)
subplot(2, 2, 4); hold on; grid on;
data_RSS = [res_Ratio_R_SS; res_Ratio_L_SS];
boxplot([res_Ratio_R_SS, res_Ratio_L_SS], 'Positions', [1, 2], 'Labels', {'Right', 'Left(Flip)'}, 'Symbol', '', 'Colors', col_BoxLine);
scatter(x_1, res_Ratio_R_SS, 20, col_R, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(x_2, res_Ratio_L_SS, 20, col_L, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('Ratio (\lambda_1/\lambda_2)'); title('SS Phase: Trajectory Ratio'); xlim([0.5, 2.5]);
if max(data_RSS)>min(data_RSS), ylim(calc_ylim(data_RSS)); end

sgtitle('Stride-by-Stride Half-Cycle Symmetry & Ratio', 'FontSize', 14, 'FontWeight', 'bold');

%% ===================================================================
% 12. 統計表データのコンソール出力
% ===================================================================
fprintf('\n===============================================================\n');
fprintf('  Statistical Table: Stride-by-Stride Analysis (Mean ± SD)\n');
fprintf('===============================================================\n');
fprintf('%-25s | %s\n', 'Metric', 'Value (Mean +/- SD)');
fprintf('---------------------------------------------------------------\n');

calc_stats = @(data) sprintf('%.2f ± %.2f', mean(data, 'omitnan'), std(data, 'omitnan'));
calc_stats_long = @(data) sprintf('%.4f ± %.4f', mean(data, 'omitnan'), std(data, 'omitnan'));

fprintf('lambda1^2/lambda2^2 (DS1) | %s\n', calc_stats(res_Ratio_R_DS)); % Right DS
fprintf('lambda1^2/lambda2^2 (SS1) | %s\n', calc_stats(res_Ratio_R_SS)); % Right SS
fprintf('lambda1^2/lambda2^2 (DS2) | %s\n', calc_stats(res_Ratio_L_DS)); % Left DS
fprintf('lambda1^2/lambda2^2 (SS2) | %s\n', calc_stats(res_Ratio_L_SS)); % Left SS
fprintf('---------------------------------------------------------------\n');
fprintf('Dot product (DS)          | %s\n', calc_stats_long(res_Sym_Dot_DS));
fprintf('Dot product (SS)          | %s\n', calc_stats_long(res_Sym_Dot_SS));
fprintf('===============================================================\n\n');

%% ===================================================================
% 13. 計算結果のデータ保存
% ===================================================================
fprintf('Saving comprehensive dataset to %s.mat...\n', tTextSaveMatName);

% 既存の変数をすべて保存しつつ、今回計算したストライド別解析結果を追加する
save(tTextSaveMatName, ...
    'data', 'theta', 'centeredtheta', 'meanposture', ...
    'averaged_stride', 'std_stride', 'all_normalized_strides', ...
    'averaged_DS1', 'averaged_SS1', 'averaged_DS2', 'averaged_SS2', ...
    'all_norm_DS1', 'all_norm_SS1', 'all_norm_DS2', 'all_norm_SS2', ...
    'mean_orbit_R', 'mean_orbit_L', 'ratio_R', 'ratio_L', 'area_R', 'area_L', ...
    'normal_dot_product', 'score_R', 'score_L', 'coeff_R', 'coeff_L', ...
    'explained_R', 'explained_L', 'res_Ratio_R', 'res_Area_R', 'res_Ratio_L', ...
    'res_Area_L', 'res_Planarity_R', 'res_Planarity_L', 'res_DotProd', 'res_Angle', ...
    'stride_svd', 'U_with_lambda', 'V_spatial', 'S_global', 'U_global', 'V_global', 'W_proj', ...
    'FLIP_MODES_GLOBAL', 'FLIP_MODES_DS1', 'FLIP_MODES_SS1', 'FLIP_MODES_DS2', 'FLIP_MODES_SS2', ...
    'V_DS1', 'S_DS1', 'U_DS1', 'normal_DS1', 'd_DS1', 'p_DS1_3D', 'P_DS1_Proj', ...
    'V_SS1', 'S_SS1', 'U_SS1', 'normal_SS1', 'd_SS1', 'p_SS1_3D', 'P_SS1_Proj', ...
    'V_DS2', 'S_DS2', 'U_DS2', 'normal_DS2', 'd_DS2', 'p_DS2_3D', 'P_DS2_Proj', ...
    'V_SS2', 'S_SS2', 'U_SS2', 'normal_SS2', 'd_SS2', 'p_SS2_3D', 'P_SS2_Proj', ...
    'Right_Half_Data', 'Left_Half_Data_Flipped', 'Traj_Right_Proj', 'Traj_Left_Flipped_Proj', ...
    'n_DS2_F', 'd_DS2_F', 'n_SS2_F', 'd_SS2_F', ...
    'res_Sym_Dot_DS', 'res_Sym_Ang_DS', 'res_Sym_Dot_SS', 'res_Sym_Ang_SS', ...
    'res_Ratio_R_DS', 'res_Ratio_L_DS', 'res_Ratio_R_SS', 'res_Ratio_L_SS');

fprintf('Processing Complete. Comprehensive results saved.\n');

%% ===================================================================
% 14. グラフ画像の一括保存
% ===================================================================
if flg_graphSave
    fprintf('Saving figures using func_graphSave2 (Batch Mode)...\n');
    size_Fig1   = [800, 1000];
    size_Fig3D  = [700, 300];
    size_Fig2D  = [600, 200];
    size_Proj   = [250, 225];  
    size_Stats  = [1000, 800]; 

    func_graphSave2(fig1, tTextGraphName_Check_Angles, flg_graphSave, 0, size_Fig1);
    func_graphSave2(fig_3D, tTextGraphName_Symmetry_3D, flg_graphSave, 0, size_Fig3D);
    func_graphSave2(fig_2D, tTextGraphName_Symmetry_2D, flg_graphSave, 0, size_Fig2D);
    
    func_graphSave2(fig4, tTextGraphName_Proj_DS, flg_graphSave, 0, size_Proj);
    func_graphSave2(fig5, tTextGraphName_Proj_SS, flg_graphSave, 0, size_Proj);

    func_graphSave2(fig6, tTextGraphName_Stats, flg_graphSave, 0, size_Stats);
end

%% ===================================================================
% 15. 内部関数定義
% ===================================================================
function [normal, d] = calculate_plane_for_halfcycle(X_data, W_global)
    % 与えられた半周期データ(X_data)から、Global空間における平面の方程式(normal, d)を計算する。
    X_mean = mean(X_data, 1);
    X_var = X_data - X_mean;
    
    [~, ~, V_sub] = svd(X_var, 'econ');
    
    % Globalな3D空間（シナジー空間）の基底ベクトルに変換
    v1_3D = (V_sub(:,1)' * W_global);
    v2_3D = (V_sub(:,2)' * W_global);
    
    % 外積から法線ベクトルを求める
    normal = cross(v1_3D, v2_3D); 
    normal = normal / norm(normal);
    
    % 平面の原点からの距離 d を求める
    p_3D = X_mean * W_global;
    d = dot(normal, p_3D);
end