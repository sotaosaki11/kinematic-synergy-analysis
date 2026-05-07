close all; 
clear;
clc;

%% =========================================================================
% 【プログラム概要: 片麻痺歩行 半周期対称性解析 (Trajectory Map, Step 5)】
% SVD解析（Step 4）の結果を読み込み、1歩行周期を左右の半周期
% (Right Half: DS1+SS1, Left Half: DS2+SS2) に分割して運動の非対称性を評価する。
% 患側（左脚）の関節角度配列を健側（右脚）の構造に合わせて反転（Flip）させることで、
% 同一のシナジー空間上での軌道形状や、両者が形成する協調平面の類似度（内積）を直接比較する。
% さらに、各フェーズの軌道を個別に分解し、視覚的な評価を行うための図を出力する。
% =========================================================================

%% ===================================================================
% 1. 初期設定 (入出力ファイルと解析パラメータ)
% ===================================================================
SubjectName    = 'KM_ID3'; 
tTextFilePath  = './'; 

% --- 入出力ファイル名 (直前のStep 4から引き継ぐ) ---
tTextLoadDataName   = [SubjectName, '_step4_SVD']; 
tTextSaveMatName    = [SubjectName, '_step5_Trajectory'];

% --- グラフ保存名設定 ---
tTextGraphName_Check_Angles = [SubjectName, '_HalfCycle_Fig1_Check_Angles']; 
tTextGraphName_Symmetry_3D  = [SubjectName, '_HalfCycle_Fig2_Symmetry_3D'];
tTextGraphName_Symmetry_2D  = [SubjectName, '_HalfCycle_Fig3_Symmetry_2D'];
tTextGraphName_Proj_DS      = [SubjectName, '_HalfCycle_Fig4_Proj_DS'];
tTextGraphName_Proj_SS      = [SubjectName, '_HalfCycle_Fig5_Proj_SS'];
tTextGraphName_Sym_Stats    = [SubjectName, '_HalfCycle_Fig6_Symmetry_Stats']; 

% 軌道と平面の分解図（可視化用）
tTextGraphName_Decomp_DS1   = [SubjectName, '_HalfCycle_Fig7_DS1_Only'];
tTextGraphName_Decomp_SS1   = [SubjectName, '_HalfCycle_Fig8_SS1_Only'];
tTextGraphName_Decomp_DS2   = [SubjectName, '_HalfCycle_Fig9_DS2_Only'];
tTextGraphName_Decomp_SS2   = [SubjectName, '_HalfCycle_Fig10_SS2_Only'];
tTextGraphName_Planes_Only  = [SubjectName, '_HalfCycle_Fig11_Planes_Only'];

flg_graphSave = 1; % グラフ保存フラグ

% --- 色の設定 (健常者用と統一) ---
color_R_Stance           = [1 0 0];               % 赤 (右脚/健側)
color_L_Stance_Flipped   = [0 0 0.8];             % 青 (左脚/患側・反転済み)
color_Full_Traj          = [0.8 0.8 0.8];         % 薄いグレー (背景全軌道)
color_L_Raw              = [0.7 0.7 0.7];         % 少し濃いグレー (左脚/患側・生データ)

%% ===================================================================
% 2. SVDデータの読込とバリデーション
% ===================================================================
fprintf('Loading data from %s.mat...\n', tTextLoadDataName);
if ~exist([tTextLoadDataName '.mat'], 'file')
    error('ファイル %s が見つからない。Step 4を実行すること。', tTextLoadDataName);
end

load(tTextLoadDataName);

% 必要な変数の存在チェック
if ~exist('valid_strides_L', 'var') || ~exist('averaged_DS1', 'var')
    error('必要な変数が不足している。Step 4までの処理を確認すること。');
end

% --- SVD空間の射影行列と全軌道データ ---
W_proj = V_spatial(:, 1:3); 
Traj_Full_Proj = U_with_lambda(:, 1:3);

% 患側(左脚)データを健側(右脚)データと直接比較するための反転インデックス
% 左足(1)→右足(7), 左下腿(2)→右下腿(6), 左大腿(3)→右大腿(5), 体幹(4)はそのまま
flip_map_indices = [7, 6, 5, 4, 3, 2, 1];

%% ===================================================================
% 3. 半周期（Half-Cycle）データの構築と反転処理
% ===================================================================
fprintf('Constructing Average Half-Cycles...\n');

% SVD空間の並びに合わせて各フェーズのデータを整理する
data_DS1 = averaged_DS1(:, newOrder);
data_SS1 = averaged_SS1(:, newOrder);
data_DS2 = averaged_DS2(:, newOrder);
data_SS2 = averaged_SS2(:, newOrder);

% 右脚（健側）の半周期データ (DS1 + SS1)
Right_Half_Data = [data_DS1; data_SS1];

% 左脚（患側）の半周期データ (DS2 + SS2) を構築し、関節順序を右脚基準に反転させる
Left_Half_Data_Raw = [data_DS2; data_SS2];
Left_Half_Data_Flipped = Left_Half_Data_Raw(:, flip_map_indices);

[T_half, N_dim] = size(Right_Half_Data);

%% ===================================================================
% 4. [Fig 1] 各体節角度の左右対称性の確認
% ===================================================================
% 反転処理が正しく行われ、左右の関節角度軌道が比較可能な状態になっているかを時系列で確認する。
fprintf('Plotting Fig 1: Angles...\n');
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
    xline(100, 'k:', 'LineWidth', 1);
    hold off; ylabel(angleNames{j}); ylim(ylim_common); xlim([1, T_half]); grid on;
    if j == 1
        title('Half-Cycle Symmetry Check (Angles)'); 
        legend('Right (Non-paretic)', 'Left (Paretic, Flipped)', 'Location', 'best'); 
    end
end
xlabel('Normalized Time (points)'); 

%% ===================================================================
% 5. Global SVD空間（シナジー空間）への射影と拘束平面の計算
% ===================================================================
% 全歩行周期の空間基底（W_proj）を用いて、半周期データを3次元のシナジー空間に射影する。
Traj_Right_Proj = Right_Half_Data * W_proj;
Traj_Left_Flipped_Proj = Left_Half_Data_Flipped * W_proj;
Traj_Left_Raw_Proj = Left_Half_Data_Raw * W_proj; 

% 患側（反転済み）の各フェーズ（DS2, SS2）が形成する平面の法線ベクトルと原点距離を計算する。
% ※健側（DS1, SS1）の法線は Step 4 で計算済みの normal_DS1, normal_SS1 を使用する。
calc_plane = @(X_data, W) calculate_plane_for_halfcycle(X_data, W);
[n_DS2_F, d_DS2_F] = calc_plane(Left_Half_Data_Flipped(1:100, :), W_proj);
[n_SS2_F, d_SS2_F] = calc_plane(Left_Half_Data_Flipped(101:200, :), W_proj);

% 両脚が形成する平面がどれだけ平行に近いか（制御方略の対称性）を内積で評価する。
dot_DS = abs(dot(normal_DS1, n_DS2_F));
dot_SS = abs(dot(normal_SS1, n_SS2_F));
angle_DS_deg = rad2deg(acos(min(dot_DS, 1)));
angle_SS_deg = rad2deg(acos(min(dot_SS, 1)));

fprintf('Mean Plane Symmetry:\n  DS Dot=%.4f (%.1f deg)\n  SS Dot=%.4f (%.1f deg)\n', ...
        dot_DS, angle_DS_deg, dot_SS, angle_SS_deg);

%% ===================================================================
% 6. [Fig 2] 3D軌道と共変動面の対称性プロット
% ===================================================================
% シナジー空間内における左右の軌道と平面を同一空間に描画し、視覚的に対称性を比較する。
fprintf('Plotting Fig 2: 3D Symmetry (DS: Solid, SS: Dashed)...\n');
fig2 = figure(2); clf(fig2);
set(fig2, 'Name', 'Symmetry_3D', 'NumberTitle', 'off');
movegui(fig2, 'center');

hold on; grid on; pbaspect([1 1 1]); view(140, 10);

idx_DS = 1:100; idx_SS = 101:200;

% --- 1. 軌道の描画 (DS: 実線, SS: 破線) ---
h_R_DS_Traj = plot3(Traj_Right_Proj(idx_DS,1), Traj_Right_Proj(idx_DS,2), Traj_Right_Proj(idx_DS,3), '-', 'Color', color_R_Stance, 'LineWidth', 1.0);
h_R_SS_Traj = plot3(Traj_Right_Proj(idx_SS,1), Traj_Right_Proj(idx_SS,2), Traj_Right_Proj(idx_SS,3), '--', 'Color', color_R_Stance, 'LineWidth', 1.0);
h_L_DS_Traj = plot3(Traj_Left_Flipped_Proj(idx_DS,1), Traj_Left_Flipped_Proj(idx_DS,2), Traj_Left_Flipped_Proj(idx_DS,3), '-', 'Color', color_L_Stance_Flipped, 'LineWidth', 1.0);
h_L_SS_Traj = plot3(Traj_Left_Flipped_Proj(idx_SS,1), Traj_Left_Flipped_Proj(idx_SS,2), Traj_Left_Flipped_Proj(idx_SS,3), '--', 'Color', color_L_Stance_Flipped, 'LineWidth', 1.0);
h_Original  = plot3(Traj_Left_Raw_Proj(:,1), Traj_Left_Raw_Proj(:,2), Traj_Left_Raw_Proj(:,3), '-', 'Color', color_L_Raw, 'LineWidth', 1.0);

% スタート地点のマーカー
plot3(Traj_Right_Proj(1,1), Traj_Right_Proj(1,2), Traj_Right_Proj(1,3), 'o', 'MarkerFaceColor', color_R_Stance, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
plot3(Traj_Left_Flipped_Proj(1,1), Traj_Left_Flipped_Proj(1,2), Traj_Left_Flipped_Proj(1,3), 's', 'MarkerFaceColor', color_L_Stance_Flipped, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);

% --- 2. 描画範囲(ボックス)の最適化 ---
all_pts = [Traj_Right_Proj; Traj_Left_Flipped_Proj; Traj_Left_Raw_Proj];
min_pt = min(all_pts); max_pt = max(all_pts);
center_pt = (min_pt + max_pt) / 2;
max_span = max(max_pt - min_pt) * 1.3; 

box_min = center_pt - max_span/2; box_max = center_pt + max_span/2;
corners = [box_min(1), box_min(2), box_min(3); box_max(1), box_min(2), box_min(3); box_min(1), box_max(2), box_min(3); box_max(1), box_max(2), box_min(3); box_min(1), box_min(2), box_max(3); box_max(1), box_min(2), box_max(3); box_min(1), box_max(2), box_max(3); box_max(1), box_max(2), box_max(3)];
edges = [1,2; 3,4; 5,6; 7,8; 1,3; 2,4; 5,7; 6,8; 1,5; 2,6; 3,7; 4,8];

% --- 3. 平面枠線の計算と描画 ---
planes_def = {
    struct('n', normal_DS1, 'd', d_DS1,   'col', color_R_Stance,         'style', '-'),  
    struct('n', normal_SS1, 'd', d_SS1,   'col', color_R_Stance,         'style', '--'), 
    struct('n', n_DS2_F,    'd', d_DS2_F, 'col', color_L_Stance_Flipped, 'style', '-'),  
    struct('n', n_SS2_F,    'd', d_SS2_F, 'col', color_L_Stance_Flipped, 'style', '--')  
};

for p = 1:4
    n = planes_def{p}.n; d_val = planes_def{p}.d;
    pts = [];
    for i = 1:size(edges, 1)
        p1 = corners(edges(i,1), :); p2 = corners(edges(i,2), :);
        vec = p2 - p1; denom = dot(n, vec);
        if abs(denom) > 1e-6
            t = (d_val - dot(n, p1)) / denom;
            if t >= -1e-5 && t <= 1 + 1e-5, pts = [pts; p1 + t * vec]; end
        end
    end
    if ~isempty(pts)
        pts = unique(round(pts, 5), 'rows');
        if size(pts, 1) > 2
            c = mean(pts, 1); [coeff_p, ~] = pca(pts - c);
            coords_2d = (pts - c) * coeff_p(:, 1:2);
            angles = atan2(coords_2d(:,2), coords_2d(:,1));
            [~, sort_idx] = sort(angles);
            pts = pts(sort_idx, :); pts(end+1, :) = pts(1, :);
            plot3(pts(:,1), pts(:,2), pts(:,3), planes_def{p}.style, 'Color', planes_def{p}.col, 'LineWidth', 0.5);
        end
    end
end

% --- 4. 軸設定 ---
xlim([box_min(1), box_max(1)]); ylim([box_min(2), box_max(2)]); zlim([box_min(3), box_max(3)]);
tick_step = 20;
xticks(ceil(box_min(1)/tick_step)*tick_step : tick_step : floor(box_max(1)/tick_step)*tick_step);
yticks(ceil(box_min(2)/tick_step)*tick_step : tick_step : floor(box_max(2)/tick_step)*tick_step);
zticks(ceil(box_min(3)/tick_step)*tick_step : tick_step : floor(box_max(3)/tick_step)*tick_step);

xlabel('Synergy 1'); ylabel('Synergy 2'); zlabel('Synergy 3');
legend([h_R_DS_Traj, h_R_SS_Traj, h_L_DS_Traj, h_L_SS_Traj, h_Original], ...
       {'Right DS Trajectory/Plane', 'Right SS Trajectory/Plane', 'Left DS Trajectory/Plane', 'Left SS Trajectory/Plane', 'Left Trajectory Original'}, ...
       'Location', 'bestoutside', 'FontSize', 8);

%% ===================================================================
% 7. [Fig 3] 2D射影による対称性プロット (各主成分ペア)
% ===================================================================
fprintf('Plotting Fig 3: 2D Symmetry (Side-by-Side, Square)...\n');
fig3 = figure(3); clf(fig3);
set(fig3, 'Name', 'Symmetry_2D', 'NumberTitle', 'off');
movegui(fig3, 'center');

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
% 8. [Fig 4 & 5] 局所平面（Local Phase Plane）への射影プロット
% ===================================================================
% 健側のフェーズ別平面（DS/SS）を基準座標系とし、そこに両脚の軌道を射影することで
% 「平面内でのループ形状の歪み」を直接比較する。
fprintf('Plotting Fig 4 & 5: Projection onto Local Phase Planes...\n');

% 基準となる平面（健側）の基底ベクトルを定義
v1_DS_3D = (V_DS1(:,1)' * W_proj)'; v1_DS_3D = v1_DS_3D / norm(v1_DS_3D);
v2_DS_3D = (V_DS1(:,2)' * W_proj)'; v2_DS_3D = v2_DS_3D / norm(v2_DS_3D);
v1_SS_3D = (V_SS1(:,1)' * W_proj)'; v1_SS_3D = v1_SS_3D / norm(v1_SS_3D);
v2_SS_3D = (V_SS1(:,2)' * W_proj)'; v2_SS_3D = v2_SS_3D / norm(v2_SS_3D);

traj_R_DS_3D = Traj_Right_Proj(1:100, :); traj_R_SS_3D = Traj_Right_Proj(101:200, :);
traj_L_DS_3D = Traj_Left_Flipped_Proj(1:100, :); traj_L_SS_3D = Traj_Left_Flipped_Proj(101:200, :);

% 基準平面への射影（内積による2D座標への変換）
diff_R_DS = traj_R_DS_3D - p_DS1_3D; diff_L_DS = traj_L_DS_3D - p_DS1_3D;
proj_R_DS(:,1) = diff_R_DS * v1_DS_3D; proj_R_DS(:,2) = diff_R_DS * v2_DS_3D;
proj_L_DS(:,1) = diff_L_DS * v1_DS_3D; proj_L_DS(:,2) = diff_L_DS * v2_DS_3D;

diff_R_SS = traj_R_SS_3D - p_SS1_3D; diff_L_SS = traj_L_SS_3D - p_SS1_3D;
proj_R_SS(:,1) = diff_R_SS * v1_SS_3D; proj_R_SS(:,2) = diff_R_SS * v2_SS_3D;
proj_L_SS(:,1) = diff_L_SS * v1_SS_3D; proj_L_SS(:,2) = diff_L_SS * v2_SS_3D;

all_proj_vals = [proj_R_DS; proj_L_DS; proj_R_SS; proj_L_SS];
lim_max = max(abs(all_proj_vals(:))) * 1.1;
lim_plane = [-lim_max, lim_max];
ticks_plane = ceil(lim_plane(1)/20)*20 : 20 : floor(lim_plane(2)/20)*20;

% --- [Fig 4] DS Phase Projection ---
fig4 = figure(4); clf(fig4);
set(fig4, 'Name', 'Proj_DS', 'NumberTitle', 'off'); movegui(fig4, 'center');
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
set(fig5, 'Name', 'Proj_SS', 'NumberTitle', 'off'); movegui(fig5, 'center');
hold on; grid on; axis equal;
plot(0, 0, 'k+', 'MarkerSize', 10);
plot(proj_R_SS(:,1), proj_R_SS(:,2), '--', 'Color', color_R_Stance, 'LineWidth', 1.5);
plot(proj_L_SS(:,1), proj_L_SS(:,2), '--', 'Color', color_L_Stance_Flipped, 'LineWidth', 1.5);
plot(proj_R_SS(1,1), proj_R_SS(1,2), 'o', 'MarkerFaceColor', color_R_Stance, 'MarkerEdgeColor', 'k');
plot(proj_L_SS(1,1), proj_L_SS(1,2), 's', 'MarkerFaceColor', color_L_Stance_Flipped, 'MarkerEdgeColor', 'k');
xlim(lim_plane); ylim(lim_plane); xticks(ticks_plane); yticks(ticks_plane);
xlabel('Local Synergy 1'); ylabel('Local Synergy 2');

%% ===================================================================
% 9. 全ストライドに対する統計解析 (DS/SSフェーズごと)
% ===================================================================
% 個々のストライドごとに平面の内積（類似度）と軌道比率（Ratio）を計算する。
fprintf('Performing Stride-by-Stride Symmetry Analysis...\n');

n_strides_L = length(valid_strides_L);
N_PHASE = 100;
all_norm_L_DS1 = zeros(N_PHASE, size(all_norm_DS1,2), n_strides_L);
all_norm_L_SS1 = zeros(N_PHASE, size(all_norm_SS1,2), n_strides_L);

for i = 1:n_strides_L
    st = valid_strides_L(i);
    if ~isempty(st.raw_DS1), all_norm_L_DS1(:,:,i) = interp1(linspace(0,1,size(st.raw_DS1,1)), st.raw_DS1(:, newOrder), linspace(0,1,N_PHASE)', 'pchip'); end
    if ~isempty(st.raw_SS1), all_norm_L_SS1(:,:,i) = interp1(linspace(0,1,size(st.raw_SS1,1)), st.raw_SS1(:, newOrder), linspace(0,1,N_PHASE)', 'pchip'); end
end

n_strides_valid = min(size(all_norm_DS1, 3), n_strides_L);
res_Sym_Dot_DS = zeros(n_strides_valid, 1); res_Sym_Ang_DS = zeros(n_strides_valid, 1);
res_Sym_Dot_SS = zeros(n_strides_valid, 1); res_Sym_Ang_SS = zeros(n_strides_valid, 1);
res_Ratio_R_DS = zeros(n_strides_valid, 1); res_Ratio_L_DS = zeros(n_strides_valid, 1);
res_Ratio_R_SS = zeros(n_strides_valid, 1); res_Ratio_L_SS = zeros(n_strides_valid, 1);

calc_plane_local = @(X) calculate_plane_for_halfcycle(X, W_proj);

for k = 1:n_strides_valid
    % --- DS Phase Analysis ---
    dat_R_DS = all_norm_DS1(:, newOrder, k);
    dat_L_DS = all_norm_L_DS1(:, :, k); 
    dat_L_DS_F = dat_L_DS(:, flip_map_indices); 
    
    [n_R_DS_3D, ~] = calc_plane_local(dat_R_DS); [n_L_DS_3D, ~] = calc_plane_local(dat_L_DS_F);
    
    dot_ds = abs(dot(n_R_DS_3D, n_L_DS_3D)); 
    res_Sym_Dot_DS(k) = dot_ds; res_Sym_Ang_DS(k) = rad2deg(acos(min(dot_ds, 1)));

    traj_R_DS_proj = (dat_R_DS - mean(dat_R_DS)) * W_proj; traj_L_DS_proj = (dat_L_DS_F - mean(dat_L_DS_F)) * W_proj;
    
    [~, ~, lat_R_DS] = pca(traj_R_DS_proj); 
    if length(lat_R_DS) >= 2, res_Ratio_R_DS(k) = lat_R_DS(1)/lat_R_DS(2); else, res_Ratio_R_DS(k) = NaN; end
    [~, ~, lat_L_DS] = pca(traj_L_DS_proj);
    if length(lat_L_DS) >= 2, res_Ratio_L_DS(k) = lat_L_DS(1)/lat_L_DS(2); else, res_Ratio_L_DS(k) = NaN; end

    % --- SS Phase Analysis ---
    dat_R_SS = all_norm_SS1(:, newOrder, k);
    dat_L_SS = all_norm_L_SS1(:, :, k); 
    dat_L_SS_F = dat_L_SS(:, flip_map_indices); 
    
    [n_R_SS_3D, ~] = calc_plane_local(dat_R_SS); [n_L_SS_3D, ~] = calc_plane_local(dat_L_SS_F);
    
    dot_ss = abs(dot(n_R_SS_3D, n_L_SS_3D)); 
    res_Sym_Dot_SS(k) = dot_ss; res_Sym_Ang_SS(k) = rad2deg(acos(min(dot_ss, 1)));

    traj_R_SS_proj = (dat_R_SS - mean(dat_R_SS)) * W_proj; traj_L_SS_proj = (dat_L_SS_F - mean(dat_L_SS_F)) * W_proj;
    
    [~, ~, lat_R_SS] = pca(traj_R_SS_proj);
    if length(lat_R_SS) >= 2, res_Ratio_R_SS(k) = lat_R_SS(1)/lat_R_SS(2); else, res_Ratio_R_SS(k) = NaN; end
    [~, ~, lat_L_SS] = pca(traj_L_SS_proj);
    if length(lat_L_SS) >= 2, res_Ratio_L_SS(k) = lat_L_SS(1)/lat_L_SS(2); else, res_Ratio_L_SS(k) = NaN; end
end

%% ===================================================================
% 10. [Fig 6] 統計プロット
% ===================================================================
fprintf('Plotting Fig 6: Stride-by-Stride Symmetry Statistics...\n');
fig6 = figure(6); clf(fig6);
set(fig6, 'Name', 'Symmetry_Stats', 'NumberTitle', 'off'); movegui(fig6, 'center');

col_BoxLine = [0.2 0.2 0.2]; col_R = color_R_Stance; col_L = color_L_Stance_Flipped; col_Sim = [0.5, 0.5, 0.5];
rng(1); jit_width = 0.15;
x_1 = 1 + (rand(n_strides_valid,1)-0.5)*jit_width; x_2 = 2 + (rand(n_strides_valid,1)-0.5)*jit_width;
calc_ylim = @(d) [min(d) - 0.1*range(d), max(d) + 0.1*range(d)];

% 1. Dot Product
subplot(2, 2, 1); hold on; grid on;
data_Dot = [res_Sym_Dot_DS; res_Sym_Dot_SS];
boxplot([res_Sym_Dot_DS, res_Sym_Dot_SS], 'Positions', [1, 2], 'Labels', {'DS Phase', 'SS Phase'}, 'Symbol', '', 'Colors', col_BoxLine);
scatter(x_1, res_Sym_Dot_DS, 30, col_Sim, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(x_2, res_Sym_Dot_SS, 30, col_Sim, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('Dot Product'); title('Plane Similarity (Dot Product)'); xlim([0.5, 2.5]);
if max(data_Dot)>min(data_Dot), ylim(calc_ylim(data_Dot)); else, ylim([0.9, 1.05]); end

% 2. Angle
subplot(2, 2, 2); hold on; grid on;
data_Ang = [res_Sym_Ang_DS; res_Sym_Ang_SS];
boxplot([res_Sym_Ang_DS, res_Sym_Ang_SS], 'Positions', [1, 2], 'Labels', {'DS Phase', 'SS Phase'}, 'Symbol', '', 'Colors', col_BoxLine);
scatter(x_1, res_Sym_Ang_DS, 30, col_Sim, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(x_2, res_Sym_Ang_SS, 30, col_Sim, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('Angle (deg)'); title('Plane Symmetry (Angle)'); xlim([0.5, 2.5]);
if max(data_Ang)>min(data_Ang), r=calc_ylim(data_Ang); ylim([max(0, r(1)), r(2)]); else, ylim([0, 5]); end

% 3. Ratio (DS)
subplot(2, 2, 3); hold on; grid on;
data_RDS = [res_Ratio_R_DS; res_Ratio_L_DS];
boxplot([res_Ratio_R_DS, res_Ratio_L_DS], 'Positions', [1, 2], 'Labels', {'Right', 'Left(Flip)'}, 'Symbol', '', 'Colors', col_BoxLine);
scatter(x_1, res_Ratio_R_DS, 30, col_R, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(x_2, res_Ratio_L_DS, 30, col_L, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('Ratio (\lambda_1/\lambda_2)'); title('DS Phase: Trajectory Ratio'); xlim([0.5, 2.5]);
if max(data_RDS)>min(data_RDS), ylim(calc_ylim(data_RDS)); end

% 4. Ratio (SS)
subplot(2, 2, 4); hold on; grid on;
data_RSS = [res_Ratio_R_SS; res_Ratio_L_SS];
boxplot([res_Ratio_R_SS, res_Ratio_L_SS], 'Positions', [1, 2], 'Labels', {'Right', 'Left(Flip)'}, 'Symbol', '', 'Colors', col_BoxLine);
scatter(x_1, res_Ratio_R_SS, 30, col_R, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(x_2, res_Ratio_L_SS, 30, col_L, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('Ratio (\lambda_1/\lambda_2)'); title('SS Phase: Trajectory Ratio'); xlim([0.5, 2.5]);
if max(data_RSS)>min(data_RSS), ylim(calc_ylim(data_RSS)); end

sgtitle('Stride-by-Stride Half-Cycle Symmetry & Ratio', 'FontSize', 14, 'FontWeight', 'bold');

%% ===================================================================
% 11. [Fig 7-11] 視覚的分解プロット (Decomposed Trajectories & Planes)
% ===================================================================
% Fig 2のように全ての軌道と平面が重なると視認性が低下するため、
% 各フェーズ（DS/SS）ごとの軌道と平面を独立させた図を生成する。
fprintf('Plotting Fig 7-11: Decomposed Trajectories & Planes...\n');

for mode_idx = 1:5
    fig_num = 6 + mode_idx; 
    f_h = figure(fig_num); clf(f_h);
    
    switch mode_idx
        case 1, tName='DS1_Only'; target_p=1; 
        case 2, tName='SS1_Only'; target_p=2; 
        case 3, tName='DS2_Only'; target_p=3; 
        case 4, tName='SS2_Only'; target_p=4; 
        case 5, tName='Planes_Only'; target_p=0; 
    end
    set(f_h, 'Name', tName, 'NumberTitle', 'off');
    movegui(f_h, 'center');
    
    hold on; grid on; pbaspect([1 1 1]); view(140, 10);
    
    % --- 軌道の描画 ---
    if mode_idx == 1 % DS1 (Right)
        plot3(Traj_Right_Proj(idx_DS,1), Traj_Right_Proj(idx_DS,2), Traj_Right_Proj(idx_DS,3), '-', 'Color', color_R_Stance, 'LineWidth', 1.0);
        plot3(Traj_Right_Proj(1,1), Traj_Right_Proj(1,2), Traj_Right_Proj(1,3), 'o', 'MarkerFaceColor', color_R_Stance, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    elseif mode_idx == 2 % SS1 (Right)
        plot3(Traj_Right_Proj(idx_SS,1), Traj_Right_Proj(idx_SS,2), Traj_Right_Proj(idx_SS,3), '--', 'Color', color_R_Stance, 'LineWidth', 1.0);
        plot3(Traj_Right_Proj(idx_SS(1),1), Traj_Right_Proj(idx_SS(1),2), Traj_Right_Proj(idx_SS(1),3), 'o', 'MarkerFaceColor', color_R_Stance, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    elseif mode_idx == 3 % DS2 (Left)
        plot3(Traj_Left_Flipped_Proj(idx_DS,1), Traj_Left_Flipped_Proj(idx_DS,2), Traj_Left_Flipped_Proj(idx_DS,3), '-', 'Color', color_L_Stance_Flipped, 'LineWidth', 1.0);
        plot3(Traj_Left_Flipped_Proj(1,1), Traj_Left_Flipped_Proj(1,2), Traj_Left_Flipped_Proj(1,3), 's', 'MarkerFaceColor', color_L_Stance_Flipped, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    elseif mode_idx == 4 % SS2 (Left)
        plot3(Traj_Left_Flipped_Proj(idx_SS,1), Traj_Left_Flipped_Proj(idx_SS,2), Traj_Left_Flipped_Proj(idx_SS,3), '--', 'Color', color_L_Stance_Flipped, 'LineWidth', 1.0);
        plot3(Traj_Left_Flipped_Proj(idx_SS(1),1), Traj_Left_Flipped_Proj(idx_SS(1),2), Traj_Left_Flipped_Proj(idx_SS(1),3), 's', 'MarkerFaceColor', color_L_Stance_Flipped, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    end
    
    % --- 平面の描画 ---
    if mode_idx == 5
        planes_to_draw = 1:4;
    else
        planes_to_draw = target_p;
    end
    
    for p = planes_to_draw
        n = planes_def{p}.n; d_val = planes_def{p}.d;
        pts = [];
        for i = 1:size(edges, 1)
            p1 = corners(edges(i,1), :); p2 = corners(edges(i,2), :);
            vec = p2 - p1; denom = dot(n, vec);
            if abs(denom) > 1e-6
                t = (d_val - dot(n, p1)) / denom;
                if t >= -1e-5 && t <= 1 + 1e-5, pts = [pts; p1 + t * vec]; end
            end
        end
        if ~isempty(pts)
            pts = unique(round(pts, 5), 'rows');
            if size(pts, 1) > 2
                c = mean(pts, 1); [coeff_p, ~] = pca(pts - c);
                coords_2d = (pts - c) * coeff_p(:, 1:2);
                angles = atan2(coords_2d(:,2), coords_2d(:,1));
                [~, sort_idx] = sort(angles);
                pts = pts(sort_idx, :); pts(end+1, :) = pts(1, :);
                plot3(pts(:,1), pts(:,2), pts(:,3), planes_def{p}.style, ...
                      'Color', planes_def{p}.col, 'LineWidth', 0.5);
            end
        end
    end
    
    % --- 軸設定 (Fig 2と完全一致) ---
    xlim([box_min(1), box_max(1)]);
    ylim([box_min(2), box_max(2)]);
    zlim([box_min(3), box_max(3)]);
    
    xticks(ceil(box_min(1)/tick_step)*tick_step : tick_step : floor(box_max(1)/tick_step)*tick_step);
    yticks(ceil(box_min(2)/tick_step)*tick_step : tick_step : floor(box_max(2)/tick_step)*tick_step);
    zticks(ceil(box_min(3)/tick_step)*tick_step : tick_step : floor(box_max(3)/tick_step)*tick_step);
    
    xlabel('Synergy 1'); ylabel('Synergy 2'); zlabel('Synergy 3');
end


%% ===================================================================
% 12. 表データの出力 (Mean ± SD)
% ===================================================================
fprintf('\n===============================================================\n');
fprintf('  Statistical Table: Stride-by-Stride Analysis (Mean ± SD)\n');
fprintf('===============================================================\n');
fprintf('%-25s | %s\n', 'Metric', 'Value (Mean +/- SD)');
fprintf('---------------------------------------------------------------\n');

calc_stats = @(data) sprintf('%.2f ± %.2f', mean(data, 'omitnan'), std(data, 'omitnan'));
calc_stats_long = @(data) sprintf('%.4f ± %.4f', mean(data, 'omitnan'), std(data, 'omitnan'));

fprintf('lambda1^2/lambda2^2 (DS1) | %s\n', calc_stats(res_Ratio_R_DS)); 
fprintf('lambda1^2/lambda2^2 (SS1) | %s\n', calc_stats(res_Ratio_R_SS)); 
fprintf('lambda1^2/lambda2^2 (DS2) | %s\n', calc_stats(res_Ratio_L_DS)); 
fprintf('lambda1^2/lambda2^2 (SS2) | %s\n', calc_stats(res_Ratio_L_SS)); 

fprintf('---------------------------------------------------------------\n');

fprintf('Dot product (DS)          | %s\n', calc_stats_long(res_Sym_Dot_DS));
fprintf('Dot product (SS)          | %s\n', calc_stats_long(res_Sym_Dot_SS));

fprintf('===============================================================\n\n');

%% ===================================================================
% 13. データ保存
% ===================================================================
fprintf('Saving results to %s.mat...\n', tTextSaveMatName);

save(tTextSaveMatName, ...
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
    'indiv_R', 'indiv_L', ...
    'ratio_R', 'ratio_L', 'area_R', 'area_L', 'normal_dot_product', ...
    'res_Ratio_R', 'res_Ratio_L', 'res_Area_R', 'res_Area_L', ...
    'res_Planarity_R', 'res_Planarity_L', ... 
    'res_DotProd', 'res_Angle', ...
    'axis_labels_3D', 'cols_R_Leg', 'cols_L_Leg', ...
    'U_with_lambda', 'V_spatial', 'S_global', 'U_global', 'V_global', 'W_proj', ...
    'FLIP_MODES_GLOBAL', 'FLIP_MODES_DS1', 'FLIP_MODES_SS1', 'FLIP_MODES_DS2', 'FLIP_MODES_SS2', ...
    'V_DS1', 'S_DS1', 'U_DS1', 'normal_DS1', 'd_DS1', 'p_DS1_3D', 'P_DS1_Proj', ...
    'V_SS1', 'S_SS1', 'U_SS1', 'normal_SS1', 'd_SS1', 'p_SS1_3D', 'P_SS1_Proj', ...
    'V_DS2', 'S_DS2', 'U_DS2', 'normal_DS2', 'd_DS2', 'p_DS2_3D', 'P_DS2_Proj', ...
    'V_SS2', 'S_SS2', 'U_SS2', 'normal_SS2', 'd_SS2', 'p_SS2_3D', 'P_SS2_Proj', ...
    'Right_Half_Data', 'Left_Half_Data_Flipped', ...
    'Traj_Right_Proj', 'Traj_Left_Flipped_Proj', ...
    'n_DS2_F', 'd_DS2_F', 'n_SS2_F', 'd_SS2_F', ...
    'res_Sym_Dot_DS', 'res_Sym_Ang_DS', 'res_Sym_Dot_SS', 'res_Sym_Ang_SS', ...
    'res_Ratio_R_DS', 'res_Ratio_L_DS', 'res_Ratio_R_SS', 'res_Ratio_L_SS', ...
    'all_norm_L_DS1', 'all_norm_L_SS1');

fprintf('=== ALL PROCESS COMPLETED ===\n');

%% ===================================================================
% 14. グラフ画像の一括保存
% ===================================================================
if flg_graphSave
    fprintf('Saving figures using func_graphSave2 (Batch Mode)...\n');
    
    size_Fig1   = [800, 1000];
    size_Fig3D  = [700, 300];
    size_Fig2D  = [600, 200];
    size_Proj   = [250, 225];
    size_FigStat= [1000, 800];
    
    func_graphSave2(fig1,       tTextGraphName_Check_Angles, flg_graphSave, 0, size_Fig1);
    func_graphSave2(fig2,       tTextGraphName_Symmetry_3D,  flg_graphSave, 0, size_Fig3D);
    func_graphSave2(fig3,       tTextGraphName_Symmetry_2D,  flg_graphSave, 0, size_Fig2D);
    func_graphSave2(fig4,       tTextGraphName_Proj_DS,      flg_graphSave, 0, size_Proj);
    func_graphSave2(fig5,       tTextGraphName_Proj_SS,      flg_graphSave, 0, size_Proj);
    func_graphSave2(fig6,       tTextGraphName_Sym_Stats,    flg_graphSave, 0, size_FigStat);
    
    % 追加分の保存
    func_graphSave2(figure(7),  tTextGraphName_Decomp_DS1,   flg_graphSave, 0, size_Fig3D);
    func_graphSave2(figure(8),  tTextGraphName_Decomp_SS1,   flg_graphSave, 0, size_Fig3D);
    func_graphSave2(figure(9),  tTextGraphName_Decomp_DS2,   flg_graphSave, 0, size_Fig3D);
    func_graphSave2(figure(10), tTextGraphName_Decomp_SS2,   flg_graphSave, 0, size_Fig3D);
    func_graphSave2(figure(11), tTextGraphName_Planes_Only,  flg_graphSave, 0, size_Fig3D);
end

%% ===================================================================
% 15. 内部関数定義
% ===================================================================
function [normal, d] = calculate_plane_for_halfcycle(X_data, W_global)
    X_mean = mean(X_data, 1);
    X_var = X_data - X_mean;
    [~, ~, V_sub] = svd(X_var, 'econ');
    v1_3D = (V_sub(:,1)' * W_global);
    v2_3D = (V_sub(:,2)' * W_global);
    normal = cross(v1_3D, v2_3D); normal = normal / norm(normal);
    p_3D = X_mean * W_global;
    d = dot(normal, p_3D);
end