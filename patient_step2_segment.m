% =========================================================================
% 【プログラム概要: 片麻痺歩行用 セグメンテーション (Step 2)】
% 前処理済みデータ（Step 1）を読み込み、歩行イベント（接地: TD、離地: LO）を
% 基準にして1歩行周期（ストライド）ごとにデータを切り出し、時間正規化を行う。
% 
% 特徴:
% 1. データの解像度を上げるための補間処理（x5）を実行。
% 2. 右脚基準（TD_R -> TD_R）のストライド分割（SVDシナジー解析用）。
% 3. 左脚基準（TD_L -> TD_L）のストライド分割（Planar Lawの対称性評価用）。
% 4. 各ストライドを4つのフェーズ（DS1, SS1, DS2, SS2）に細分化。
% 5. 異常値（外れ値）のストライドを自動除外。
% =========================================================================

close all; 
clear;
clc;

%% =========================================================================
% 1. 初期設定（入出力ファイルと各種パラメータ）
% =========================================================================
SubjectName = 'KM_ID3'; % 被験者ID

% --- 入出力ファイル設定 ---
tTextLoadMatName = [SubjectName, '_step1_preprocessed']; % 読み込みファイル（Step 1の出力）
outputMatName    = [SubjectName, '_step2_segmented'];    % 本プログラムの保存先

% --- 補間処理（Interpolation）設定 ---
% フレーム間のデータ点数を擬似的に増やし、イベントタイミングのズレを減らす
interp_ratio  = 5;        % 補間倍率（5倍に引き伸ばす）
interp_method = 'pchip';  % 区分的3次エルミート補間（オーバーシュートを防ぐ手法）

% --- 時間正規化（Normalization）設定 ---
N_POINTS_GLOBAL = 200;   % 1ストライド全体をリサンプリングするデータ点数（SVD用）
N_POINTS_PHASE  = 100;   % 各歩行フェーズをリサンプリングするデータ点数

% --- 外れ値（Outlier）除去設定 ---
remove_outliers = true;  % 外れ値ストライドの除外フラグ
outlier_threshold = 2.5; % 中央値の2.5倍以上/以下を異常値として除外

% --- グラフ描画・保存設定 ---
show_stride_idx = 5;     % Fig 2で詳細を描画する代表ストライドの番号
flg_graphSave = 1;       % グラフ画像の保存フラグ

%% =========================================================================
% 2. データの読み込みと補間処理
% =========================================================================
fprintf('Loading data: "%s.mat" ...\n', tTextLoadMatName);
if exist([tTextLoadMatName, '.mat'], 'file')
    load([tTextLoadMatName, '.mat']); 
else
    error('ファイルが見つからない: %s.mat\nstep1を実行すること。', tTextLoadMatName);
end

% 仰角（theta）の変動成分（centeredtheta）が存在しない場合は再計算する
if ~exist('centeredtheta', 'var') || ~exist('mean_posture', 'var')
    fprintf('Warning: centeredtheta not found. Recalculating...\n');
    mean_posture = mean(theta, 1);
    centeredtheta = theta - mean_posture;
end

[n_frames, n_cols] = size(theta);

% 時間軸を interp_ratio 倍に引き伸ばし、各データを補間して滑らかにする
fprintf('Interpolating Data (x%d)...\n', interp_ratio);
t_original = 1:n_frames;
t_interp   = linspace(1, n_frames, n_frames * interp_ratio);

% PCA/SVD解析用（平均減算データ）の補間
theta_centered_interp = interp1(t_original, centeredtheta, t_interp, interp_method);

% Planar Lawなどの絶対角度表示用の補間
theta_interp = interp1(t_original, theta, t_interp, interp_method);

% 接地フラグ（0 or 1）の補間（線形補間後に0.5を閾値として二値化）
flags_temp = interp1(t_original, contact_flags, t_interp, 'linear');
contact_flags_interp = double(flags_temp > 0.5); 

% グラフ描画用の共通設定（Y軸のスケール統一など）
plotOrder = [7, 5, 3, 1, 2, 4, 6]; 
angleNames_SVD = {'L Foot', 'L Shank', 'L Thigh', 'Trunk', 'R Thigh', 'R Shank', 'R Foot'};
max_abs = max(abs(theta_centered_interp), [], 'all');
global_ylim = [-max_abs*1.1, max_abs*1.1]; 

%% =========================================================================
% 3. [Fig 1] 補間処理の結果チェック
% =========================================================================
% 補間によって波形が歪んでいないか、オリジナルデータと比較表示する。
fprintf('Generating Fig 1: Interpolation Check...\n');
fig1 = figure(1); clf(fig1);
set(fig1, 'Name', 'Fig1_InterpCheck', 'NumberTitle', 'off');
set(gcf, 'Position', [50, 50, 800, 1000]); 

for i = 1:7
    col_idx = plotOrder(i); 
    
    % オリジナルデータ
    subplot(7, 2, (i*2)-1); 
    plot(t_original, centeredtheta(:, col_idx), 'b.', 'MarkerSize', 3);
    if i == 1, title('Original'); end 
    ylabel([angleNames_SVD{i} ' (deg)']); grid on; xlim tight;
    
    % 補間後データ
    subplot(7, 2, i*2); 
    plot(t_interp, theta_centered_interp(:, col_idx), 'r.', 'MarkerSize', 1);
    if i == 1, title('Interpolated'); end
    grid on; xlim tight; ylabel('deg');
end
if flg_graphSave, func_graphSave2(fig1, [outputMatName, '_Fig1_InterpCheck'], flg_graphSave, 0, [960,600]); end

%% =========================================================================
% 4. [Right Stride] 右脚基準のストライド分割とフェーズ抽出
% =========================================================================
% 右足の接地(TD_R)から次の接地(TD_R)までを1ストライドとし、さらに左右の
% 離地・接地イベントを用いて4つの歩行フェーズに細分化する。
fprintf('Segmenting Right Strides (TD_R -> TD_R)...\n');
flag_L = contact_flags_interp(:, 1);
flag_R = contact_flags_interp(:, 2);

% 右足の接地イベント(0から1に切り替わる瞬間)を特定
TD_R_indices = find(diff([0; flag_R]) == 1);

num_strides_R_raw = length(TD_R_indices) - 1;
if num_strides_R_raw < 1, error('No Right strides found.'); end

% 切り出したストライドデータを格納する構造体の初期化
strides_R_accumulated = struct('raw_stride',{}, 'raw_stride_theta',{}, ...
    'raw_DS1',{}, 'raw_SS1',{}, 'raw_DS2',{}, 'raw_SS2',{}, ...
    'dur_Total',{}, 'dur_DS1',{}, 'dur_SS1',{}, 'dur_DS2',{}, 'dur_SS2',{}, ...
    'pct_LO',{});

valid_mask_R = true(num_strides_R_raw, 1);

for k = 1:num_strides_R_raw
    idx_start = TD_R_indices(k);
    idx_end   = TD_R_indices(k+1); 
    
    range_idxs = idx_start:idx_end;
    fL = flag_L(range_idxs); fR = flag_R(range_idxs);
    
    % --- フェーズ検出 ---
    % DS1 (両脚支持1): 右接地 ～ 左離地
    idx_L_LO_local = find(diff([1; fL]) == -1, 1, 'first');
    
    % SS1 (単脚支持1): 左離地 ～ 左接地
    if ~isempty(idx_L_LO_local)
        idx_L_TD_local = find(diff([0; fL(idx_L_LO_local:end)]) == 1, 1, 'first') + (idx_L_LO_local - 1);
    else, idx_L_TD_local = []; end
    
    % DS2 (両脚支持2): 左接地 ～ 右離地
    if ~isempty(idx_L_TD_local)
        idx_R_LO_local = find(diff([1; fR(idx_L_TD_local:end)]) == -1, 1, 'first') + (idx_L_TD_local - 1);
    else, idx_R_LO_local = []; end
    
    % イベントが一つでも欠損しているストライドは無効とする
    if isempty(idx_L_LO_local) || isempty(idx_L_TD_local) || isempty(idx_R_LO_local)
        valid_mask_R(k) = false; continue; 
    end
    
    idx_L_LO = idx_start + idx_L_LO_local - 1;
    idx_L_TD = idx_start + idx_L_TD_local - 1;
    idx_R_LO = idx_start + idx_R_LO_local - 1;
    
    % --- データ格納 ---
    st = struct();
    st.raw_stride = theta_centered_interp(idx_start:idx_end, :);     % SVD解析用 (平均減算済み)
    st.raw_stride_theta = theta_interp(idx_start:idx_end, :);        % Planar Law用 (絶対角度)
    
    st.raw_DS1    = theta_centered_interp(idx_start:idx_L_LO, :);
    st.raw_SS1    = theta_centered_interp(idx_L_LO:idx_L_TD, :);
    st.raw_DS2    = theta_centered_interp(idx_L_TD:idx_R_LO, :);
    st.raw_SS2    = theta_centered_interp(idx_R_LO:idx_end, :);
    
    st.dur_Total = size(st.raw_stride, 1);
    st.dur_DS1   = size(st.raw_DS1, 1);
    st.dur_SS1   = size(st.raw_SS1, 1);
    st.dur_DS2   = size(st.raw_DS2, 1);
    st.dur_SS2   = size(st.raw_SS2, 1);
    
    % Planar Law用: ストライド全体に対する右足離地(LO)のタイミング割合
    st.pct_LO = (idx_R_LO - idx_start) / st.dur_Total;
    
    strides_R_accumulated(k) = st;
end

% --- 外れ値の除去 (Right) ---
fprintf('Cleaning Right Strides...\n');
if remove_outliers
    v_idx = find(valid_mask_R);
    if ~isempty(v_idx)
        med_stride = median([strides_R_accumulated(v_idx).dur_Total]);
        med_DS1 = median([strides_R_accumulated(v_idx).dur_DS1]);
        med_SS1 = median([strides_R_accumulated(v_idx).dur_SS1]);
        med_DS2 = median([strides_R_accumulated(v_idx).dur_DS2]);
        med_SS2 = median([strides_R_accumulated(v_idx).dur_SS2]);
        
        % 中央値から大きく外れた時間長を持つフェーズを除外
        check_val = @(val, med) (val < med * (1/outlier_threshold)) || (val > med * outlier_threshold);
        
        for k = 1:num_strides_R_raw
            if valid_mask_R(k)
                is_bad = false;
                st = strides_R_accumulated(k);
                if check_val(st.dur_Total, med_stride), is_bad = true; end
                if check_val(st.dur_DS1, med_DS1), is_bad = true; end
                if check_val(st.dur_SS1, med_SS1), is_bad = true; end
                if check_val(st.dur_DS2, med_DS2), is_bad = true; end
                if check_val(st.dur_SS2, med_SS2), is_bad = true; end
                if is_bad, valid_mask_R(k) = false; end
            end
        end
    end
end
valid_indices_R = find(valid_mask_R);
valid_strides_R = strides_R_accumulated(valid_indices_R); 
num_valid_R = length(valid_indices_R);
fprintf('-> Valid Right Strides: %d\n', num_valid_R);

%% =========================================================================
% 5. [Left Stride] 左脚基準のストライド分割 (Planar Lawの対称性評価用)
% =========================================================================
% 歩行の非対称性を正確に評価するため、左脚の接地(TD_L)を基準としたストライドも別途抽出する。
fprintf('Segmenting Left Strides (TD_L -> TD_L)...\n');
TD_L_indices = find(diff([0; flag_L]) == 1);
LO_L_indices = find(diff([1; flag_L]) == -1);

strides_L_accumulated = struct('raw_stride',{}, 'raw_stride_theta',{}, ...
    'raw_DS1',{}, 'raw_SS1',{}, 'pct_LO', {}, 'dur_Total', {});
num_strides_L_raw = length(TD_L_indices) - 1;
valid_mask_L = true(num_strides_L_raw, 1);

if num_strides_L_raw > 0
    for k = 1:num_strides_L_raw
        idx_start = TD_L_indices(k);
        idx_end   = TD_L_indices(k+1);
        
        % ストライド内の左離地(LO_L)を特定
        cand_LO = LO_L_indices(LO_L_indices > idx_start & LO_L_indices < idx_end);
        
        if isempty(cand_LO)
            valid_mask_L(k) = false; continue;
        end
        idx_LO = cand_LO(1);
        
        st = struct();
        st.raw_stride = theta_centered_interp(idx_start:idx_end, :);
        st.raw_stride_theta = theta_interp(idx_start:idx_end, :); 
        st.dur_Total = idx_end - idx_start;
        st.pct_LO = (idx_LO - idx_start) / st.dur_Total;
        
        % Trajectory Map解析等で用いる簡易的なフェーズ分割（L_Stance = DS + SS）
        idx_TO_R_global = find(diff([1; flag_R]) == -1);
        cand_TO_R = idx_TO_R_global(idx_TO_R_global > idx_start & idx_TO_R_global < idx_end);
        idx_HS_R_global = find(diff([0; flag_R]) == 1);
        cand_HS_R = idx_HS_R_global(idx_HS_R_global > idx_start & idx_HS_R_global < idx_end);
        
        if ~isempty(cand_TO_R) && ~isempty(cand_HS_R)
             st.raw_DS1 = theta_centered_interp(idx_start:cand_TO_R(1), :);
             st.raw_SS1 = theta_centered_interp(cand_TO_R(1):cand_HS_R(1), :);
        else
             st.raw_DS1 = []; st.raw_SS1 = []; 
        end

        strides_L_accumulated(k) = st;
    end
    
    % --- 外れ値の除去 (Left) ---
    if remove_outliers
        v_idx_L = find(valid_mask_L);
        if ~isempty(v_idx_L)
            med_stride_L = median([strides_L_accumulated(v_idx_L).dur_Total]);
            check_val = @(val, med) (val < med * (1/outlier_threshold)) || (val > med * outlier_threshold);
            for k = 1:num_strides_L_raw
                if valid_mask_L(k)
                    if check_val(strides_L_accumulated(k).dur_Total, med_stride_L)
                        valid_mask_L(k) = false;
                    end
                end
            end
        end
    end
end
valid_indices_L = find(valid_mask_L);
valid_strides_L = strides_L_accumulated(valid_indices_L); 
fprintf('-> Valid Left Strides: %d\n', length(valid_strides_L));


%% =========================================================================
% 6. 時間正規化と平均化 (Right Stride基準)
% =========================================================================
% 抽出した各ストライドのデータ長を、統一したデータ点数（Global: 200, Phase: 100）にリサンプリングする。
all_normalized_strides = zeros(N_POINTS_GLOBAL, n_cols, num_valid_R);
all_normalized_strides_RAW = zeros(N_POINTS_GLOBAL, n_cols, num_valid_R); 

all_norm_DS1 = zeros(N_POINTS_PHASE, n_cols, num_valid_R);
all_norm_SS1 = zeros(N_POINTS_PHASE, n_cols, num_valid_R);
all_norm_DS2 = zeros(N_POINTS_PHASE, n_cols, num_valid_R);
all_norm_SS2 = zeros(N_POINTS_PHASE, n_cols, num_valid_R);

clean_raw_stride = cell(num_valid_R, 1);
clean_raw_DS1 = cell(num_valid_R, 1); clean_raw_SS1 = cell(num_valid_R, 1);
clean_raw_DS2 = cell(num_valid_R, 1); clean_raw_SS2 = cell(num_valid_R, 1);

dur_vec_Total = zeros(num_valid_R, 1);
dur_vec_DS1 = zeros(num_valid_R, 1); dur_vec_SS1 = zeros(num_valid_R, 1);
dur_vec_DS2 = zeros(num_valid_R, 1); dur_vec_SS2 = zeros(num_valid_R, 1); 

for i = 1:num_valid_R
    st = valid_strides_R(i);
    
    % ストライド全体 (平均減算済み)
    seg = st.raw_stride; clean_raw_stride{i} = seg;
    all_normalized_strides(:,:,i) = interp1(linspace(0,1,size(seg,1))', seg, linspace(0,1,N_POINTS_GLOBAL)', 'pchip');
    
    % ストライド全体 (絶対角度)
    seg_r = st.raw_stride_theta;
    all_normalized_strides_RAW(:,:,i) = interp1(linspace(0,1,size(seg_r,1))', seg_r, linspace(0,1,N_POINTS_GLOBAL)', 'pchip');
    
    % 各フェーズの時間正規化
    seg = st.raw_DS1; clean_raw_DS1{i} = seg;
    all_norm_DS1(:,:,i) = interp1(linspace(0,1,size(seg,1))', seg, linspace(0,1,N_POINTS_PHASE)', 'pchip');
    
    seg = st.raw_SS1; clean_raw_SS1{i} = seg;
    all_norm_SS1(:,:,i) = interp1(linspace(0,1,size(seg,1))', seg, linspace(0,1,N_POINTS_PHASE)', 'pchip');
    
    seg = st.raw_DS2; clean_raw_DS2{i} = seg;
    all_norm_DS2(:,:,i) = interp1(linspace(0,1,size(seg,1))', seg, linspace(0,1,N_POINTS_PHASE)', 'pchip');
    
    seg = st.raw_SS2; clean_raw_SS2{i} = seg;
    all_norm_SS2(:,:,i) = interp1(linspace(0,1,size(seg,1))', seg, linspace(0,1,N_POINTS_PHASE)', 'pchip');
    
    % フェーズ時間の記録
    dur_vec_Total(i) = st.dur_Total;
    dur_vec_DS1(i) = st.dur_DS1; dur_vec_SS1(i) = st.dur_SS1; 
    dur_vec_DS2(i) = st.dur_DS2; dur_vec_SS2(i) = st.dur_SS2; 
end

% 全ストライドの平均波形と標準偏差の算出
averaged_stride = mean(all_normalized_strides, 3);
std_stride      = std(all_normalized_strides, 0, 3); 

averaged_DS1 = mean(all_norm_DS1, 3); averaged_SS1 = mean(all_norm_SS1, 3);
averaged_DS2 = mean(all_norm_DS2, 3); averaged_SS2 = mean(all_norm_SS2, 3);

% 各フェーズが切り替わるタイミングの平均割合の算出
mean_pct_DS1_end = mean(dur_vec_DS1 ./ dur_vec_Total);
mean_pct_SS1_end = mean((dur_vec_DS1 + dur_vec_SS1) ./ dur_vec_Total);
mean_pct_DS2_end = mean((dur_vec_DS1 + dur_vec_SS1 + dur_vec_DS2) ./ dur_vec_Total);
mean_phase_transition_pcts = [mean_pct_DS1_end, mean_pct_SS1_end, mean_pct_DS2_end];

%% =========================================================================
% 7. [Fig 2] 単一ストライドの角度軌道詳細プロット
% =========================================================================
fprintf('Generating Fig 2: Single Stride Breakdown...\n');
if show_stride_idx > num_valid_R, target_i = num_valid_R; else, target_i = show_stride_idx; end
target_st = valid_strides_R(target_i);

len_Total = target_st.dur_Total;
pct_DS1 = (target_st.dur_DS1 / len_Total) * 100;
pct_SS1 = ((target_st.dur_DS1 + target_st.dur_SS1) / len_Total) * 100;
pct_DS2 = ((target_st.dur_DS1 + target_st.dur_SS1 + target_st.dur_DS2) / len_Total) * 100;
t_percent = linspace(0, 100, len_Total);

fig2 = figure(2); clf(fig2);
set(fig2, 'Name', 'Fig2_SingleStride', 'NumberTitle', 'off');
set(gcf, 'Position', [50, 50, 750, 800]); 

plot_pairs = {1, [], 'Trunk'; 3, 2, 'Thigh'; 5, 4, 'Shank'; 7, 6, 'Foot'};
for i = 1:4
    subplot(4, 1, i); hold on;
    idx_L = plot_pairs{i, 1}; idx_R = plot_pairs{i, 2};
    
    if isempty(idx_R)
        plot(t_percent, target_st.raw_stride(:, idx_L), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Trunk');
    else
        plot(t_percent, target_st.raw_stride(:, idx_R), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Right');
        plot(t_percent, target_st.raw_stride(:, idx_L), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Left');
    end
    
    % フェーズ境界線
    xline(pct_DS1, 'k-', 'HandleVisibility', 'off'); 
    xline(pct_SS1, 'k-', 'HandleVisibility', 'off'); 
    xline(pct_DS2, 'k-', 'HandleVisibility', 'off');
    
    ylim(global_ylim); xlim([0 100]); grid on; ylabel([plot_pairs{i, 3} ' [deg]']);
    
    if i == 4, lgd = legend('Location', 'southeast'); lgd.Position = [0.74, 0.02, 0.20, 0.05]; lgd.Box = 'off'; else, legend off; end
    if i == 1
        y_max = max(ylim);
        text(pct_DS1/2, y_max, 'DS1', 'Horiz','center', 'Vert','top');
        text((pct_DS1+pct_SS1)/2, y_max, 'SS1', 'Horiz','center', 'Vert','top');
        text((pct_SS1+pct_DS2)/2, y_max, 'DS2', 'Horiz','center', 'Vert','top');
        text((pct_DS2+100)/2, y_max, 'SS2', 'Horiz','center', 'Vert','top');
    end
    if i == 4, xlabel('Gait Cycle (%)'); else, xticklabels({}); end
    hold off;
end
if flg_graphSave, func_graphSave2(fig2, [outputMatName, '_Fig2_SingleStride'], flg_graphSave, 0, [480,300]); end

%% =========================================================================
% 8. [Fig 3~7] 各フェーズの平均角度軌道プロット
% =========================================================================
plots_config = {
    3, 'Whole Stride', N_POINTS_GLOBAL, clean_raw_stride, all_normalized_strides, averaged_stride, [400, 400];
    4, 'Phase: DS1',   N_POINTS_PHASE,  clean_raw_DS1,    all_norm_DS1,    averaged_DS1, [300, 400];
    5, 'Phase: SS1',   N_POINTS_PHASE,  clean_raw_SS1,    all_norm_SS1,    averaged_SS1, [300, 400];
    6, 'Phase: DS2',   N_POINTS_PHASE,  clean_raw_DS2,    all_norm_DS2,    averaged_DS2, [300, 400];
    7, 'Phase: SS2',   N_POINTS_PHASE,  clean_raw_SS2,    all_norm_SS2,    averaged_SS2, [300, 400];
};

for p = 1:5
    fig_num = plots_config{p, 1}; 
    title_str = plots_config{p, 2}; 
    n_pts = plots_config{p, 3};
    norm_d = plots_config{p, 5}; 
    avg_d = plots_config{p, 6};
    fig_size = plots_config{p, 7};
    
    std_d = std(norm_d, 0, 3); 
    x_ax = 1:n_pts;
    
    fprintf('Generating Fig %d: %s (Normalized Only)...\n', fig_num, title_str);
    
    f = figure(fig_num); clf(f);
    set(f, 'Name', ['Fig', num2str(fig_num)], 'NumberTitle', 'off');
    set(f, 'Position', [50 + (p-1)*20, 50, 450, 1000]); 
    
    for i = 1:7
        col = plotOrder(i); 
        
        subplot(7, 1, i); hold on;
        
        mu = avg_d(:, col)'; 
        sigma = std_d(:, col)';
        
        % 平均 ± 標準偏差の塗りつぶし描画
        fill([x_ax, fliplr(x_ax)], [mu+sigma, fliplr(mu-sigma)], ...
            [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        plot(x_ax, mu, 'k-', 'LineWidth', 2);
        
        hold off; 
        
        ylabel(angleNames_SVD{i}); 
        ylim(global_ylim); 
        grid on; 
        xlim([1 n_pts]); 
        
        if i < 7
            set(gca, 'XTickLabel', []); 
        else
            xlabel('Time (frame)');     
        end
        
        if i==7
            lgd = legend('Mean \pm SD', 'Mean', 'Location', 'southeast');
            lgd.Position = [0.60, 0.015, 0.35, 0.05]; 
            lgd.Box = 'off';
        end
    end
    
    if flg_graphSave
        func_graphSave2(f, [outputMatName, '_Fig', num2str(fig_num)], flg_graphSave, 0, fig_size);
    end
end

%% =========================================================================
% 9. [Fig 8] 各歩行フェーズの平均割合の帯グラフ
% =========================================================================
fprintf('Generating Fig 8: Phase Proportions...\n');
mean_vals = [mean(dur_vec_DS1), mean(dur_vec_SS1), mean(dur_vec_DS2), mean(dur_vec_SS2)];
phase_pct = (mean_vals / sum(mean_vals)) * 100; 

fig8 = figure(8); clf(fig8);
set(fig8, 'Name', 'Fig8_PhaseProp', 'NumberTitle', 'off');
set(gcf, 'Position', [300, 300, 800, 300]); 
b = barh(1, phase_pct, 'stacked', 'EdgeColor', 'none'); 

color_DS = [0.0, 0.45, 0.74]; color_SS = [0.6, 0.85, 1.0];   
b(1).FaceColor = color_DS; b(2).FaceColor = color_SS;
b(3).FaceColor = color_DS; b(4).FaceColor = color_SS;

phase_labels = {'DS1', 'SS1', 'DS2', 'SS2'}; cum_pct = [0, cumsum(phase_pct)];
for k = 1:4
    x_pos = (cum_pct(k) + cum_pct(k+1)) / 2;
    txt_col = 'black'; if mod(k, 2) == 1, txt_col = 'white'; end
    text(x_pos, 1, sprintf('%s\n%.1f%%', phase_labels{k}, phase_pct(k)), ...
        'HorizontalAlignment', 'center', 'Color', txt_col, 'FontWeight', 'bold', 'FontSize', 12);
end
xlabel('Gait Cycle [%]', 'FontSize', 12); xlim([0 100]); yticks([]); box off;    
if flg_graphSave, func_graphSave2(fig8, [outputMatName, '_Fig8_PhaseProp'], flg_graphSave, 0, [1200,400]); end

%% =========================================================================
% 10. データ保存
% =========================================================================
% 後続の解析ステップ（Planar Law, SVD）で必要となる全変数を保存する。
fprintf('Saving Data to "%s.mat" ...\n', outputMatName);

save(outputMatName, ...
    ... % [From Step 1]
    'SubjectName', 'samplingFrequency', 'data', 'rawData', ...
    'theta', 'centeredtheta', 'mean_posture', 'contact_flags', ...
    ...
    ... % [From Step 2: Main Results]
    'valid_strides_R', 'valid_strides_L', 'num_valid_R', ...
    'all_normalized_strides', 'all_normalized_strides_RAW', ...
    'averaged_stride', 'std_stride', ...
    ...
    ... % [From Step 2: Phase Data]
    'averaged_DS1', 'averaged_SS1', 'averaged_DS2', 'averaged_SS2', ...
    'all_norm_DS1', 'all_norm_SS1', 'all_norm_DS2', 'all_norm_SS2', ...
    'mean_phase_transition_pcts', ...
    ...
    ... % [From Step 2: Interpolated Data]
    'theta_centered_interp', 'theta_interp', 'contact_flags_interp');

fprintf('=== PROCESS COMPLETED ===\n');