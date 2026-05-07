close all; 
clear;     
clc;

% 機器設定ファイルの読み込み
make_status_bertec; 
make_status_mocap_29; 

%% =========================================================================
% 【プログラム概要: セグメンテーション (Segmentation & Normalization)】
% 前処理済みの連続データから、床反力（GRF）を基準にして1歩行周期（ストライド）
% ごとにデータを切り出し、時間正規化（0〜100% Gait Cycle）を行うプログラム。
% キネマティックシナジー解析（PCA/SVD）用の平均減算データと、
% 平面則（Planar Law）解析用の絶対角度データの両方を処理し保存する。
% =========================================================================

%% ===================================================================
% 1. 初期設定（入出力ファイルとパラメータ）
% ===================================================================
SubjectName = 'IBA';

% --- ファイル設定 ---
tTextLoadMatName = [SubjectName, '_12_step1_preprocessed']; % 読み込みファイル
tTextSaveMatName = [SubjectName, '_12_step2_segmented'];    % 最終保存ファイル

% --- グラフ保存設定 ---
tTextGraphName_SingleStride= [SubjectName, '_Seg_Fig1_SingleStride']; 
tTextGraphName_Check_Angle = [SubjectName, '_Seg_Fig2_AllAngles'];          
tTextGraphName_DS1         = [SubjectName, '_Seg_Fig3_DS1'];
tTextGraphName_SS1         = [SubjectName, '_Seg_Fig4_SS1'];
tTextGraphName_DS2         = [SubjectName, '_Seg_Fig5_DS2'];
tTextGraphName_SS2         = [SubjectName, '_Seg_Fig6_SS2'];
tTextGraphName_GRF_Check   = [SubjectName, '_Seg_Fig7_GRF_Check'];
tTextGraphName_PhaseProp   = [SubjectName, '_Seg_Fig8_PhaseProp'];

flg_graphSave = 1; % グラフ保存フラグ

% --- セグメンテーション用パラメータ ---
THRESHOLD_PERCENT = 0.07; % 接地(HS)/離地(TO)を判定する床反力の閾値（体重の7%）
N_POINTS_GLOBAL   = 200;  % 1ストライド全体をリサンプリングするデータ点数
N_POINTS_PHASE    = 100;  % 各歩行フェーズ（DS1, SS1など）をリサンプリングするデータ点数
REMOVE_OUTLIERS   = true; % 外れ値ストライドの除外フラグ
OUTLIER_THRESHOLD = 1.5;  % 外れ値判定の閾値（中央値の1.5倍以上/以下を除外）
Target_Stride_Index = 2;  % Fig1で詳細表示する代表ストライドの番号

GRF_CH_R = eNUM_Fz1; % 右足の鉛直方向床反力チャンネル
GRF_CH_L = eNUM_Fz2; % 左足の鉛直方向床反力チャンネル

% 体重に基づく閾値(N)の算出
if exist('common', 'var') && isfield(common, 'subject') && isfield(common.subject', 'WeightN')
    weight_val = common.subject.WeightN;
else
    weight_val = 60 * 9.8; 
    fprintf('Warning: common data not found. Using default weight 60kg.\n');
end
threshold_N = weight_val * THRESHOLD_PERCENT;

fprintf('--- Segmentation Parameters ---\n');
fprintf('Threshold: %.1f N\n', threshold_N);

%% ===================================================================
% 2. データの読込と全試行の結合処理
% ===================================================================
fprintf('\nLoading Merged Data: %s.mat...\n', tTextLoadMatName);

if ~isfile([tTextLoadMatName, '.mat'])
    error('File not found: %s. Please run preprocess.m first.', tTextLoadMatName);
end

% 前処理データのロード
loadedStruct = load(tTextLoadMatName);
mergedData = loadedStruct.mergedData; 

% 複数試行のデータを縦に結合し、1つの連続データとして扱う
fprintf('Concatenating data from all trials...\n');
theta = vertcat(mergedData.theta);                      % 全試行の絶対角度 (Planar Law用)
centeredtheta = vertcat(mergedData.centeredtheta);      % 全試行の変動成分 (シナジー抽出用)
meanposture = mean(vertcat(mergedData.meanposture), 1); % 平均姿勢

data_combined = struct();
force_cell = arrayfun(@(x) x.data.force, mergedData, 'UniformOutput', false);
mocap_cell = arrayfun(@(x) x.data.mocap, mergedData, 'UniformOutput', false);
data_combined.force = vertcat(force_cell{:});
data_combined.mocap = vertcat(mocap_cell{:});
data = data_combined; 

%% ===================================================================
% 3. 歩行イベントの検出とストライド抽出（右足基準）
% ===================================================================
% ここでは「右足接地(HS) 〜 次の右足接地(HS)」を1ストライドとして切り出す。
strides_accumulated = [];
raw_strides_total_count = 0;

for fIdx = 1:length(mergedData)
    
    CurrentTrial = mergedData(fIdx).trialName;
    fprintf('Processing Trial (Right Base): %s...\n', CurrentTrial);
    
    force_data_local    = mergedData(fIdx).data.force;
    centeredtheta_local = mergedData(fIdx).centeredtheta;
    theta_local         = mergedData(fIdx).theta; 
    
    % --- イベント検出 (床反力ベース) ---
    % 閾値を超えた瞬間をHeel Strike(HS)、下回った瞬間をToe Off(TO)とする
    is_stance_R = force_data_local(:, GRF_CH_R) > threshold_N;
    ev_HS_R = find(diff(is_stance_R) == 1) + 1;
    ev_TO_R = find(diff(is_stance_R) == -1) + 1;

    is_stance_L = force_data_local(:, GRF_CH_L) > threshold_N;
    ev_HS_L = find(diff(is_stance_L) == 1) + 1;
    ev_TO_L = find(diff(is_stance_L) == -1) + 1;

    % --- ストライド抽出 ---
    local_strides = struct('raw_stride',{}, 'raw_stride_theta', {}, ... 
                          'raw_DS1',{}, 'raw_SS1',{}, 'raw_DS2',{}, 'raw_SS2',{}, ...
                          'dur_Total',{}, 'dur_DS1',{}, 'dur_SS1',{}, 'dur_DS2',{}, 'dur_SS2',{}, ...
                          'idx_Start',{}, 'source_trial',{}, 'raw_grf_R', {}, 'raw_grf_L', {}, 'pct_LO', {});
    local_count = 0;

    for i = 1:length(ev_HS_R)-1
        idx_HS_R_start = ev_HS_R(i);   % 1. 右接地
        idx_HS_R_next  = ev_HS_R(i+1); % 5. 次の右接地
        
        % ストライド内の各イベントを時系列順に特定
        cand_TO_L = ev_TO_L(ev_TO_L > idx_HS_R_start & ev_TO_L < idx_HS_R_next);
        if isempty(cand_TO_L), continue; end; idx_TO_L = cand_TO_L(1); % 2. 左離地
        
        cand_HS_L = ev_HS_L(ev_HS_L > idx_TO_L & ev_HS_L < idx_HS_R_next);
        if isempty(cand_HS_L), continue; end; idx_HS_L = cand_HS_L(1); % 3. 左接地
        
        cand_TO_R = ev_TO_R(ev_TO_R > idx_HS_L & ev_TO_R < idx_HS_R_next);
        if isempty(cand_TO_R), continue; end; idx_TO_R = cand_TO_R(1); % 4. 右離地
        
        local_count = local_count + 1;
        
        % 角度データの保存
        local_strides(local_count).raw_stride       = centeredtheta_local(idx_HS_R_start:idx_HS_R_next, :); % PCA用
        local_strides(local_count).raw_stride_theta = theta_local(idx_HS_R_start:idx_HS_R_next, :);         % Planar Law用
        
        % 4つの歩行フェーズ（両脚支持・単脚支持）に分割
        local_strides(local_count).raw_DS1 = centeredtheta_local(idx_HS_R_start:idx_TO_L, :);
        local_strides(local_count).raw_SS1 = centeredtheta_local(idx_TO_L:idx_HS_L, :);
        local_strides(local_count).raw_DS2 = centeredtheta_local(idx_HS_L:idx_TO_R, :);
        local_strides(local_count).raw_SS2 = centeredtheta_local(idx_TO_R:idx_HS_R_next, :);
        
        % 各フェーズのフレーム数（時間長）を記録
        local_strides(local_count).dur_Total = idx_HS_R_next - idx_HS_R_start;
        local_strides(local_count).dur_DS1   = idx_TO_L - idx_HS_R_start;
        local_strides(local_count).dur_SS1   = idx_HS_L - idx_TO_L;
        local_strides(local_count).dur_DS2   = idx_TO_R - idx_HS_L;
        local_strides(local_count).dur_SS2   = idx_HS_R_next - idx_TO_R;
        
        % 属性情報の保存
        local_strides(local_count).idx_Start    = idx_HS_R_start;
        local_strides(local_count).source_trial = CurrentTrial;
        local_strides(local_count).raw_grf_R    = force_data_local(idx_HS_R_start:idx_HS_R_next, GRF_CH_R);
        local_strides(local_count).raw_grf_L    = force_data_local(idx_HS_R_start:idx_HS_R_next, GRF_CH_L);
        
        % Planar Law用：右足基準における右足離地(LO)のタイミング (0〜1の割合)
        local_strides(local_count).pct_LO = (idx_TO_R - idx_HS_R_start) / local_strides(local_count).dur_Total;
    end
    
    if local_count > 0
        strides_accumulated = [strides_accumulated, local_strides];
    end
    raw_strides_total_count = raw_strides_total_count + local_count;
end

strides_temp = strides_accumulated;
raw_strides_count = raw_strides_total_count;

if raw_strides_count == 0, error('No strides found.'); end
fprintf('\nTotal Right Strides Loaded: %d\n', raw_strides_count);

%% ===================================================================
% 4. ストライド選別（外れ値の除去）
% ===================================================================
% 旋回時やノイズによる異常な長さのストライドを中央値を利用して除外する。
valid_mask = true(raw_strides_count, 1);
if REMOVE_OUTLIERS
    med_total = median([strides_temp.dur_Total]);
    med_DS1   = median([strides_temp.dur_DS1]); 
    med_SS1   = median([strides_temp.dur_SS1]);
    med_DS2   = median([strides_temp.dur_DS2]); 
    med_SS2   = median([strides_temp.dur_SS2]);
    
    check_val = @(val, med) (val < med * (1/OUTLIER_THRESHOLD)) || (val > med * OUTLIER_THRESHOLD);
    
    for k = 1:raw_strides_count
        is_bad = false;
        if check_val(strides_temp(k).dur_Total, med_total), is_bad = true; end
        if check_val(strides_temp(k).dur_DS1,   med_DS1),   is_bad = true; end
        if check_val(strides_temp(k).dur_SS1,   med_SS1),   is_bad = true; end
        if check_val(strides_temp(k).dur_DS2,   med_DS2),   is_bad = true; end
        if check_val(strides_temp(k).dur_SS2,   med_SS2),   is_bad = true; end
        if is_bad, valid_mask(k) = false; end
    end
end
valid_indices = find(valid_mask);
valid_stride_count = length(valid_indices);
fprintf('-> Valid Right Strides: %d\n', valid_stride_count);
if valid_stride_count == 0, error('No valid strides left.'); end

%% ===================================================================
% 5. 時間正規化処理 (100% Gait Cycle へのリサンプリング)
% ===================================================================
% 全てのストライドを同じデータ点数（N_POINTS_GLOBAL）に補間して長さを揃える。
[~, N_angles] = size(strides_temp(1).raw_stride); 

% 初期化
all_norm_Global     = zeros(N_POINTS_GLOBAL, N_angles, valid_stride_count); % PCA用
all_norm_Global_RAW = zeros(N_POINTS_GLOBAL, N_angles, valid_stride_count); % Planar Law用

all_norm_DS1 = zeros(N_POINTS_PHASE, N_angles, valid_stride_count);
all_norm_SS1 = zeros(N_POINTS_PHASE, N_angles, valid_stride_count);
all_norm_DS2 = zeros(N_POINTS_PHASE, N_angles, valid_stride_count);
all_norm_SS2 = zeros(N_POINTS_PHASE, N_angles, valid_stride_count);

clean_raw_Full = cell(valid_stride_count, 1);
clean_raw_DS1  = cell(valid_stride_count, 1);
clean_raw_SS1  = cell(valid_stride_count, 1);
clean_raw_DS2  = cell(valid_stride_count, 1);
clean_raw_SS2  = cell(valid_stride_count, 1);

time_vec_global = linspace(0, 1, N_POINTS_GLOBAL)';
time_vec_phase  = linspace(0, 1, N_POINTS_PHASE)';

% フェーズ比率計算用配列
dur_vec_Total = zeros(valid_stride_count, 1);
dur_vec_DS1   = zeros(valid_stride_count, 1);
dur_vec_SS1   = zeros(valid_stride_count, 1);
dur_vec_DS2   = zeros(valid_stride_count, 1);

for i = 1:valid_stride_count
    st = strides_temp(valid_indices(i));
    
    clean_raw_Full{i} = st.raw_stride;
    clean_raw_DS1{i}  = st.raw_DS1; 
    clean_raw_SS1{i}  = st.raw_SS1;
    clean_raw_DS2{i}  = st.raw_DS2; 
    clean_raw_SS2{i}  = st.raw_SS2;
    
    % 正規化 (pchip: 区分的3次エルミート補間)
    all_norm_Global(:,:,i)     = interp1(linspace(0,1,size(st.raw_stride,1)), st.raw_stride, time_vec_global, 'pchip');
    all_norm_Global_RAW(:,:,i) = interp1(linspace(0,1,size(st.raw_stride_theta,1)), st.raw_stride_theta, time_vec_global, 'pchip');
    
    % フェーズごとの正規化
    all_norm_DS1(:,:,i) = interp1(linspace(0,1,size(st.raw_DS1,1)), st.raw_DS1, time_vec_phase, 'pchip');
    all_norm_SS1(:,:,i) = interp1(linspace(0,1,size(st.raw_SS1,1)), st.raw_SS1, time_vec_phase, 'pchip');
    all_norm_DS2(:,:,i) = interp1(linspace(0,1,size(st.raw_DS2,1)), st.raw_DS2, time_vec_phase, 'pchip');
    all_norm_SS2(:,:,i) = interp1(linspace(0,1,size(st.raw_SS2,1)), st.raw_SS2, time_vec_phase, 'pchip');

    % フェーズ時間の記録
    dur_vec_Total(i) = st.dur_Total;
    dur_vec_DS1(i)   = st.dur_DS1;
    dur_vec_SS1(i)   = st.dur_SS1;
    dur_vec_DS2(i)   = st.dur_DS2;
end

% 統計量の計算（平均と標準偏差）
averaged_stride_global = mean(all_norm_Global, 3); std_stride_global = std(all_norm_Global, 0, 3);
averaged_DS1 = mean(all_norm_DS1, 3); std_DS1 = std(all_norm_DS1, 0, 3);
averaged_SS1 = mean(all_norm_SS1, 3); std_SS1 = std(all_norm_SS1, 0, 3);
averaged_DS2 = mean(all_norm_DS2, 3); std_DS2 = std(all_norm_DS2, 0, 3);
averaged_SS2 = mean(all_norm_SS2, 3); std_SS2 = std(all_norm_SS2, 0, 3);

% GRF平均用ダミー変数（将来的な拡張用）
averaged_grf_R = zeros(400, 1); averaged_grf_L = zeros(400, 1);
std_grf_R = zeros(400, 1); std_grf_L = zeros(400, 1);
averaged_grf = averaged_grf_R; 

% --- フェーズ遷移タイミングの平均 (0.0 - 1.0) の計算 ---
mean_pct_DS1_end = mean(dur_vec_DS1 ./ dur_vec_Total);                               % Left TO
mean_pct_SS1_end = mean((dur_vec_DS1 + dur_vec_SS1) ./ dur_vec_Total);               % Left HS
mean_pct_DS2_end = mean((dur_vec_DS1 + dur_vec_SS1 + dur_vec_DS2) ./ dur_vec_Total); % Right TO

% [Left_TO, Left_HS, Right_TO] のタイミング比率
mean_phase_transition_pcts = [mean_pct_DS1_end, mean_pct_SS1_end, mean_pct_DS2_end];
fprintf('Mean Phase Transitions: L_TO=%.1f%%, L_HS=%.1f%%, R_TO=%.1f%%\n', ...
        mean_phase_transition_pcts(1)*100, mean_phase_transition_pcts(2)*100, mean_phase_transition_pcts(3)*100);

%% ===================================================================
% 6. グラフ作成の共通設定
% ===================================================================
plotOrder = [7, 5, 3, 1, 2, 4, 6]; 
angleNames_SVD = {'L Foot', 'L Shank', 'L Thigh', 'Trunk', 'R Thigh', 'R Shank', 'R Foot'};
max_abs = max(abs(averaged_stride_global), [], 'all'); 
global_ylim = [-max_abs*1.1, max_abs*1.1]; % Y軸の統一

%% ===================================================================
% 7. 各種グラフの生成と保存
% ===================================================================

% --- [Fig 1] 単一ストライドの角度軌道詳細 ---
fprintf('Generating Fig 1: Single Stride Breakdown...\n');
if Target_Stride_Index > valid_stride_count
    display_idx_in_valid = valid_stride_count;
else
    display_idx_in_valid = Target_Stride_Index; 
end
target_st = strides_temp(valid_indices(display_idx_in_valid));

len_Total = target_st.dur_Total;
pct_DS1 = (target_st.dur_DS1 / len_Total) * 100;
pct_SS1 = ((target_st.dur_DS1 + target_st.dur_SS1) / len_Total) * 100;
pct_DS2 = ((target_st.dur_DS1 + target_st.dur_SS1 + target_st.dur_DS2) / len_Total) * 100;
t_percent = linspace(0, 100, size(target_st.raw_stride, 1));

fig1 = figure('Name', 'Fig1_SingleStride', 'Position', [50, 50, 750, 800]);
plot_pairs = { 1, [], 'Trunk'; 3, 2, 'Thigh'; 5, 4, 'Shank'; 7, 6, 'Foot' };

for i = 1:4
    subplot(4, 1, i); hold on;
    idx_L = plot_pairs{i, 1}; idx_R = plot_pairs{i, 2}; title_str = plot_pairs{i, 3};
    
    if isempty(idx_R)
        h_trunk = plot(t_percent, target_st.raw_stride(:, idx_L), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Trunk');
        h_leg = []; 
    else
        h_R = plot(t_percent, target_st.raw_stride(:, idx_R), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Right');
        h_L = plot(t_percent, target_st.raw_stride(:, idx_L), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Left');
        h_leg = [h_R, h_L];
    end
    
    xline(pct_DS1, 'k-'); xline(pct_SS1, 'k-'); xline(pct_DS2, 'k-');
    ylim(global_ylim); xlim([0 100]); grid on; ylabel([title_str ' (deg)']);
    
    if i == 4
        if isempty(h_leg)
            lgd = legend(h_trunk, 'Location', 'southeast');
        else
            lgd = legend(h_leg, 'Location', 'southeast');
        end
        lgd.Position = [0.74, 0.02, 0.20, 0.05]; lgd.Box = 'off';
    else
        legend off;
    end
    
    if i == 1
        y_max = max(ylim);
        text(pct_DS1/2, y_max, 'DS1', 'Horiz','center', 'Vert','top');
        text((pct_DS1+pct_SS1)/2, y_max, 'SS1', 'Horiz','center', 'Vert','top');
        text((pct_SS1+pct_DS2)/2, y_max, 'DS2', 'Horiz','center', 'Vert','top');
        text((pct_DS2+100)/2, y_max, 'SS2', 'Horiz','center', 'Vert','top');
    end
    
    if i < 4
        set(gca, 'XTickLabel', []); 
    else
        xlabel('Gait Cycle (%)');     
    end
    hold off;
end
if flg_graphSave, func_graphSave2(fig1, tTextGraphName_SingleStride, flg_graphSave, 0, [480,300]); end

% --- [Fig 2] 全7体節の平均角度軌道 ---
fprintf('Generating Fig 2: All Angles (Normalized Only)...\n');
fig2 = figure('Name', 'Fig2_AllAngles_Global', 'Position', [100, 50, 450, 1000]); 
plot_x_global = 1:N_POINTS_GLOBAL;

for j = 1:7
    col_idx = plotOrder(j);
    subplot(7, 1, j); hold on;
    mu = averaged_stride_global(:, col_idx)'; 
    sigma = std_stride_global(:, col_idx)';
    
    fill([plot_x_global, fliplr(plot_x_global)], [mu+sigma, fliplr(mu-sigma)], ...
        [0.8 0.8 0.8], 'EdgeColor','none','FaceAlpha',0.6);
    plot(plot_x_global, mu, 'k-', 'LineWidth', 2);
    
    hold off; 
    ylabel(angleNames_SVD{j}); 
    ylim(global_ylim); xlim([1 N_POINTS_GLOBAL]); grid on; 
    
    if j < 7
        set(gca, 'XTickLabel', []); 
    else
        xlabel('Time (frame)');     
    end

    if j==7
        lgd = legend('Mean \pm SD', 'Mean', 'Location', 'southeast');
        lgd.Position = [0.60, 0.015, 0.35, 0.05]; 
        lgd.Box = 'off';
    end
end
if flg_graphSave, func_graphSave2(fig2, tTextGraphName_Check_Angle, flg_graphSave, 0, [400,400]); end

% --- [Fig 3~6] 各フェーズ(DS1/SS1/DS2/SS2)ごとの角度軌道 ---
figs_config = {
    3, 'DS1', clean_raw_DS1, all_norm_DS1, averaged_DS1, std_DS1, tTextGraphName_DS1;
    4, 'SS1', clean_raw_SS1, all_norm_SS1, averaged_SS1, std_SS1, tTextGraphName_SS1;
    5, 'DS2', clean_raw_DS2, all_norm_DS2, averaged_DS2, std_DS2, tTextGraphName_DS2;
    6, 'SS2', clean_raw_SS2, all_norm_SS2, averaged_SS2, std_SS2, tTextGraphName_SS2;
};
plot_x_phase = 1:N_POINTS_PHASE;

for p = 1:4
    f_num  = figs_config{p, 1};
    p_name = figs_config{p, 2};
    avg_d  = figs_config{p, 5};
    std_d  = figs_config{p, 6};
    s_name = figs_config{p, 7};
    
    ff = figure(f_num); set(ff, 'Name', ['Phase ' p_name], 'Position', [100+p*20, 50, 450, 1000]);
    
    for j = 1:7
        col_idx = plotOrder(j);
        subplot(7, 1, j); hold on;
        mu = avg_d(:, col_idx)'; 
        sigma = std_d(:, col_idx)';
        
        fill([plot_x_phase, fliplr(plot_x_phase)], [mu+sigma, fliplr(mu-sigma)], ...
            [0.8 0.8 0.8], 'EdgeColor','none','FaceAlpha',0.6);
        plot(plot_x_phase, mu, 'k-', 'LineWidth', 2);
        
        hold off; 
        ylabel(angleNames_SVD{j}); 
        ylim(global_ylim); xlim([1 N_POINTS_PHASE]); grid on;
        
        if j < 7
            set(gca, 'XTickLabel', []); 
        else
            xlabel('Time (frame)');     
        end
        
        if j==7
            lgd = legend('Mean \pm SD', 'Mean', 'Location', 'southeast');
            lgd.Position = [0.60, 0.015, 0.35, 0.05];
            lgd.Box = 'off';
        end
    end
    if flg_graphSave, func_graphSave2(ff, s_name, flg_graphSave, 0, [300,400]); end
end

% --- [Fig 7] 床反力によるセグメンテーション結果の確認 ---
fprintf('Generating Fig 7: GRF Segmentation Check (Stride #%d)...\n', display_idx_in_valid);

if isfield(target_st, 'raw_grf_R')
    raw_F_R = target_st.raw_grf_R;
    raw_F_L = target_st.raw_grf_L;
    
    grf_len = length(raw_F_R);
    t_norm_grf = linspace(0, 100, grf_len);
    
    p_DS1_end = (target_st.dur_DS1 / target_st.dur_Total) * 100;
    p_SS1_end = ((target_st.dur_DS1 + target_st.dur_SS1) / target_st.dur_Total) * 100;
    p_DS2_end = ((target_st.dur_DS1 + target_st.dur_SS1 + target_st.dur_DS2) / target_st.dur_Total) * 100;
    
    fig7 = figure('Name', 'Fig7_GRF_Segmentation', 'Position', [150, 150, 900, 500]);
    hold on;
    
    y_min = -50; 
    y_max = max([max(raw_F_R), max(raw_F_L)]) * 1.1;
    ylim([y_min, y_max]);
    
    fill([0, p_DS1_end, p_DS1_end, 0], [y_min, y_min, y_max, y_max], ...
        [0.8 0.9 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'DS1');
    fill([p_DS1_end, p_SS1_end, p_SS1_end, p_DS1_end], [y_min, y_min, y_max, y_max], ...
        [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'SS1');
    fill([p_SS1_end, p_DS2_end, p_DS2_end, p_SS1_end], [y_min, y_min, y_max, y_max], ...
        [0.8 1 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'DS2');
    fill([p_DS2_end, 100, 100, p_DS2_end], [y_min, y_min, y_max, y_max], ...
        [0.8 0.9 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'SS2');
    
    yline(threshold_N, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Threshold');
    
    plot(t_norm_grf, raw_F_R, 'r-', 'LineWidth', 2.0, 'DisplayName', 'Right GRF');
    plot(t_norm_grf, raw_F_L, 'b-', 'LineWidth', 2.0, 'DisplayName', 'Left GRF');
    
    xlabel('Gait Cycle (%)'); ylabel('Force (N)'); xlim([0, 100]); grid on;
    title(['GRF Check: ' target_st.source_trial], 'Interpreter', 'none');
    legend('Location', 'best');
    hold off;
    
    if flg_graphSave, func_graphSave2(fig7, tTextGraphName_GRF_Check, flg_graphSave, 0, [800, 550]); end
else
    warning('GRF data missing in stride struct. Cannot generate Fig 7.');
end

% --- [Fig 8] 各フェーズの平均割合（帯グラフ） ---
fprintf('Generating Fig 8: Average Phase Proportions...\n');

vals_DS1 = [strides_temp(valid_indices).dur_DS1];
vals_SS1 = [strides_temp(valid_indices).dur_SS1];
vals_DS2 = [strides_temp(valid_indices).dur_DS2];
vals_SS2 = [strides_temp(valid_indices).dur_SS2];

mean_vals = [mean(vals_DS1), mean(vals_SS1), mean(vals_DS2), mean(vals_SS2)];
phase_pct = (mean_vals / sum(mean_vals)) * 100;

fig8 = figure('Name', 'Fig8_PhaseProportions', 'Position', [300, 300, 800, 300]); 
b = barh(1, phase_pct, 'stacked', 'EdgeColor', 'none');

color_DS = [0.0, 0.45, 0.74];  
color_SS = [0.6, 0.85, 1.0];   
b(1).FaceColor = color_DS; b(2).FaceColor = color_SS;
b(3).FaceColor = color_DS; b(4).FaceColor = color_SS;

phase_labels = {'DS1', 'SS1', 'DS2', 'SS2'};
cum_pct = [0, cumsum(phase_pct)];

for k = 1:4
    x_pos = (cum_pct(k) + cum_pct(k+1)) / 2;
    text(x_pos, 1, sprintf('%s\n%.1f%%', phase_labels{k}, phase_pct(k)), ...
        'Horizontal', 'center', 'Color', 'white', 'FontWeight', 'bold', 'FontSize', 12);
    if mod(k, 2) == 0, text(x_pos, 1, ...
        sprintf('%s\n%.1f%%', phase_labels{k}, phase_pct(k)), ...
        'Horizontal', 'center', 'Color', 'black', 'FontWeight', 'bold', 'FontSize', 12);
    end
end
xlabel('Gait Cycle [%]', 'FontSize', 12); xlim([0 100]); yticks([]); box off;    
if flg_graphSave
    func_graphSave2(fig8, tTextGraphName_PhaseProp, flg_graphSave, 0, [1200,400]);
end

%% ===================================================================
% 8. Planar Law解析用：左足基準(HS_L -> HS_L)のストライド抽出
% ===================================================================
% 歩行の非対称性を考慮し、左足を基準としたストライドも抽出する。
fprintf('\n--- Generating Left Leg Based Strides (HS_L -> HS_L) for Planar Law ---\n');

strides_L_accumulated = [];

for fIdx = 1:length(mergedData)
    force_local  = mergedData(fIdx).data.force;
    ctheta_local = mergedData(fIdx).centeredtheta;
    theta_local  = mergedData(fIdx).theta;
    
    is_stance_L = force_local(:, GRF_CH_L) > threshold_N;
    ev_HS_L = find(diff(is_stance_L) == 1) + 1;
    ev_TO_L = find(diff(is_stance_L) == -1) + 1;
    
    % ストライド抽出 (HS_L -> HS_L)
    for i = 1:length(ev_HS_L)-1
        idx_Start = ev_HS_L(i);
        idx_End   = ev_HS_L(i+1);
        
        cand_TO_L = ev_TO_L(ev_TO_L > idx_Start & ev_TO_L < idx_End);
        if isempty(cand_TO_L), continue; end
        idx_TO_L_in_stride = cand_TO_L(1);
        
        st = struct();
        st.raw_stride       = ctheta_local(idx_Start:idx_End, :);
        st.raw_stride_theta = theta_local(idx_Start:idx_End, :);
        st.dur_Total        = idx_End - idx_Start;
        st.source_trial     = mergedData(fIdx).trialName;
        
        % 左足基準におけるSwing開始点 (TO_L) の相対位置
        st.pct_LO = (idx_TO_L_in_stride - idx_Start) / st.dur_Total;
        
        if isempty(strides_L_accumulated)
            strides_L_accumulated = st;
        else
            strides_L_accumulated(end+1) = st;
        end
    end
end

% 左足基準データの外れ値除去
valid_strides_L = [];
if ~isempty(strides_L_accumulated)
    valid_mask_L = true(length(strides_L_accumulated), 1);
    if REMOVE_OUTLIERS
        med_total_L = median([strides_L_accumulated.dur_Total]);
        check_val = @(val, med) (val < med * (1/OUTLIER_THRESHOLD)) || (val > med * OUTLIER_THRESHOLD);
        for k = 1:length(strides_L_accumulated)
            if check_val(strides_L_accumulated(k).dur_Total, med_total_L)
                valid_mask_L(k) = false;
            end
        end
    end
    valid_strides_L = strides_L_accumulated(valid_mask_L);
end

fprintf('Valid Left Strides Extracted: %d\n', length(valid_strides_L));
valid_strides_R = strides_temp(valid_indices);

%% ===================================================================
% 9. データの保存
% ===================================================================
% 後続のシナジー解析（PCA/SVD）や平面則解析で利用する全ての変数を保存する。
averaged_stride            = averaged_stride_global;
std_stride                 = std_stride_global;
all_normalized_strides     = all_norm_Global;     % シナジー解析用（平均減算）
all_normalized_strides_RAW = all_norm_Global_RAW; % 平面則解析用（絶対角度）

fprintf('Saving results to %s.mat...\n', tTextSaveMatName);

save(tTextSaveMatName, ...
    'data', 'theta', 'centeredtheta', ...
    'averaged_stride', 'std_stride', 'all_normalized_strides', ...
    'all_normalized_strides_RAW', ... 
    'mean_phase_transition_pcts', ... 
    'averaged_grf', 'averaged_grf_R', 'averaged_grf_L', ...
    'std_grf_R', 'std_grf_L', ...
    'all_norm_DS1', 'all_norm_SS1', 'all_norm_DS2', 'all_norm_SS2', ...
    'averaged_DS1', 'averaged_SS1', 'averaged_DS2', 'averaged_SS2', ...
    'valid_strides_R', 'valid_strides_L'); 

fprintf('Processing Complete.\n');