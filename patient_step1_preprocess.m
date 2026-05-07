% =========================================================================
% SCRIPT: preprocess_paretic.m
% DESCRIPTION:
%   脳卒中片麻痺歩行解析用 前処理・解析スクリプト (Step 1)
%   
%   【概要】
%   モーションキャプチャ等で取得したCSVデータを読み込み、以下の前処理を行います。
%     1. データのトリミング（不要な前後フレームのカット）と時間軸のゼロリセット
%     2. 欠損値の補間とローパスフィルタ（LPF）によるノイズ除去
%     3. グローバル座標から矢状面（Sagittal Plane）への座標変換と平均姿勢モデルの構築
%     4. 踵の曲率（Curvature）とつま先位置を用いた歩行イベント（TD: 接地, LO: 離地）の自動検出
%     5. 各種検証用グラフの描画と、次ステップへ引き継ぐためのMATファイルの保存
% =========================================================================
close all; 
clear;
clc;

%% =========================================================================
% --- 1. 解析条件・パラメータ設定 ---
% =========================================================================
% [1-1] 被験者情報と入出力ファイルの設定
% -------------------------------------------------------------------------
SubjectName      = 'KM_ID3';                                     % 被験者ID
tTextFileNameCSV = '2023-06-04-15-17_添木なし　21.28sec.csv';    % 解析対象のCSVファイル名

% 次の解析ステップ（segment等）へ変数を引き継ぐためのMATファイル名
tTextSaveMatName = [SubjectName, '_step1_preprocessed']; 

% [1-2] 前処理パラメータ
% -------------------------------------------------------------------------
frame_cropEnds   = 200;  % 計測開始時・終了時の不安定なデータを除外するためのカットフレーム数
cutoff_freq      = 20;   % ローパスフィルタのカットオフ周波数 [Hz]
lpf_order        = 2;    % バターワースフィルタの次数

% [1-3] 歩行イベント（TD/LO）自動検出用パラメータ
% -------------------------------------------------------------------------
trig_ratio_swing_phase = 0.5; % 遊脚相とみなすための脚角度の閾値（正規化角度の割合）
search_range_TD        = 1.0; % TD（接地）を探索する範囲の係数

% [1-4] グラフ出力・検証用設定
% -------------------------------------------------------------------------
flg_graphSave    = 1;    % グラフを画像として保存するかどうか (1: 保存する, 0: 保存しない)
cols_tiltFwd_csv = [14, 32, 48, 65, 81, 97, 114]; % 前傾角度が格納されているCSVの列番号

% 検証用グラフ（Fig 7～10）のウィンドウ表示サイズと位置 [左端位置, 下端位置, 幅, 高さ]
% ※ 使用するPCのモニター環境に合わせて適宜調整してください。
fig_pos_verification = [100, 100, 800, 250]; 


%% =========================================================================
% --- 2. データのインポートと初期化 ---
% =========================================================================
fprintf('Reading file: "%s" ...\n', tTextFileNameCSV);
try
    % メタデータ（最初の3行）とヘッダー（1行）の計4行をスキップして数値データのみを読み込む
    fullDataMatrix = readmatrix(tTextFileNameCSV, 'NumHeaderLines', 4);
    fprintf(' -> 読み込み完了\n');
catch e
    error('エラー: "%s" が見つからないか、読み込めません。\n詳細: %s', tTextFileNameCSV, e.message);
end

% 列構成の定義とデータの切り出し
col_time  = 1;                          % 1列目: 時間データ
cols_data = 5:size(fullDataMatrix, 2);  % 5列目以降: 各種マーカーや角度などの解析データ（2〜4列目の不要列を除外）

rawData.time = fullDataMatrix(:, col_time);
rawData.data = fullDataMatrix(:, cols_data);
clear fullDataMatrix; % メモリ節約のため一時変数をクリア


%% =========================================================================
% --- 3. データのトリミングと時間軸のゼロリセット ---
% =========================================================================
fprintf('Trimming data (%d frames from both ends)...\n', frame_cropEnds);
[row, ~] = size(rawData.time);

% データ長がカットするフレーム数より短い場合はエラーを出力
if (row < frame_cropEnds * 2)
    error('データ長が短すぎます。frame_cropEnds の値を小さくしてください。');
end

% 指定したフレーム数分、データの前後を削除（ノイズや定常状態でない部分の除外）
rawData.time([1:frame_cropEnds, end-frame_cropEnds+1:end], :) = [];
rawData.data([1:frame_cropEnds, end-frame_cropEnds+1:end], :) = [];

% 最初のフレームの時間を0秒にリセット
rawData.time = rawData.time - rawData.time(1);
fprintf(' -> 時間リセット完了 (Duration: %.2f s)\n', rawData.time(end));


%% =========================================================================
% --- 4. サンプリング周波数の推定とフィルタリング処理 ---
% =========================================================================
% 時間間隔からサンプリング周波数を計算
samplingPeriod_est = mean(diff(rawData.time));
samplingFrequency  = 1 / samplingPeriod_est;
fprintf('Estimated Sampling Freq: %.2f Hz\n', samplingFrequency);

% ローパスフィルタ（バターワース）の設計
fprintf('Applying Low-Pass Filter (Cutoff: %d Hz)...\n', cutoff_freq);
[b, a] = butter(lpf_order, (cutoff_freq * 2) / samplingFrequency);

% フィルタを適用しない列（接地センサーなど、二値的または非連続なデータ）の指定
cols_contact = [61, 110]; 
all_cols_idx = 1:size(rawData.data, 2);
cols_filterable = setdiff(all_cols_idx, cols_contact);

% --- 欠損値 (NaN) の処理 ---
% 1. トラッキング切れなどの一時的な欠損を線形補間で埋める
rawData.data = fillmissing(rawData.data, 'linear');
% 2. 補間しきれない端点の値などは0で埋め、後続処理でのエラーを回避する
rawData.data(isnan(rawData.data)) = 0;

% フィルタ処理の実行
data.time = rawData.time;
data.data = zeros(size(rawData.data));
data.data(:, cols_filterable) = filtfilt(b, a, rawData.data(:, cols_filterable));
data.data(:, cols_contact)    = rawData.data(:, cols_contact); 
fprintf(' -> フィルタ処理完了\n');


%% =========================================================================
% --- 5. [Fig 1] LPF適用後の角度データ確認 ---
% =========================================================================
fprintf('Generating Fig 1: Angle Comparison...\n');

% プロット対象の列インデックス定義
pitch_indices   = [12, 30, 46, 63, 79, 95, 112]; 
tiltFwd_indices = cols_tiltFwd_csv - 1; 
plot_tilt       = ~isempty(cols_tiltFwd_csv) && length(cols_tiltFwd_csv) == 7;

fig1 = figure(1); clf(fig1);
set(fig1, 'Name', 'Fig1_Angle_Comparison', 'NumberTitle', 'off');

for i = 1:7 
    if plot_tilt
        subplot(7, 2, (i*2)-1); 
    else
        subplot(7, 1, i); 
    end
    plot(data.time, data.data(:, pitch_indices(i)), 'b-');
    ylabel('Angle [deg]'); grid on; xlim tight;
    if i == 7, xlabel('Time [s]'); end
    
    if plot_tilt
        subplot(7, 2, i*2); 
        plot(data.time, data.data(:, tiltFwd_indices(i)), 'r-');
        ylabel('Angle [deg]'); grid on; xlim tight;
        if i == 7, xlabel('Time [s]'); end
    end
end


%% =========================================================================
% --- 6. マーカー軌道データの抽出と進行方向（矢状面）の定義 ---
% =========================================================================
fprintf('Extracting Marker Data for Trajectory Analysis...\n');
pelvis_yaw_deg = rawData.data(:, 11);

% 各部位の3次元座標（X, Y, Z）を取得
Pelvis_XYZ = rawData.data(:, 152:154); CoM_XYZ    = rawData.data(:, 137:139);
HipL_XYZ   = rawData.data(:,155:157);  HipR_XYZ   = rawData.data(:,158:160);
KneeL_XYZ  = rawData.data(:,173:175);  KneeR_XYZ  = rawData.data(:,176:178);
AnkleL_XYZ = rawData.data(:,191:193);  AnkleR_XYZ = rawData.data(:,194:196);
ToeL_XYZ   = rawData.data(:,197:199);  ToeR_XYZ   = rawData.data(:,200:202);
MTPL_XYZ   = rawData.data(:,203:205);  MTPR_XYZ   = rawData.data(:,206:208);
HeelL_XYZ  = rawData.data(:,215:217);  HeelR_XYZ  = rawData.data(:,218:220);

% 重心(CoM)の移動軌跡から、歩行の主となる進行方向（Sagittal Plane: 矢状面）を計算
deltaX_CoM = CoM_XYZ(end, 1) - CoM_XYZ(1, 1);
deltaY_CoM = CoM_XYZ(end, 2) - CoM_XYZ(1, 2);
com_sagittal_angle_rad = atan2(deltaY_CoM, deltaX_CoM);
fprintf(' -> Walking Direction (CoM): %.2f [rad]\n', com_sagittal_angle_rad);


%% =========================================================================
% --- 7. グローバル座標系から矢状面座標系への投影変換 ---
% =========================================================================
fprintf('Projecting data to Global Sagittal Plane...\n');

% 進行方向を基準軸(X軸)とするための回転行列成分
theta_rot = com_sagittal_angle_rad; 
c_th = cos(theta_rot); s_th = sin(theta_rot);
rot_X = @(X, Y) X .* c_th + Y .* s_th; % X-Y平面上での回転関数

% 各マーカーのX座標を矢状面へ投影
Pelvis_X_sag = rot_X(Pelvis_XYZ(:,1), Pelvis_XYZ(:,2)); CoM_X_sag = rot_X(CoM_XYZ(:,1), CoM_XYZ(:,2));
HipL_X_sag   = rot_X(HipL_XYZ(:,1), HipL_XYZ(:,2));     HipR_X_sag = rot_X(HipR_XYZ(:,1), HipR_XYZ(:,2));
KneeL_X_sag  = rot_X(KneeL_XYZ(:,1), KneeL_XYZ(:,2));   KneeR_X_sag = rot_X(KneeR_XYZ(:,1), KneeR_XYZ(:,2));
AnkleL_X_sag = rot_X(AnkleL_XYZ(:,1), AnkleL_XYZ(:,2)); AnkleR_X_sag = rot_X(AnkleR_XYZ(:,1), AnkleR_XYZ(:,2));
ToeL_X_sag   = rot_X(ToeL_XYZ(:,1), ToeL_XYZ(:,2));     ToeR_X_sag = rot_X(ToeR_XYZ(:,1), ToeR_XYZ(:,2));
MTPL_X_sag   = rot_X(MTPL_XYZ(:,1), MTPL_XYZ(:,2));     MTPR_X_sag = rot_X(MTPR_XYZ(:,1), MTPR_XYZ(:,2));
HeelL_X_sag  = rot_X(HeelL_XYZ(:,1), HeelL_XYZ(:,2));   HeelR_X_sag = rot_X(HeelR_XYZ(:,1), HeelR_XYZ(:,2));

% Z座標（高さ）はそのまま使用
Pelvis_Z = Pelvis_XYZ(:,3); CoM_Z = CoM_XYZ(:,3);
HipL_Z = HipL_XYZ(:,3); HipR_Z = HipR_XYZ(:,3);
KneeL_Z = KneeL_XYZ(:,3); KneeR_Z = KneeR_XYZ(:,3);
AnkleL_Z = AnkleL_XYZ(:,3); AnkleR_Z = AnkleR_XYZ(:,3);
ToeL_Z = ToeL_XYZ(:,3); ToeR_Z = ToeR_XYZ(:,3);
MTPL_Z = MTPL_XYZ(:,3); MTPR_Z = MTPR_XYZ(:,3);
HeelL_Z = HeelL_XYZ(:,3); HeelR_Z = HeelR_XYZ(:,3);

% 平均歩行速度の算出
distance_sagittal_m = (CoM_X_sag(end) - CoM_X_sag(1)) / 1000;
duration_s          = data.time(end) - data.time(1);
avg_speed_mps       = distance_sagittal_m / duration_s;
dt = diff(data.time); 
inst_velocity_mps = diff(CoM_X_sag / 1000) ./ dt; 
std_speed_mps = std(inst_velocity_mps);
fprintf(' -> Average Walking Speed (Sagittal): %.3f ± %.3f [m/s] (Mean ± SD)\n', avg_speed_mps, std_speed_mps);


%% =========================================================================
% --- 8. [Fig 2] CoMの水平面（XY平面）移動軌跡の描画 ---
% =========================================================================
fprintf('Generating Fig 2: CoM XY Trajectory...\n');
idx_arrows = 1:50:length(rawData.time); 
arrow_len  = 0.3;

fig2 = figure(2); set(fig2, 'Name', 'Fig2_CoM_XY', 'NumberTitle', 'off');
plot(CoM_XYZ(:,1), CoM_XYZ(:,2), 'k-', 'LineWidth', 1); hold on;

% 骨盤の向きを示す矢印を描画
quiver(CoM_XYZ(idx_arrows,1), CoM_XYZ(idx_arrows,2), ...
    cos(com_sagittal_angle_rad + deg2rad(pelvis_yaw_deg(idx_arrows)-mean(pelvis_yaw_deg))), ...
    sin(com_sagittal_angle_rad + deg2rad(pelvis_yaw_deg(idx_arrows)-mean(pelvis_yaw_deg))), ...
    arrow_len, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
xlabel('Global X [mm]'); ylabel('Global Y [mm]'); 
grid on; axis equal; hold off;


%% =========================================================================
% --- 9. [Fig 3] 重心相対での矢状面軌道プロット ---
% =========================================================================
fprintf('Generating Fig 3: Sagittal Trajectory (CoM Relative)...\n');
fig3 = figure(3); set(fig3, 'Name', 'Fig3_Sagittal_Trajectory', 'NumberTitle', 'off');

col_gold = [0.85, 0.6, 0.1]; col_ankle = [0.5, 0, 0.5];

% 描画用インライン関数
plot_traj = @(X, Z, col, name) plot(X - CoM_X_sag, Z, 'Color', col, 'LineStyle', '-', 'DisplayName', name);
plot_pt   = @(X, Z) plot(mean(X - CoM_X_sag), mean(Z), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'HandleVisibility', 'off');
add_txt   = @(X, Z, txt) text(mean(X - CoM_X_sag)+0.02, mean(Z)+0.02, txt, 'FontSize', 8);

% 左下肢軌道
ax1 = subplot(2, 1, 1); hold on;
plot_traj(ToeL_X_sag, ToeL_Z, 'r', 'Toe'); plot_traj(MTPL_X_sag, MTPL_Z, 'g', 'MTP');
plot_traj(HeelL_X_sag, HeelL_Z, 'b', 'Heel'); plot_traj(KneeL_X_sag, KneeL_Z, 'm', 'Knee');
plot_traj(HipL_X_sag, HipL_Z, 'c', 'Hip'); plot_traj(AnkleL_X_sag, AnkleL_Z, col_ankle, 'Ankle');
plot_traj(CoM_X_sag, CoM_Z, 'k', 'CoM'); plot_traj(Pelvis_X_sag, Pelvis_Z, col_gold, 'Pelvis');
plot_pt(ToeL_X_sag, ToeL_Z); add_txt(ToeL_X_sag, ToeL_Z, 'Toe');
plot_pt(HeelL_X_sag, HeelL_Z); add_txt(HeelL_X_sag, HeelL_Z, 'Heel');
plot_pt(AnkleL_X_sag, AnkleL_Z); add_txt(AnkleL_X_sag, AnkleL_Z, 'Ankle');
plot_pt(KneeL_X_sag, KneeL_Z); add_txt(KneeL_X_sag, KneeL_Z, 'Knee');
plot_pt(HipL_X_sag, HipL_Z); add_txt(HipL_X_sag, HipL_Z, 'Hip');
plot_pt(Pelvis_X_sag, Pelvis_Z); add_txt(Pelvis_X_sag, Pelvis_Z, 'Pelvis');
grid on; axis equal; ylabel('Height Z [mm]');
legend('show', 'Location', 'best');

% 右下肢軌道
ax2 = subplot(2, 1, 2); hold on;
plot_traj(ToeR_X_sag, ToeR_Z, 'r', 'Toe'); plot_traj(MTPR_X_sag, MTPR_Z, 'g', 'MTP');
plot_traj(HeelR_X_sag, HeelR_Z, 'b', 'Heel'); plot_traj(KneeR_X_sag, KneeR_Z, 'm', 'Knee');
plot_traj(HipR_X_sag, HipR_Z, 'c', 'Hip'); plot_traj(AnkleR_X_sag, AnkleR_Z, col_ankle, 'Ankle');
plot_traj(CoM_X_sag, CoM_Z, 'k', 'CoM'); plot_traj(Pelvis_X_sag, Pelvis_Z, col_gold, 'Pelvis');
plot_pt(ToeR_X_sag, ToeR_Z); add_txt(ToeR_X_sag, ToeR_Z, 'Toe');
plot_pt(HeelR_X_sag, HeelR_Z); add_txt(HeelR_X_sag, HeelR_Z, 'Heel');
plot_pt(AnkleR_X_sag, AnkleR_Z); add_txt(AnkleR_X_sag, AnkleR_Z, 'Ankle');
plot_pt(KneeR_X_sag, KneeR_Z); add_txt(KneeR_X_sag, KneeR_Z, 'Knee');
plot_pt(HipR_X_sag, HipR_Z); add_txt(HipR_X_sag, HipR_Z, 'Hip');
plot_pt(Pelvis_X_sag, Pelvis_Z); add_txt(Pelvis_X_sag, Pelvis_Z, 'Pelvis');
grid on; axis equal; ylabel('Height Z [mm]'); xlabel('Sagittal X [mm]');
linkaxes([ax1, ax2], 'x');


%% =========================================================================
% --- 10. [Fig 4] 平均姿勢モデル（Stick Figure）の構築と検証 ---
% =========================================================================
fprintf('\n=== Calculating Average Posture Model ===\n');

% [1] 各マーカーの平均基準位置を計算
mL_Hip_X = mean(HipL_X_sag-CoM_X_sag); mL_Hip_Z = mean(HipL_Z);
mR_Hip_X = mean(HipR_X_sag-CoM_X_sag); mR_Hip_Z = mean(HipR_Z);
m_MidHip_X = (mL_Hip_X + mR_Hip_X)/2;  m_MidHip_Z = (mL_Hip_Z + mR_Hip_Z)/2;
m_Pelvis_X = mean(Pelvis_X_sag-CoM_X_sag); m_Pelvis_Z = mean(Pelvis_Z);
mL_Knee_X  = mean(KneeL_X_sag-CoM_X_sag);  mL_Knee_Z  = mean(KneeL_Z);
mL_Ankle_X = mean(AnkleL_X_sag-CoM_X_sag); mL_Ankle_Z = mean(AnkleL_Z);
mL_Toe_X   = mean(ToeL_X_sag-CoM_X_sag);   mL_Toe_Z   = mean(ToeL_Z);
mL_Heel_X  = mean(HeelL_X_sag-CoM_X_sag);  mL_Heel_Z  = mean(HeelL_Z);
mR_Knee_X  = mean(KneeR_X_sag-CoM_X_sag);  mR_Knee_Z  = mean(KneeR_Z);
mR_Ankle_X = mean(AnkleR_X_sag-CoM_X_sag); mR_Ankle_Z = mean(AnkleR_Z);
mR_Toe_X   = mean(ToeR_X_sag-CoM_X_sag);   mR_Toe_Z   = mean(ToeR_Z);
mR_Heel_X  = mean(HeelR_X_sag-CoM_X_sag);  mR_Heel_Z  = mean(HeelR_Z);

% [2] セグメント角度と長さの平均値を算出
rad2deg = @(r) r * 180/pi; deg2rad = @(d) d * pi/180;
calc_mean_ang = @(DistX, DistZ, ProxX, ProxZ) mean(rad2deg(atan2(DistX - ProxX, -(DistZ - ProxZ))));
calc_mean_len = @(X1, Z1, X2, Z2) mean(sqrt((X1 - X2).^2 + (Z1 - Z2).^2));

MidHip_X_sag = (HipL_X_sag + HipR_X_sag) ./ 2; MidHip_Z = (HipL_Z + HipR_Z) ./ 2;

ang_Trunk   = calc_mean_ang(MidHip_X_sag, MidHip_Z, Pelvis_X_sag, Pelvis_Z);
len_Trunk   = calc_mean_len(MidHip_X_sag, MidHip_Z, Pelvis_X_sag, Pelvis_Z);
ang_Thigh_L = calc_mean_ang(KneeL_X_sag, KneeL_Z, HipL_X_sag, HipL_Z); 
len_Thigh_L = calc_mean_len(KneeL_X_sag, KneeL_Z, HipL_X_sag, HipL_Z);
ang_Shank_L = calc_mean_ang(AnkleL_X_sag, AnkleL_Z, KneeL_X_sag, KneeL_Z); 
len_Shank_L = calc_mean_len(AnkleL_X_sag, AnkleL_Z, KneeL_X_sag, KneeL_Z);
ang_Thigh_R = calc_mean_ang(KneeR_X_sag, KneeR_Z, HipR_X_sag, HipR_Z); 
len_Thigh_R = calc_mean_len(KneeR_X_sag, KneeR_Z, HipR_X_sag, HipR_Z);
ang_Shank_R = calc_mean_ang(AnkleR_X_sag, AnkleR_Z, KneeR_X_sag, KneeR_Z); 
len_Shank_R = calc_mean_len(AnkleR_X_sag, AnkleR_Z, KneeR_X_sag, KneeR_Z);

ang_Ank_Heel_L = calc_mean_ang(HeelL_X_sag, HeelL_Z, AnkleL_X_sag, AnkleL_Z); 
len_Ank_Heel_L = calc_mean_len(AnkleL_X_sag, AnkleL_Z, HeelL_X_sag, HeelL_Z);
ang_Ank_Toe_L  = calc_mean_ang(ToeL_X_sag, ToeL_Z, AnkleL_X_sag, AnkleL_Z); 
len_Ank_Toe_L  = calc_mean_len(AnkleL_X_sag, AnkleL_Z, ToeL_X_sag, ToeL_Z);
ang_Ank_Heel_R = calc_mean_ang(HeelR_X_sag, HeelR_Z, AnkleR_X_sag, AnkleR_Z); 
len_Ank_Heel_R = calc_mean_len(AnkleR_X_sag, AnkleR_Z, HeelR_X_sag, HeelR_Z);
ang_Ank_Toe_R  = calc_mean_ang(ToeR_X_sag, ToeR_Z, AnkleR_X_sag, AnkleR_Z); 
len_Ank_Toe_R  = calc_mean_len(AnkleR_X_sag, AnkleR_Z, ToeR_X_sag, ToeR_Z);

calc_mean_ang_foot = @(ToeX, ToeZ, HeelX, HeelZ) mean(rad2deg(atan2(ToeX - HeelX, -(ToeZ - HeelZ)))) - 90;
ang_Foot_L = calc_mean_ang_foot(ToeL_X_sag, ToeL_Z, HeelL_X_sag, HeelL_Z);
ang_Foot_R = calc_mean_ang_foot(ToeR_X_sag, ToeR_Z, HeelR_X_sag, HeelR_Z);

% [3] 算出した平均角度・長さを用いてStick Figureを再構築 (Forward Kinematics)
fk_next = @(prevX, prevZ, len, ang_deg) deal(prevX + len * sin(deg2rad(ang_deg)), prevZ - len * cos(deg2rad(ang_deg)));
[rec_Pelvis_X, rec_Pelvis_Z]   = fk_next(m_MidHip_X, m_MidHip_Z, len_Trunk, ang_Trunk + 180);
[rec_L_Knee_X, rec_L_Knee_Z]   = fk_next(mL_Hip_X, mL_Hip_Z, len_Thigh_L, ang_Thigh_L);
[rec_L_Ankle_X, rec_L_Ankle_Z] = fk_next(rec_L_Knee_X, rec_L_Knee_Z, len_Shank_L, ang_Shank_L);
[rec_L_Heel_X, rec_L_Heel_Z]   = fk_next(rec_L_Ankle_X, rec_L_Ankle_Z, len_Ank_Heel_L, ang_Ank_Heel_L);
[rec_L_Toe_X, rec_L_Toe_Z]     = fk_next(rec_L_Ankle_X, rec_L_Ankle_Z, len_Ank_Toe_L, ang_Ank_Toe_L);
[rec_R_Knee_X, rec_R_Knee_Z]   = fk_next(mR_Hip_X, mR_Hip_Z, len_Thigh_R, ang_Thigh_R);
[rec_R_Ankle_X, rec_R_Ankle_Z] = fk_next(rec_R_Knee_X, rec_R_Knee_Z, len_Shank_R, ang_Shank_R);
[rec_R_Heel_X, rec_R_Heel_Z]   = fk_next(rec_R_Ankle_X, rec_R_Ankle_Z, len_Ank_Heel_R, ang_Ank_Heel_R);
[rec_R_Toe_X, rec_R_Toe_Z]     = fk_next(rec_R_Ankle_X, rec_R_Ankle_Z, len_Ank_Toe_R, ang_Ank_Toe_R);

fprintf('Generating Fig 4: Posture Comparison...\n');
fig4 = figure(4); set(fig4, 'Name', 'Fig4_Avg_Posture_Comp', 'NumberTitle', 'off');
hold on; grid on; axis equal;
xlabel('Sagittal X [mm]'); ylabel('Height Z [mm]');

% 従来の座標平均モデル（破線）
col_old = [0.7 0.7 0.7]; lw_old = 1.5;
plot([mL_Hip_X, mR_Hip_X], [mL_Hip_Z, mR_Hip_Z], '--', 'Color', col_old, 'LineWidth', lw_old);
plot([m_MidHip_X, m_Pelvis_X], [m_MidHip_Z, m_Pelvis_Z], '--', 'Color', col_old, 'LineWidth', lw_old);
plot([mL_Hip_X, mL_Knee_X, mL_Ankle_X, mL_Toe_X], [mL_Hip_Z, mL_Knee_Z, mL_Ankle_Z, mL_Toe_Z], '--', 'Color', col_old, 'LineWidth', lw_old);
plot([mL_Ankle_X, mL_Heel_X], [mL_Ankle_Z, mL_Heel_Z], '--', 'Color', col_old, 'LineWidth', lw_old);
plot([mR_Hip_X, mR_Knee_X, mR_Ankle_X, mR_Toe_X], [mR_Hip_Z, mR_Knee_Z, mR_Ankle_Z, mR_Toe_Z], '--', 'Color', col_old, 'LineWidth', lw_old);
plot([mR_Ankle_X, mR_Heel_X], [mR_Ankle_Z, mR_Heel_Z], '--', 'Color', col_old, 'LineWidth', lw_old);

% 新手法での再構築モデル（実線）
lw_new = 2.0;
plot([mL_Hip_X, m_MidHip_X, mR_Hip_X], [mL_Hip_Z, m_MidHip_Z, mR_Hip_Z], 'k-', 'LineWidth', 2);
plot([m_MidHip_X, rec_Pelvis_X], [m_MidHip_Z, rec_Pelvis_Z], 'k-', 'LineWidth', 3);
plot([mL_Hip_X, rec_L_Knee_X, rec_L_Ankle_X, rec_L_Toe_X], [mL_Hip_Z, rec_L_Knee_Z, rec_L_Ankle_Z, rec_L_Toe_Z], 'b-', 'LineWidth', lw_new);
plot([rec_L_Ankle_X, rec_L_Heel_X], [rec_L_Ankle_Z, rec_L_Heel_Z], 'b-', 'LineWidth', lw_new);
plot([mR_Hip_X, rec_R_Knee_X, rec_R_Ankle_X, rec_R_Toe_X], [mR_Hip_Z, rec_R_Knee_Z, rec_R_Ankle_Z, rec_R_Toe_Z], 'r-', 'LineWidth', lw_new);
plot([rec_R_Ankle_X, rec_R_Heel_X], [rec_R_Ankle_Z, rec_R_Heel_Z], 'r-', 'LineWidth', lw_new);

lgd = legend('Old (Coord Mean)', '','','','','', 'New (Angle Mean)', 'Trunk', 'Left', '', 'Right', '');
set(lgd, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off'); 
hold off;


%% =========================================================================
% --- 11. 仰角データ (theta) の時系列再構築 ---
% =========================================================================
% 生データの変動を平均姿勢に足し合わせることで、各セグメントの角度変動を計算
fprintf('Reconstructing Elevation Angles...\n');
theta = zeros(length(data.time), 7);
get_ts = @(idx) data.data(:, pitch_indices(idx));

theta(:, 1) = ang_Trunk   + (get_ts(1) - mean(get_ts(1)));
theta(:, 2) = ang_Thigh_R + (get_ts(5) - mean(get_ts(5)));
theta(:, 3) = ang_Thigh_L + (get_ts(2) - mean(get_ts(2)));
theta(:, 4) = ang_Shank_R + (get_ts(6) - mean(get_ts(6)));
theta(:, 5) = ang_Shank_L + (get_ts(3) - mean(get_ts(3)));
theta(:, 6) = ang_Foot_R  + (get_ts(7) - mean(get_ts(7)));
theta(:, 7) = ang_Foot_L  + (get_ts(4) - mean(get_ts(4)));

% 後の解析ステップ用に、平均姿勢とそこからの偏差を保存
mean_posture  = mean(theta, 1);
centeredtheta = theta - mean_posture;


%% =========================================================================
% --- 12. 仮想脚 (Hip-Ankle) 仰角の算出 ---
% =========================================================================
% 股関節から足関節を結ぶ仮想的な脚の角度（位相計算等に使用）を計算
fprintf('Calculating Virtual Leg (Hip-Ankle) Angles...\n');

len_Thigh_L = sqrt((mL_Knee_X - mL_Hip_X)^2 + (mL_Knee_Z - mL_Hip_Z)^2);
len_Shank_L = sqrt((mL_Ankle_X - mL_Knee_X)^2 + (mL_Ankle_Z - mL_Knee_Z)^2);
len_Thigh_R = sqrt((mR_Knee_X - mR_Hip_X)^2 + (mR_Knee_Z - mR_Hip_Z)^2);
len_Shank_R = sqrt((mR_Ankle_X - mR_Knee_X)^2 + (mR_Ankle_Z - mR_Knee_Z)^2);

t_TR = deg2rad(theta(:, 2)); t_TL = deg2rad(theta(:, 3));
t_SR = deg2rad(theta(:, 4)); t_SL = deg2rad(theta(:, 5));

vec_HA_L_X = len_Thigh_L * sin(t_TL) + len_Shank_L * sin(t_SL);
vec_HA_L_Z = -(len_Thigh_L * cos(t_TL) + len_Shank_L * cos(t_SL));
vec_HA_R_X = len_Thigh_R * sin(t_TR) + len_Shank_R * sin(t_SR);
vec_HA_R_Z = -(len_Thigh_R * cos(t_TR) + len_Shank_R * cos(t_SR));

theta_HipAnkle_L = rad2deg(atan2(vec_HA_L_X, -vec_HA_L_Z));
theta_HipAnkle_R = rad2deg(atan2(vec_HA_R_X, -vec_HA_R_Z));


%% =========================================================================
% --- 13. 曲率解析および微分処理 (歩行イベント検出用) ---
% =========================================================================
fprintf('Performing Curvature & MTP Z Analysis...\n');
tN_curv = 5; tCutoff_curv = 15; tInterval = 1/samplingFrequency; tGain = 1/500;
data_FrameNum = (1:length(data.time))';
func_normalize_amp = @(v) (v - min(v)) / (max(v) - min(v));

% 曲率計算用のインライン関数
calc_curv = @(rawZ_mm) func_normalize_amp(curvature_numeric([tGain*data_FrameNum, func_normalize_amp(func_LPF(rawZ_mm/1000, tN_curv, tCutoff_curv, samplingFrequency))]));

% [1] TD(接地)イベント検出用: 踵(Heel)の曲率とLPF
heelZ_LPF_L = func_LPF(HeelL_Z/1000, tN_curv, tCutoff_curv, samplingFrequency);
heelZ_LPF_R = func_LPF(HeelR_Z/1000, tN_curv, tCutoff_curv, samplingFrequency);
heel_curv_L = calc_curv(HeelL_Z);
heel_curv_R = calc_curv(HeelR_Z);

% [2] 保存用: MTPのZ座標と曲率、および2階微分（加速度相当）
mtpZ_LPF_L  = func_LPF(MTPL_Z/1000, tN_curv, tCutoff_curv, samplingFrequency);
mtpZ_LPF_R  = func_LPF(MTPR_Z/1000, tN_curv, tCutoff_curv, samplingFrequency);
mtpZ_curv_L = calc_curv(MTPL_Z);
mtpZ_curv_R = calc_curv(MTPR_Z);
mtpZ_DDt_L  = func_diff5pVec(func_diff5pVec(mtpZ_LPF_L, tInterval), tInterval); 
mtpZ_DDt_R  = func_diff5pVec(func_diff5pVec(mtpZ_LPF_R, tInterval), tInterval); 

% [3] LO(離地)イベント検出用: 重心に対するつま先(Toe)の相対X座標
toe_X_rel_L = ToeL_X_sag - CoM_X_sag;
toe_X_rel_R = ToeR_X_sag - CoM_X_sag;


%% =========================================================================
% --- 14. 歩行イベント（TD / LO）の自動検出 ---
% =========================================================================
fprintf('Detecting Gait Events...\n');

% [TD検出] 仮想脚角度のピーク周辺で、踵の曲率が最大となるポイントを探索
[Locs_TD_L, Ranges_TD_L] = detect_TD_Heel(theta_HipAnkle_L, heel_curv_L, func_normalize_amp, trig_ratio_swing_phase, search_range_TD);
[Locs_TD_R, Ranges_TD_R] = detect_TD_Heel(theta_HipAnkle_R, heel_curv_R, func_normalize_amp, trig_ratio_swing_phase, search_range_TD);

% [LO検出] 立脚相の中で、つま先の重心相対位置が最も後方（最小値）になるポイントを探索
[Locs_LO_L, Ranges_LO_L] = detect_LO_SampleDog(theta_HipAnkle_L, toe_X_rel_L, func_normalize_amp, trig_ratio_swing_phase);
[Locs_LO_R, Ranges_LO_R] = detect_LO_SampleDog(theta_HipAnkle_R, toe_X_rel_R, func_normalize_amp, trig_ratio_swing_phase);


%% =========================================================================
% --- 15. [Fig 5 & 6] 脚角度の変動とイベント探索範囲の可視化 ---
% =========================================================================
fprintf('Generating Fig 5 & 6: Leg Angle & Search Ranges (Split)...\n');

% グラフ描画補助関数
plot_events_overlay = @(x_ev, y_data, mk, col) plot(x_ev, y_data(x_ev), mk, 'MarkerFaceColor', col, 'MarkerSize', 8, 'Color', 'k');
func_plot_shaded    = @(ranges, y_min, y_max, col) arrayfun(@(k) patch([ranges(k,1), ranges(k,2), ranges(k,2), ranges(k,1)], [y_min, y_min, y_max, y_max], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off'), 1:size(ranges,1));

% --- Fig 5: Left Leg ---
fig5 = figure(5); set(fig5, 'Name', 'Fig5_Left_LegAngle_SearchRanges', 'NumberTitle', 'off');
hold on;
plot(theta_HipAnkle_L, 'k-', 'LineWidth', 1.2); 
y_range = max(theta_HipAnkle_L) - min(theta_HipAnkle_L);
y_min = min(theta_HipAnkle_L) - 0.1*y_range; 
y_max = max(theta_HipAnkle_L) + 0.1*y_range;
ylim([y_min, y_max]);

% TDとLOの探索範囲をシェードで描画
func_plot_shaded(Ranges_TD_L, y_min, y_max, 'g');
func_plot_shaded(Ranges_LO_L, y_min, y_max, 'y');

% 正規化角度0.5の位置（閾値）を描画
thresh_val_L = min(theta_HipAnkle_L) + 0.5 * (max(theta_HipAnkle_L) - min(theta_HipAnkle_L));
yline(thresh_val_L, '--k', 'LineWidth', 1.5); 

ylabel('\theta_{Leg} [deg]'); grid on; xlim tight; xlabel('Frame');

h_g_L = patch(NaN,NaN,'g','FaceAlpha',0.2,'EdgeColor','none');
h_y_L = patch(NaN,NaN,'y','FaceAlpha',0.2,'EdgeColor','none');
h_line_L = plot(NaN,NaN,'--k','LineWidth',1.5); 
lgd5 = legend([h_g_L, h_y_L, h_line_L], 'TD Search', 'LO Search', 'Threshold (0.5)');
set(lgd5, 'Location', 'eastoutside', 'Box', 'off');
ax5 = gca; 

% --- Fig 6: Right Leg ---
fig6 = figure(6); set(fig6, 'Name', 'Fig6_Right_LegAngle_SearchRanges', 'NumberTitle', 'off');
hold on;
plot(theta_HipAnkle_R, 'k-', 'LineWidth', 1.2); 
y_range = max(theta_HipAnkle_R) - min(theta_HipAnkle_R);
y_min = min(theta_HipAnkle_R) - 0.1*y_range; 
y_max = max(theta_HipAnkle_R) + 0.1*y_range;
ylim([y_min, y_max]);

func_plot_shaded(Ranges_TD_R, y_min, y_max, 'g');
func_plot_shaded(Ranges_LO_R, y_min, y_max, 'y');

thresh_val_R = min(theta_HipAnkle_R) + 0.5 * (max(theta_HipAnkle_R) - min(theta_HipAnkle_R));
yline(thresh_val_R, '--k', 'LineWidth', 1.5);

ylabel('\theta_{Leg} [deg]'); grid on; xlim tight; xlabel('Frame');

h_g_R = patch(NaN,NaN,'g','FaceAlpha',0.2,'EdgeColor','none');
h_y_R = patch(NaN,NaN,'y','FaceAlpha',0.2,'EdgeColor','none');
h_line_R = plot(NaN,NaN,'--k','LineWidth',1.5); 
lgd6 = legend([h_g_R, h_y_R, h_line_R], 'TD Search', 'LO Search', 'Threshold (0.5)');
set(lgd6, 'Location', 'eastoutside', 'Box', 'off');

ax6 = gca;
linkaxes([ax5, ax6], 'x');


%% =========================================================================
% --- 16. [Fig 7 & 8] TDイベントの検出検証 ---
% =========================================================================
fprintf('Generating Fig 7 & 8: TD Verification (Bottom Only)...\n');

% --- Fig 7: Left TD ---
fig7 = figure(7); set(fig7, 'Name', 'Fig7_Left_TD_Verification', 'NumberTitle', 'off');
set(gcf, 'Position', fig_pos_verification); 

yyaxis left;
p_L = plot(heelZ_LPF_L, 'b-'); ylabel('Heel Height Z [m]'); 
y_min_L = min(heelZ_LPF_L)-0.01; y_max_L = max(heelZ_LPF_L)+0.01; ylim([y_min_L, y_max_L]);
hold on;
p_TD = plot(NaN, NaN, '^', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'Color', 'k'); 
if ~isempty(Locs_TD_L), plot_events_overlay(Locs_TD_L, heelZ_LPF_L, '^', 'r'); end

yyaxis right;
p_Curv = plot(heel_curv_L, 'g-', 'LineWidth', 1.2); ylabel('Curvature (Norm)'); ylim([0 1.1]); hold on;
func_plot_shaded(Ranges_TD_L, 0, 1.2, 'g');
p_Shade = patch(NaN,NaN,'g','FaceAlpha',0.2,'EdgeColor','none'); 

grid on; xlim tight; xlabel('Frame');
ax = gca; ax.YAxis(1).Color = 'b'; ax.YAxis(2).Color = [0 0.5 0];

lgd = legend([p_L, p_TD, p_Curv, p_Shade], 'Heel Z', 'TD Event', 'Curvature', 'TD Search Range');
set(lgd, 'Location', 'eastoutside', 'Box', 'off');

% --- Fig 8: Right TD ---
fig8 = figure(8); set(fig8, 'Name', 'Fig8_Right_TD_Verification', 'NumberTitle', 'off');
set(gcf, 'Position', fig_pos_verification); 

yyaxis left;
p_R = plot(heelZ_LPF_R, 'b-'); ylabel('Heel Height Z [m]'); 
y_min_R = min(heelZ_LPF_R)-0.01; y_max_R = max(heelZ_LPF_R)+0.01; ylim([y_min_R, y_max_R]);
hold on;
p_TD = plot(NaN, NaN, '^', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'Color', 'k'); 
if ~isempty(Locs_TD_R), plot_events_overlay(Locs_TD_R, heelZ_LPF_R, '^', 'r'); end

yyaxis right;
p_Curv = plot(heel_curv_R, 'g-', 'LineWidth', 1.2); ylabel('Curvature (Norm)'); ylim([0 1.1]); hold on;
func_plot_shaded(Ranges_TD_R, 0, 1.2, 'g');
p_Shade = patch(NaN,NaN,'g','FaceAlpha',0.2,'EdgeColor','none');

grid on; xlim tight; xlabel('Frame');
ax = gca; ax.YAxis(1).Color = 'b'; ax.YAxis(2).Color = [0 0.5 0];

lgd = legend([p_R, p_TD, p_Curv, p_Shade], 'Heel Z', 'TD Event', 'Curvature', 'TD Search Range');
set(lgd, 'Location', 'eastoutside', 'Box', 'off');


%% =========================================================================
% --- 17. [Fig 9 & 10] LOイベントの検出検証 ---
% =========================================================================
fprintf('Generating Fig 9 & 10: LO Verification (Bottom Only, Yellow Shade)...\n');

% --- Fig 9: Left LO ---
fig9 = figure(9); set(fig9, 'Name', 'Fig9_Left_LO_Verification', 'NumberTitle', 'off');
set(gcf, 'Position', fig_pos_verification); 

p_Toe = plot(toe_X_rel_L, 'b-'); ylabel('Toe X (Relative) [mm]'); hold on;
y_range_x = max(toe_X_rel_L) - min(toe_X_rel_L);
y_min_ax = min(toe_X_rel_L) - 0.1*y_range_x; y_max_ax = max(toe_X_rel_L) + 0.1*y_range_x;
ylim([y_min_ax, y_max_ax]);
func_plot_shaded(Ranges_LO_L, y_min_ax, y_max_ax, 'y');
p_Shade = patch(NaN,NaN,'y','FaceAlpha',0.2,'EdgeColor','none');

p_LO = plot(NaN, NaN, 'v', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'Color', 'k');
if ~isempty(Locs_LO_L), plot_events_overlay(Locs_LO_L, toe_X_rel_L, 'v', 'b'); end

grid on; xlim tight; xlabel('Frame');
lgd = legend([p_Toe, p_LO, p_Shade], 'Toe X', 'LO Event', 'LO Search Range');
set(lgd, 'Location', 'eastoutside', 'Box', 'off');

% --- Fig 10: Right LO ---
fig10 = figure(10); set(fig10, 'Name', 'Fig10_Right_LO_Verification', 'NumberTitle', 'off');
set(gcf, 'Position', fig_pos_verification); 

p_Toe = plot(toe_X_rel_R, 'b-'); ylabel('Toe X (Relative) [mm]'); hold on;
y_range_x = max(toe_X_rel_R) - min(toe_X_rel_R);
y_min_ax = min(toe_X_rel_R) - 0.1*y_range_x; y_max_ax = max(toe_X_rel_R) + 0.1*y_range_x;
ylim([y_min_ax, y_max_ax]);
func_plot_shaded(Ranges_LO_R, y_min_ax, y_max_ax, 'y');
p_Shade = patch(NaN,NaN,'y','FaceAlpha',0.2,'EdgeColor','none');

p_LO = plot(NaN, NaN, 'v', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'Color', 'k');
if ~isempty(Locs_LO_R), plot_events_overlay(Locs_LO_R, toe_X_rel_R, 'v', 'b'); end

grid on; xlim tight; xlabel('Frame');
lgd = legend([p_Toe, p_LO, p_Shade], 'Toe X', 'LO Event', 'LO Search Range');
set(lgd, 'Location', 'eastoutside', 'Box', 'off');


%% =========================================================================
% --- 18. [Fig 11] 歩行パターンの比較（生データ vs 検出結果） ---
% =========================================================================
fprintf('Generating Fig 11: Gait Pattern Comparison...\n');
fig11 = figure(11); set(fig11, 'Name', 'Fig11_Gait_Pattern_Comparison', 'NumberTitle', 'off');

% Subplot 1: センサー生データによる接地判定
ax11_1 = subplot(2, 1, 1); hold on; 
idx_ContactL = 61; idx_ContactR = 110; threshold = 500;    
color_Left = 'b'; color_Right = 'r'; 

binary_L = rawData.data(:, idx_ContactL) > threshold;
idx_start_L = find(diff([0; binary_L; 0]) == 1);
idx_end_L   = find(diff([0; binary_L; 0]) == -1);
if length(idx_start_L) > length(idx_end_L), idx_end_L(end+1) = length(binary_L); end
for k = 1:length(idx_start_L)
    if idx_start_L(k) <= length(rawData.time) && idx_end_L(k) <= length(rawData.time)
        patch([idx_start_L(k), idx_end_L(k), idx_end_L(k), idx_start_L(k)], [1, 1, 2, 2], color_Left, 'EdgeColor', 'none'); 
    end
end

binary_R = rawData.data(:, idx_ContactR) > threshold;
idx_start_R = find(diff([0; binary_R; 0]) == 1);
idx_end_R   = find(diff([0; binary_R; 0]) == -1);
if length(idx_start_R) > length(idx_end_R), idx_end_R(end+1) = length(binary_R); end
for k = 1:length(idx_start_R)
    if idx_start_R(k) <= length(rawData.time) && idx_end_R(k) <= length(rawData.time)
        patch([idx_start_R(k), idx_end_R(k), idx_end_R(k), idx_start_R(k)], [0, 0, 1, 1], color_Right, 'EdgeColor', 'none'); 
    end
end
grid on; xlim([1 length(data.time)]); ylim([-0.5, 2.5]); 
set(gca, 'YTick', [0.5, 1.5], 'YTickLabel', {'Right Foot', 'Left Foot'});
hold off;

% Subplot 2: 今回のアルゴリズムによる検出結果
ax11_2 = subplot(2, 1, 2); hold on;
for j = 1:length(Locs_TD_L)
    t_s = Locs_TD_L(j);
    next_LO = Locs_LO_L(Locs_LO_L > t_s);
    if ~isempty(next_LO), t_e = next_LO(1); else, t_e = length(data.time); end
    patch([t_s, t_e, t_e, t_s], [1, 1, 2, 2], 'b', 'EdgeColor', 'none');
end
for j = 1:length(Locs_TD_R)
    t_s = Locs_TD_R(j);
    next_LO = Locs_LO_R(Locs_LO_R > t_s);
    if ~isempty(next_LO), t_e = next_LO(1); else, t_e = length(data.time); end
    patch([t_s, t_e, t_e, t_s], [0, 0, 1, 1], 'r', 'EdgeColor', 'none');
end
grid on; xlim([1 length(data.time)]); ylim([-0.5, 2.5]);
set(gca, 'YTick', [0.5, 1.5], 'YTickLabel', {'Right Foot', 'Left Foot'});
xlabel('Frame');
hold off;
linkaxes([ax11_1, ax11_2], 'x');


%% =========================================================================
% --- 19. [Fig 12] 検出されたイベント位置の総合的な確認 ---
% =========================================================================
fprintf('Generating Fig 12: Event Location Check...\n');
fig12 = figure(12); set(fig12, 'Name', 'Fig12_Event_Location_Check', 'NumberTitle', 'off');

plot_traj_event = @(X, Z, TD_idx, LO_idx, col, name) ...
    sub_plot_traj_check(X, Z, CoM_X_sag, CoM_Z, TD_idx, LO_idx, col, name);

ax12_1 = subplot(2, 1, 1); hold on;
plot_traj_event(ToeL_X_sag,   ToeL_Z,   Locs_TD_L, Locs_LO_L, 'r', 'Toe');
plot_traj_event(MTPL_X_sag,   MTPL_Z,   Locs_TD_L, Locs_LO_L, 'g', 'MTP');
plot_traj_event(HeelL_X_sag,  HeelL_Z,  Locs_TD_L, Locs_LO_L, 'b', 'Heel');
plot_traj_event(AnkleL_X_sag, AnkleL_Z, Locs_TD_L, Locs_LO_L, [0.5 0 0.5], 'Ankle');
plot_traj_event(KneeL_X_sag,  KneeL_Z,  Locs_TD_L, Locs_LO_L, 'm', 'Knee');
plot_traj_event(HipL_X_sag,   HipL_Z,   Locs_TD_L, Locs_LO_L, 'c', 'Hip');
plot(Pelvis_X_sag - CoM_X_sag, Pelvis_Z, 'Color', [0.85, 0.6, 0.1], 'LineWidth', 0.5);
grid on; axis equal; ylabel('Height Z [mm]'); 
legend('Location', 'eastoutside');

ax12_2 = subplot(2, 1, 2); hold on;
plot_traj_event(ToeR_X_sag,   ToeR_Z,   Locs_TD_R, Locs_LO_R, 'r', 'Toe');
plot_traj_event(MTPR_X_sag,   MTPR_Z,   Locs_TD_R, Locs_LO_R, 'g', 'MTP');
plot_traj_event(HeelR_X_sag,  HeelR_Z,  Locs_TD_R, Locs_LO_R, 'b', 'Heel');
plot_traj_event(AnkleR_X_sag, AnkleR_Z, Locs_TD_R, Locs_LO_R, [0.5 0 0.5], 'Ankle');
plot_traj_event(KneeR_X_sag,  KneeR_Z,  Locs_TD_R, Locs_LO_R, 'm', 'Knee');
plot_traj_event(HipR_X_sag,   HipR_Z,   Locs_TD_R, Locs_LO_R, 'c', 'Hip');
plot(Pelvis_X_sag - CoM_X_sag, Pelvis_Z, 'Color', [0.85, 0.6, 0.1], 'LineWidth', 0.5);
grid on; axis equal; ylabel('Height Z [mm]'); xlabel('Sagittal X (Relative to CoM) [mm]'); 
linkaxes([ax12_1, ax12_2], 'x');


%% =========================================================================
% --- 20. グラフ画像の保存 ---
% =========================================================================
fNames = {
    [tTextSaveMatName, '_1_LPF_Angle_Comparison'],
    [tTextSaveMatName, '_2_CoM_XY_Trajectory'],         
    [tTextSaveMatName, '_3_Global_Sagittal_Trajectory'],
    [tTextSaveMatName, '_4_Average_Stick_Figure'],
    [tTextSaveMatName, '_5_Left_LegAngle_SearchRanges'],   
    [tTextSaveMatName, '_6_Right_LegAngle_SearchRanges'],  
    [tTextSaveMatName, '_7_Left_TD_Verification'],         
    [tTextSaveMatName, '_8_Right_TD_Verification'],        
    [tTextSaveMatName, '_9_Left_LO_Verification'],         
    [tTextSaveMatName, '_10_Right_LO_Verification'],       
    [tTextSaveMatName, '_11_Gait_Pattern_Comparison'],     
    [tTextSaveMatName, '_12_Event_Location_Check']         
};
figs = [fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9, fig10, fig11, fig12];

if flg_graphSave
    fprintf('Saving figures...\n');
    for i = 1:4 
        func_graphSave2(figs(i), fNames{i}, 1, 0, [600, 500]);
    end
    for i = 5:6 
        func_graphSave2(figs(i), fNames{i}, 1, 0, [600, 200]);
    end
    for i = 7:10 
        func_graphSave2(figs(i), fNames{i}, 1, 0, [800, 250]);
    end
    for i = 11:length(figs) 
        func_graphSave2(figs(i), fNames{i}, 1, 0, [600, 400]);
    end
end


%% =========================================================================
% --- 21. 接地フラグの作成とMATファイルへの保存 ---
% =========================================================================
fprintf('Generating Binary Contact Flags & Saving...\n');
contact_flags = zeros(length(theta), 2);
num_frames = length(theta);

% TDからLOの区間を「1(接地)」としてフラグ立て
for k = 1:length(Locs_TD_L)
    t_s = Locs_TD_L(k);
    next_LO = Locs_LO_L(Locs_LO_L > t_s);
    if ~isempty(next_LO), t_e = next_LO(1); else, t_e = num_frames; end
    t_e = min(t_e, num_frames);
    contact_flags(t_s:t_e, 1) = 1;
end
for k = 1:length(Locs_TD_R)
    t_s = Locs_TD_R(k);
    next_LO = Locs_LO_R(Locs_LO_R > t_s);
    if ~isempty(next_LO), t_e = next_LO(1); else, t_e = num_frames; end
    t_e = min(t_e, num_frames);
    contact_flags(t_s:t_e, 2) = 1;
end

% 必要な変数をまとめて次ステップ用に保存
save(tTextSaveMatName, ...
    "SubjectName", ...
    "data", "rawData", "samplingFrequency", ...
    "theta", "centeredtheta", "mean_posture", ... 
    "contact_flags", ... 
    "theta_HipAnkle_L", "theta_HipAnkle_R", ...
    "mtpZ_curv_L", "mtpZ_curv_R", "mtpZ_DDt_L", "mtpZ_DDt_R", ...
    "Locs_TD_L", "Locs_LO_L", "Locs_TD_R", "Locs_LO_R"); 

fprintf('=== ALL PROCESS COMPLETED ===\n');


%% =========================================================================
% --- LOCAL FUNCTIONS (ローカル関数群) ---
% =========================================================================

% -------------------------------------------------------------------------
% 踵の曲率から接地(TD)を検出する関数
% -------------------------------------------------------------------------
function [Locs_TD, Ranges] = detect_TD_Heel(leg_angle_raw, heel_curv, normalize_fn, th_swing, search_range_ratio)
    Locs_TD = []; Ranges = [];
    ang_norm = normalize_fn(leg_angle_raw);
    logic_Swing = (ang_norm >= th_swing);
    d_logic = diff([0; logic_Swing; 0]);
    onsets = find(d_logic == 1); 
    offsets = find(d_logic == -1) - 1;
    
    if ~isempty(offsets) && ~isempty(onsets) && offsets(1) < onsets(1), offsets(1) = []; end
    len = min(length(onsets), length(offsets));
    
    for i = 1:len
        [~, idx_peak_rel] = max(ang_norm(onsets(i):offsets(i)));
        idx_peak = onsets(i) + idx_peak_rel - 1;
        full_len = offsets(i) - idx_peak;
        if full_len < 1, continue; end
        
        search_len = max(1, floor(full_len * search_range_ratio));
        r_start = idx_peak; r_end = idx_peak + search_len;
        Ranges = [Ranges; r_start, r_end];
        
        [~, idx_td_rel] = max(heel_curv(r_start:r_end));
        Locs_TD(end+1) = r_start + idx_td_rel - 1;
    end
end

% -------------------------------------------------------------------------
% つま先の相対位置から離地(LO)を検出する関数
% -------------------------------------------------------------------------
function [Locs_LO, Ranges] = detect_LO_SampleDog(leg_angle_raw, toe_x_rel, normalize_fn, th_swing)
    Locs_LO = []; Ranges = [];
    ang_norm = normalize_fn(leg_angle_raw);
    logic_Swing = (ang_norm >= th_swing);
    d_logic = diff([0; logic_Swing; 0]);
    swing_starts = find(d_logic == 1); 
    swing_ends = find(d_logic == -1) - 1;

    if isempty(swing_ends) && isempty(swing_starts), return; end

    % ケース1: データ冒頭の立脚相区間（最初のSwing開始前）
    if ~isempty(swing_starts)
        if isempty(swing_ends) || swing_starts(1) < swing_ends(1)
            r_start = 1;
            r_end   = swing_starts(1);
            Ranges = [Ranges; r_start, r_end];
            
            segment = toe_x_rel(r_start:r_end);
            [~, idx_min_rel] = min(segment);
            
            % 端点の場合は検出対象外
            if idx_min_rel > 1 && idx_min_rel < length(segment)
                Locs_LO(end+1) = r_start + idx_min_rel - 1;
            end
        end
    end

    % ケース2: 通常の立脚相区間
    for i = 1:length(swing_ends)
        s_idx = swing_ends(i);
        next_starts = swing_starts(swing_starts > s_idx);
        if ~isempty(next_starts), e_idx = next_starts(1); else, e_idx = length(leg_angle_raw); end
        
        if e_idx > s_idx
            Ranges = [Ranges; s_idx, e_idx];
            
            segment = toe_x_rel(s_idx:e_idx);
            [~, idx_min_rel] = min(segment);
            
            % 端点の場合は検出対象外
            if idx_min_rel > 1 && idx_min_rel < length(segment)
                Locs_LO(end+1) = s_idx + idx_min_rel - 1;
            end
        end
    end
end

% -------------------------------------------------------------------------
% 軌道とイベントの検証用プロット関数
% -------------------------------------------------------------------------
function sub_plot_traj_check(X, Z, CoM_X, CoM_Z, TD_idx, LO_idx, col, name)
    X_rel = X - CoM_X;
    plot(X_rel, Z, 'Color', col, 'DisplayName', name);
    if ~isempty(TD_idx)
        plot(X_rel(TD_idx), Z(TD_idx), 'x', 'Color', 'r', 'LineWidth', 1.5, 'MarkerSize', 10, 'HandleVisibility', 'off'); 
    end
    if ~isempty(LO_idx)
        plot(X_rel(LO_idx), Z(LO_idx), 'x', 'Color', 'b', 'LineWidth', 1.5, 'MarkerSize', 10, 'HandleVisibility', 'off'); 
    end
end