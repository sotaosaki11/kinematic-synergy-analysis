clear;
close all;

disp("=================================");
disp("ムービーを停止するときには[Ctrl]+C");
disp("=================================");

%% =========================================================================
% 【プログラム概要: スティックピクチャ動画・平均姿勢画像の生成】
% 前処理済みデータ（step1）の仰角（Elevation Angle）を用いて順運動学計算を行い、
% 歩行動作のアニメーション（動画）および1試行の平均姿勢（静止画）を生成・保存する。
% =========================================================================

%% ===================================================================
% 1. 初期設定（入出力ファイル・描画パラメータ）
% ===================================================================
% 描画時の線の太さを一括設定する
set(0, 'DefaultLineLineWidth', 3.0); 

SubjectName = 'SUKEGAWA'; 
InputMatFileName = [SubjectName, '_12_step1_preprocessed.mat'];

targetTrialIdx = 1; % 動画化する試行のインデックスを指定

% --- 動画（Movie）生成の設定 ---
tMovieStep = 4;          % 描画のスキップフレーム数（再生速度の調整用）
flg_kinematicsMovie = 1; % 動画保存フラグ（1:保存する, 0:保存しない）

% --- 平均姿勢（静止画）保存の設定 ---
flg_saveMeanPosture = 1;      % 静止画保存フラグ（1:保存する, 0:保存しない）
cfg_MeanFigSize = [900, 500]; % 保存する画像サイズ [width, height]

%% ===================================================================
% 2. データの読み込みと対象試行の抽出
% ===================================================================
if ~exist(InputMatFileName, 'file')
    error('ファイル %s が見つからない。preprocess.m を先に実行すること。', InputMatFileName);
end

fprintf('Loading data from: %s ...\n', InputMatFileName);
load(InputMatFileName); 

if targetTrialIdx > length(mergedData)
    error('指定されたインデックス(%d)はデータ数(%d)を超えている。', targetTrialIdx, length(mergedData));
end

% 対象試行のデータを抽出する
TargetData = mergedData(targetTrialIdx);
tTheta = TargetData.theta;           % 各体節の時間変化する仰角
tMeanTheta = TargetData.meanposture; % 各体節の時間平均角度（平均姿勢）
TargetTrialName = TargetData.trialName; 

fprintf('Target Trial: %s (Index: %d)\n', TargetTrialName, targetTrialIdx);

outputFileName = ['Movie_', SubjectName, '_', TargetTrialName];

%% ===================================================================
% 3. スティックピクチャのモデル定義読み込み
% ===================================================================
% 各リンク（体節）の長さや接続関係が定義されたスクリプトを読み込む
make_status_7angle; 

%% ===================================================================
% 4. 順運動学計算（動画用：時系列姿勢の算出）
% ===================================================================
% 各体節の仰角から、回転行列を用いて関節のグローバル座標を時系列で算出する。
fprintf('--- Starting Video Processing ---\n');
[NUM_t, ~] = size(tTheta);

% 角度をdegreeからradianに変換し、モデルの基準軸に対する相対角度を求める
tTheta_rad = deg2rad(tTheta - ones(NUM_t, 1) * stickfigure.refAxes');

stickfigure_movie = stickfigure; 
for i = 1:NUM_LINKS
    stickfigure_movie.link{i}.Pglobal = zeros(2, NUM_t);
end

% 各フレーム・各リンクのグローバル座標を計算する
for j = 1:NUM_t
    for i = 1:NUM_LINKS
        if stickfigure_movie.linkDef(i, 1) == 0
            % 親リンクがない（骨盤など起点となるリンク）場合
            stickfigure_movie.link{i}.Pglobal(:, j) = stickfigure_movie.P_0 + func_rot_matrix(tTheta_rad(j, i)) * stickfigure_movie.link{i}.Plocal;
        else
            % 親リンクが存在する場合、親のグローバル座標に自身の相対座標を足す
            parentIdx = stickfigure_movie.linkDef(i, 1);
            stickfigure_movie.link{i}.Pglobal(:, j) = stickfigure_movie.link{parentIdx}.Pglobal(:, j) + func_rot_matrix(tTheta_rad(j, i)) * stickfigure_movie.link{i}.Plocal;
        end
    end
end

%% ===================================================================
% 5. 動画の生成と保存
% ===================================================================
% 描画関数(func_plotKinematics)に渡すためのデータ行列を作成する。
% 形式: [起点X, 起点Z, 0, リンク1X, リンク1Z, 0, ...]
tZeros = zeros(NUM_t, 1);
tData = [ones(NUM_t, 1) * stickfigure_movie.P_0', tZeros]; 

for i = 1:NUM_LINKS
    tData = [tData, stickfigure_movie.link{i}.Pglobal', tZeros];
end

fig1 = figure(1);
set(fig1, 'Color', 'white', 'Name', 'Movie Preview'); 
flg_subplot = 0;

fprintf('Generating animation: %s ...\n', outputFileName);
func_plotKinematics(tData, stickfigure_movie.links, stickfigure_movie.linkColors, flg_subplot, flg_kinematicsMovie, outputFileName, fig1, 0.2, tMovieStep, [0 90;0 90;0 90; 0 90]);
fprintf('Video processing complete.\n\n');

%% ===================================================================
% 6. 順運動学計算および静止画の保存（平均姿勢用）
% ===================================================================
if flg_saveMeanPosture == 1
    fprintf('--- Starting Mean Posture Saving ---\n');
    
    NUM_t_mean = 1; % 平均姿勢は1フレームのみで計算する
    tMeanTheta_rad = deg2rad(tMeanTheta - stickfigure.refAxes');
    
    stickfigure_mean = stickfigure; 
    for i = 1:NUM_LINKS
        stickfigure_mean.link{i}.Pglobal = zeros(2, NUM_t_mean);
    end
    
    % 平均角度に基づくグローバル座標の計算
    j = 1;
    for i = 1:NUM_LINKS
        if stickfigure_mean.linkDef(i, 1) == 0
            stickfigure_mean.link{i}.Pglobal(:, j) = stickfigure_mean.P_0 + func_rot_matrix(tMeanTheta_rad(j, i)) * stickfigure_mean.link{i}.Plocal;
        else
            parentIdx = stickfigure_mean.linkDef(i, 1);
            stickfigure_mean.link{i}.Pglobal(:, j) = stickfigure_mean.link{parentIdx}.Pglobal(:, j) + func_rot_matrix(tMeanTheta_rad(j, i)) * stickfigure_mean.link{i}.Plocal;
        end
    end
    
    % 描画データの結合
    tZeros_mean = zeros(NUM_t_mean, 1);
    tData_mean = [stickfigure_mean.P_0', tZeros_mean]; 
    for i = 1:NUM_LINKS
        tData_mean = [tData_mean, stickfigure_mean.link{i}.Pglobal', tZeros_mean];
    end
    
    % 描画関数内部の min() 関数等で行列の次元（reshape）エラーを回避するため、
    % 便宜的に1フレームのデータを2行（2フレーム分）に複製して渡す。
    tData_mean = [tData_mean; tData_mean]; 
    
    fig2 = figure(2); 
    set(fig2, 'Color', 'white', 'Name', 'Mean Posture');
    
    % 静止画の描画
    func_plotKinematics(tData_mean, stickfigure_mean.links, stickfigure_mean.linkColors, 0, 0, '', fig2, 0.2, 1, [0 90;0 90;0 90; 0 90]);
    
    title(['Mean Posture: ', TargetTrialName], 'Interpreter', 'none', 'FontSize', 14);
    
    outputFileNameMean = ['MeanPosture_', SubjectName, '_', TargetTrialName];
    fprintf('Saving mean posture figure to: %s.* ...\n', outputFileNameMean);
    
    % 画像の保存
    func_graphSave2(fig2, outputFileNameMean, 1, 0, cfg_MeanFigSize);
    
    fprintf('Mean posture saving complete.\n');
end

fprintf('\nAll tasks finished.\n');