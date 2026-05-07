close all; % すべてのグラフウィンドウを閉じる
clear;     % ワークスペース上の変数をすべて消去する

%% =========================================================================
% 【プログラム概要: 前処理 (Preprocessing)】
% モーションキャプチャ（QTM）と床反力計（Bertec）の生データを読み込み、
% 同期処理、ノイズ除去（ローパスフィルタ）、床反力の単位変換、および
% 各体節の仰角（Elevation Angle）の算出を行うプログラム。
% 処理した全試行のデータは最終的に1つの .mat ファイルとして保存される。
% =========================================================================

%% ===================================================================
% 0. 初期設定（解析対象・出力ファイルの設定）
% ===================================================================
SubjectName = 'IBA'; % 被験者名

% 解析対象となる試行（ファイル名の一部）のリスト
TrialList = {
    '12_1_1', ...
    '12_2_1', ...
    '12_3_1'
    };

% 各種パラメータ・定数が定義された外部スクリプトを読み込む
make_status_bertec;   % 床反力計（Bertec）に関する設定（キャリブレーション行列やサンプリング周波数など）
make_status_mocap_29; % モーションキャプチャのマーカー位置や設定に関する定義

tTextFilePath = './';            % データファイルが保存されているディレクトリパス
flg_graphSave = 1;               % グラフ画像を保存するかどうかのフラグ（1:保存する, 0:保存しない）

% 保存される最終的な出力ファイル名（拡張子 .mat は保存時に付与される）
tTextSaveMatName = [SubjectName, '_12_step1_preprocessed'];

% 全試行のデータをまとめて格納するための構造体（Struct）を初期化する
mergedData = struct();

%% ===================================================================
% メインループ: リスト化された全試行ファイルを順番に処理
% ===================================================================
for fIdx = 1:length(TrialList)
    
    TargetTrial = TrialList{fIdx};
    tTextFileName = [SubjectName, '_', TargetTrial]; % 例: 'IBA_12_1_1'
    
    fprintf('\n========== Processing: %s (%d/%d) ==========\n', tTextFileName, fIdx, length(TrialList));

    % 読み込む入力ファイル名の設定（※拡張子や後置詞はロード関数内で処理される想定）
    tTextFileNameMocap  = tTextFileName;             % モーションキャプチャ用ファイル名
    tTextFileNameAnalog = [tTextFileName, '_a'];     % アナログ（床反力）用ファイル名

    %% -------------------------------------------------------------------
    % 1. データの読み込みと同期処理
    % -------------------------------------------------------------------
    fprintf('Step 1: Loading and Syncing Data...\n');

    % カスタム関数を用いてアナログデータとマーカーデータを読み込む
    tDataA = func_loadQTManalog(tTextFilePath, tTextFileNameAnalog);
    rawData.analog = tDataA.data;

    tDataM = func_loadQTMmarker(tTextFilePath, tTextFileNameMocap);
    rawData.mocap = tDataM.data;

    % アナログデータとマーカーデータのフレーム数を揃える（不要な前後フレームのクロップ）
    [row, ~] = size(rawData.analog);
    if row > Bertec.frame_cropEnds * 2
        rawData.analog(row-Bertec.frame_cropEnds:row,:) = [];
        rawData.analog(1:Bertec.frame_cropEnds-1,:) = [];
        rawData.mocap(row-Bertec.frame_cropEnds:row,:) = [];
        rawData.mocap(1:Bertec.frame_cropEnds-1,:) = [];
    end

    % 時間軸データを抽出し、データ本体の行列からは時間列（1〜2列目）を削除する
    rawData.time = rawData.analog(:,2);
    rawData.analog(:,1:2) = []; 
    rawData.mocap(:,1:2)  = []; 

    % 同期確認用のグラフ作成（左右の鉛直方向床反力 Fz をプロット）
    fig1 = figure(1); clf; 
    plot(rawData.time, rawData.analog(:, eNUM_Fz1), 'r'); hold on
    plot(rawData.time, rawData.analog(:, eNUM_Fz2), 'b');
    xlabel('time[s]'); ylabel('Fz[V]'); grid on; xlim tight; ylim padded;
    title(['Sync Check: ', tTextFileName], 'Interpreter', 'none');
    func_font_resize(20); func_line_resize(2);
    func_graphSave(fig1, [tTextFileName, '_sync_check'], flg_graphSave);

    %% -------------------------------------------------------------------
    % 2. アナログ電圧から床反力（N）への単位変換
    % -------------------------------------------------------------------
    fprintf('Step 2: Converting Analog to Force...\n');
    tData.analog = rawData.analog;
    
    % 力(Fx,Fy,Fz)とモーメント(Mx,My,Mz)を格納する変数を初期化する（左右合わせて12成分）
    rawData.force = zeros(size(tData.analog,1), 12);
    
    % Bertecのキャリブレーション行列を用いて、電圧値[V]を物理量[N, N・m]に変換する
    rawData.force(:,eNUM_Fx1:eNUM_Mz1) = tData.analog(:,eNUM_Fx1:eNUM_Mz1) * Bertec.carib_matrix;
    rawData.force(:,eNUM_Fx2:eNUM_Mz2) = tData.analog(:,eNUM_Fx2:eNUM_Mz2) * Bertec.carib_matrix;

    % 変換後の床反力データの確認グラフ
    fig2 = figure(2); clf;
    plot(rawData.time, rawData.force(:,eNUM_Fz1)); hold on;
    plot(rawData.time, rawData.force(:,eNUM_Fz2));
    title(['Force Raw: ', tTextFileName], 'Interpreter', 'none');
    xlabel('time[s]'); ylabel('GRF[N]'); legend('Fz(Right)','Fz(Left)');
    grid on; xlim tight; ylim padded; hold off;
    func_font_resize(20); func_line_resize(2);
    func_graphSave(fig2, [tTextFileName, '_force'], flg_graphSave);

    %% -------------------------------------------------------------------
    % 3. ローパスフィルタ（LPF）によるノイズ除去
    % -------------------------------------------------------------------
    fprintf('Step 3: Applying LPF...\n');
    data.time = rawData.time;

    % フィルタの設計（2次バターワースフィルタ、カットオフ周波数20Hz）
    cutoff = 20; 
    frequency = 300; % サンプリング周波数[Hz]
    order = 2;
    [b,a] = butter(order, (cutoff*2)/frequency); 
    
    % 床反力（Force）にゼロ位相フィルタを適用する（位相遅れを防ぐためfiltfiltを使用）
    data.force = filtfilt(b, a, rawData.force);  

    % LPF適用前後の比較グラフ（床反力）
    plotRange = 10; % プロットする範囲（秒）
    timeWidth = min(length(rawData.time), Bertec.condition.samplingF * plotRange);
    fig3 = figure(3); clf;
    plot(rawData.time(1:timeWidth), rawData.force(1:timeWidth, eNUM_Fz1)); hold on;
    plot(data.time(1:timeWidth), data.force(1:timeWidth, eNUM_Fz1));
    title(['GRF LPF: ', tTextFileName], 'Interpreter', 'none');
    xlabel('time[s]'); ylabel('GRF[N]'); legend('Raw', 'LPF');
    grid on; xlim tight; ylim padded; hold off;
    func_graphSave(fig3, [tTextFileName, '_LPF_GRF'], flg_graphSave);

    % モーションキャプチャのマーカー座標データにも同様のLPFを適用する
    [b,a] = butter(order, (cutoff*2)/frequency); 
    data.mocap = filtfilt(b, a, rawData.mocap);  

    % LPF適用前後の比較グラフ（マーカー：右つま先Z座標の例）
    fig4 = figure(4); clf;
    plot(rawData.time(1:timeWidth), rawData.mocap(1:timeWidth, kinematics.marker.R_Toe*3)); hold on;
    plot(data.time(1:timeWidth), data.mocap(1:timeWidth, kinematics.marker.R_Toe*3));
    title(['Toe Z LPF Check: ', tTextFileName], 'Interpreter', 'none');
    xlabel('time[s]'); ylabel('position[m]'); legend('Raw', 'LPF');
    grid on; xlim tight; ylim padded; hold off;
    func_maximizeFigure(fig4); func_font_resize(20); func_line_resize(2);

    %% -------------------------------------------------------------------
    % 4. 各体節の仰角（Elevation Angle）の算出
    % -------------------------------------------------------------------
    % 鉛直軸を基準とした各リンクの角度を計算する。
    % 以降の主成分分析（PCA）やシナジー解析で重要になる指標。
    fprintf('Step 4: Calculating Elevation Angles...\n');
    mocap_data = data.mocap;

    % 体幹の基準となる中点（左右肩の中点、左右腰の中点）の算出
    Mean_Shoulder_X = (mocap_data(:,10) + mocap_data(:,13)) / 2;
    Mean_Shoulder_Z = (mocap_data(:,12) + mocap_data(:,15)) / 2;
    Mean_Hip_X      = (mocap_data(:,46) + mocap_data(:,49)) / 2;
    Mean_Hip_Z      = (mocap_data(:,48) + mocap_data(:,51)) / 2;
    
    % 各関節マーカーのX, Z座標の抽出
    R_Knee_X = mocap_data(:,52); R_Knee_Z = mocap_data(:,54);
    L_Knee_X = mocap_data(:,55); L_Knee_Z = mocap_data(:,57);
    R_Ankle_X = mocap_data(:,58); R_Ankle_Z = mocap_data(:,60);
    L_Ankle_X = mocap_data(:,61); L_Ankle_Z = mocap_data(:,63);
    R_Toe_X   = mocap_data(:,64); R_Toe_Z   = mocap_data(:,66);
    L_Toe_X   = mocap_data(:,67); L_Toe_Z   = mocap_data(:,69);
    R_Heel_X  = mocap_data(:,70); R_Heel_Z  = mocap_data(:,72);
    L_Heel_X  = mocap_data(:,73); L_Heel_Z  = mocap_data(:,75);

    % 逆正接関数(atan2d)を用いて、鉛直軸に対する各体節の仰角[deg]を算出する
    theta = zeros(size(mocap_data,1), 7);
    theta(:,1) = -atan2d(Mean_Shoulder_X - Mean_Hip_X, Mean_Shoulder_Z - Mean_Hip_Z); % 体幹 (Trunk)
    theta(:,2) = -atan2d(Mean_Hip_X - R_Knee_X, Mean_Hip_Z - R_Knee_Z);               % 右大腿 (Right Thigh)
    theta(:,3) = -atan2d(Mean_Hip_X - L_Knee_X, Mean_Hip_Z - L_Knee_Z);               % 左大腿 (Left Thigh)
    theta(:,4) = -atan2d(R_Knee_X - R_Ankle_X, R_Knee_Z - R_Ankle_Z);                 % 右下腿 (Right Shank)
    theta(:,5) = -atan2d(L_Knee_X - L_Ankle_X, L_Knee_Z - L_Ankle_Z);                 % 左下腿 (Left Shank)
    theta(:,6) = atan2d(R_Toe_Z - R_Heel_Z, R_Toe_X - R_Heel_X);                      % 右足部 (Right Foot)
    theta(:,7) = atan2d(L_Toe_Z - L_Heel_Z, L_Toe_X - L_Heel_X);                      % 左足部 (Left Foot)

    %% -------------------------------------------------------------------
    % 5. 仰角のセンタリング（平均値の減算）
    % -------------------------------------------------------------------
    % 歩行中の姿勢の「変動成分」のみを抽出するために、時間平均を引く。
    fprintf('Step 5: Centering Elevation Angles...\n');
    theta_data = theta;
    meanposture = mean(theta_data, 1);       % 1試行中の各体節の平均角度
    centeredtheta = theta_data - meanposture; % 平均姿勢からの変動分

    %% -------------------------------------------------------------------
    % 6. 処理データの構造体への格納
    % -------------------------------------------------------------------
    % 後続の解析プログラムで扱いやすいよう、試行ごとのデータをまとめる。
    mergedData(fIdx).trialName = TargetTrial;
    mergedData(fIdx).theta = theta;
    mergedData(fIdx).centeredtheta = centeredtheta;
    mergedData(fIdx).meanposture = meanposture;
    mergedData(fIdx).data = data;
    
end

%% ===================================================================
% Final Output: 全試行のデータをまとめて .mat ファイルとして保存
% ===================================================================
fprintf('\nSaving ALL data to: %s.mat\n', tTextSaveMatName);

% 'mergedData' 変数を指定したファイル名で保存する
save(tTextSaveMatName, 'mergedData');

fprintf('All processing complete.\n');