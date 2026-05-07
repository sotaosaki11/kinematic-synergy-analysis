clear
close all
disp("=================================");
disp("ムービーを停止するときには[Ctrl]+C");
disp("=================================");

%% =================================================================================
% movie設定
tMovieStep = 2;          % movieの早さ． step分だけ描画速度が早くなる． 1～20くらいを指定．
flg_kinematicsMovie = 1; % movieファイルを書き出したいときは1に
outputFileName = 'Gait_Animation_With_Contact'; % 出力ファイル名

% ★ユーザー設定：対象の被験者ID
SubjectName = 'KM_ID3';
%% =================================================================================

%% データ読込
% preprocess_paretic_step1.m で作成した _step1_preprocessed.mat を読み込む
tTextLoadMatName = [SubjectName, '_step1_preprocessed'];

fprintf('Loading data from "%s.mat" ...\n', tTextLoadMatName);
if exist([tTextLoadMatName, '.mat'], 'file')
    load(tTextLoadMatName);
else
    error('ファイルが見つかりません: %s.mat\npreprocess_paretic_step1.m を先に実行してください。', tTextLoadMatName);
end

% 安全対策: thetaがない場合
if ~exist('theta', 'var')
    error('変数 theta が見つかりません。');
end

% 安全対策: contact_flagsがない場合はゼロで作成
if ~exist('contact_flags', 'var')
    warning('contact_flagsが見つかりません。接地表示なしで再生します。');
    contact_flags = zeros(size(theta, 1), 2);
end

tTheta = theta; % 読み込んだデータを関数実行用配列にコピー

%% 定義ファイル読み込み
% ※このファイル(make_status_7angle.m)が同じフォルダにある必要があります
make_status_7angle; 


%% ここから下はいじらなくて良い
%% リンク仰角データに従って，リンク構造の時系列データを生成
[NUM_t, ignore] = size(tTheta);
% ankle以外は垂直軸を0度としているため，水平軸を0度に変換
tTheta = tTheta - ones(NUM_t, 1) * stickfigure.refAxes';
tTheta = deg2rad(tTheta); % deg を radに変換

for i=1:NUM_LINKS
    stickfigure.link{i}.Pglobal = zeros(2, NUM_t);
end

for j=1:NUM_t
    for i=1:NUM_LINKS
        if stickfigure.linkDef(i, 1) == 0
            stickfigure.link{i}.Pglobal(:, j) = stickfigure.P_0 + func_rot_matrix(tTheta(j, i)) * stickfigure.link{i}.Plocal;
        else
            stickfigure.link{i}.Pglobal(:, j) = stickfigure.link{ stickfigure.linkDef(i, 1) }.Pglobal(:, j) + func_rot_matrix(tTheta(j, i)) * stickfigure.link{i}.Plocal;
        end
    end
end

%% 生成したリンク端点の時系列データを元にstickfigureを描画
tZeros = zeros(NUM_t, 1);
tData = ones(NUM_t,1) * stickfigure.P_0';
tData = [tData tZeros];
for i=1:NUM_LINKS
    tData = [tData stickfigure.link{i}.Pglobal' tZeros];
end

fig=figure(1);
flg_subplot = 0;

% 関数呼び出し
% func_plotKinematics_Contact(1:data, 2:links, 3:colors, 4:flg_subplot, 5:flg_CreateMovie, 6:fileName, 7:figure, 8:margin, 9:moviestep, 10:viewangle, 11:contact_flags)

func_plotKinematics_Contact( tData , stickfigure.links, stickfigure.linkColors, flg_subplot, flg_kinematicsMovie, outputFileName, fig, 0.2, tMovieStep, [0 90;0 90;0 90; 0 90], contact_flags);

%% =========================================================================
% --- LOCAL FUNCTIONS ---
% =========================================================================

function func_stickplot3_Contact(Link, Node, Color, inFig, inPlotCommandText, inContact)
% func_stickplot3_Contact (修正版)
% 
% inContact: [1x2] フラグ (1:Left Contact, 2:Right Contact)
%

if nargin < 6
    inContact = [0, 0];
end
if nargin < 5
    inPlotCommandText = '';
end
   
[M, ~] = size(Link);

figure(inFig);

% 色の定義
col_Red  = [1, 0, 0];  % 赤 (左足用: Left)
col_Blue = [0, 0, 1];  % 青 (右足用: Right)
width_Normal  = 0.5;   % 通常時の太さ
width_Contact = 4.0;   % 接地時の太さ

for i=1:M
    nodeA = Link(i, 1);
    nodeB = Link(i, 2);
    X = [ Node(nodeA,:); Node(nodeB,:)];
    
    % デフォルトの色と太さ (make_status_7angleで定義された色を使用)
    tColor = Color(i, :);
    tLineWidth = width_Normal;
    
    % --- 接地判定ロジック (修正済み) ---
    
    % Link 6 = Right Foot (thetaの6列目)
    % 右足接地フラグは inContact(2)
    if i == 6 && inContact(2) == 1
        tColor = col_Red;        
        tLineWidth = width_Contact;
    end
    
    % Link 7 = Left Foot (thetaの7列目)
    % 左足接地フラグは inContact(1)
    if i == 7 && inContact(1) == 1
        tColor = col_Blue;         
        tLineWidth = width_Contact;
    end
    % ---------------------------------------------

    plot3( X(:,1), X(:,2), X(:,3), 'Color', tColor, 'LineWidth', tLineWidth);
    hold on;
end

% --- ラベル表示 ---
str_L = 'Right : Red';
str_R = 'Left: Blue';

% テキスト表示位置の調整 (左上)
text(0.05, 0.95, str_L, 'Units', 'normalized', 'FontSize', 12, ...
     'Color', 'r', 'FontWeight', 'bold');
text(0.05, 0.90, str_R, 'Units', 'normalized', 'FontSize', 12, ...
     'Color', 'b', 'FontWeight', 'bold');

eval(inPlotCommandText);
hold off;

end

function func_plotKinematics_Contact(inData, inLinks, inColors, flg_subplot, flg_CreateMovie, inFileName, inFigure, inMargin, inStep, inViewagl, inContactFlags)
% func_plotKinematics_Contact
% 接地フラグ(contact_flags)を受け取り、フレームごとに処理を回す関数

if nargin < 11
    error('エラー: 接地フラグ(contact_flags)を第11引数に渡してください。');
end
if nargin<10
    inViewagl = [-24 18; 0 90; 90 0; 0 0;];
end
if nargin<9
    inStep = 4;
end
if nargin<8
    inMargin = 0.2; % 0.2m の表示マージン
end
if nargin<7
    inFigure=figure(1);
end
if nargin<6
    inFileName = 'kinematics_plot';
end

%% プロットレンジ決定================
tMargin = inMargin;
tMin = min(reshape(min(inData)', 3, [])') - tMargin;
tMax = max(reshape(max(inData)', 3, [])') + tMargin;
plotRange = reshape( [tMin;tMax], 1, 6);

%% movie設定 ==========================
outputFileName = inFileName;
fig=figure(inFigure);
set(fig,'DoubleBuffer','on');

% movie
if flg_CreateMovie
    v = VideoWriter([outputFileName '.mp4'], 'MPEG-4');
	v.Quality = 100;
	v.FrameRate = 30;
    open(v);
end

%% ===============================
[row, ~] = size(inData);
tStep = inStep;
tEnd = row;
tCounter = 1; 
viewagl = inViewagl;

for k = 1:tStep:tEnd
    % 現在のフレームの接地情報を取得 [Left, Right]
    if k <= size(inContactFlags, 1)
        currentContact = inContactFlags(k, :);
    else
        currentContact = [0, 0];
    end

    %===============================
	% Plot関数実行
    if flg_subplot
        for i=1:4
            tH(i) = subplot(2,2,i);
            func_stickplot3_Contact(inLinks, reshape(inData(k, :), 3, [])', inColors, fig, '', currentContact);
            view(viewagl(i,:));
            set(gca, 'XLim', plotRange(1:2), 'YLim', plotRange(3:4), 'ZLim', plotRange(5:6), 'DataAspectRatio', [1 1 1]);
        end
    else
        func_stickplot3_Contact(inLinks, reshape(inData(k, :), 3, [])', inColors, fig, '', currentContact);
        axis(plotRange);
        view(viewagl(1,:));
        set(gca, 'XLim', plotRange(1:2), 'YLim', plotRange(3:4), 'ZLim', plotRange(5:6), 'DataAspectRatio', [1 1 1]);
    end
    
    if flg_CreateMovie == 1
		Fframe = getframe(fig);
        writeVideo(v,Fframe);
    else
        drawnow;
    end
    tCounter = tCounter+1;
end

%% フィギュア===============================
if flg_CreateMovie == 1
	close(v);
    fprintf('動画保存完了: %s.mp4\n', outputFileName);
end
end