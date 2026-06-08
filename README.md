# 生体力学的歩行解析と運動学シナジー抽出

## 概要

本リポジトリは、修士論文 **「運動学シナジーに基づく片麻痺歩行における脚内・脚間協調構造の解析」** のために開発した MATLAB スクリプト群を収録している。

本プロジェクトは、生のモーションキャプチャデータから高次の協調指標までを扱う一連のデータ処理パイプラインを実装し、ヒトの中枢神経系が「運動学シナジー」を通じて複雑な運動制御をいかに単純化しているかを定量化する。

## 解析ワークフロー

各プログラムはモジュール式の「バケツリレー」方式で設計されており、各ステップがデータを処理して次の段階へ保存する。

### [Step 1] 前処理 (`preprocess.m`)

* **役割:** 生データのクリーニングと座標変換。
* **主な機能:**
    * **ノイズ低減:** 生のマーカーデータに 2 次のバターワース・ローパスフィルタ（20Hz）を適用する。
    * **角度算出:** 7 セグメント（体幹・大腿・下腿・足部）の鉛直軸に対する **仰角 ($\theta$)** を計算する。
    * **イベント検出:** つま先マーカーの **X 座標の極小値** を用いた頑健なアルゴリズムで「離地（Lift-Off）」を、踵の曲率で「接地（Touch-Down）」を検出する。

### [Step 2] 分割と正規化 (`segment.m`)

* **役割:** ストライド分割と時間標準化。
* **主な機能:**
    * **位相分解:** 接地フラグに基づき、各歩行周期を 4 つの位相（DS1、SS1、DS2、SS2）に細分化する。
    * **スプライン補間:** 時系列データを 1 ストライドあたり 200 点（各位相あたり 100 点）に正規化し、解像度を高めて統計的平均化を可能にする。
    * **外れ値除去:** 継続時間の閾値に基づき、不規則なストライド（つまずきなど）を自動的に除外する。

### [Step 3] 平面法則解析 (`planar_law.m`)

* **役割:** 「セグメント間協調の平面法則」の定量化。
* **主な機能:**
    * **幾何モデリング:** 大腿・下腿・足部の仰角が形成する共変動平面を算出する。
    * **対称性指標:** 麻痺側と非麻痺側の法線ベクトルの内積を計算し、協調の非対称性を評価する。

### [Step 4] SVD による シナジー抽出 (`SVD.m`)

* **役割:** 特異値分解を用いた次元削減。
* **主な機能:**
    * **協調モデリング:** **空間モード**（セグメント間の重み付け）と **時間基底**（活性化パターン）を抽出する。
    * **整合性ロジック:** 被験者間・群間で物理的整合性（屈曲・伸展の方向）を保つための手動符号反転ロジックを含む。

### [Step 5] 対称性解析 (`symmetry_analysis.m`)

* **役割:** 歩行回復の臨床的評価。
* **主な機能:**
    * **半周期比較:** 関節インデックスを反転させることで麻痺側データを非麻痺側のシナジー空間へ写像し、直接比較を行う。
    * **軌道可視化:** 2D/3D の投影マップを生成し、患者における協調パターンの変形を可視化する。

### [Visualizer] 動作アニメーション (`make_movie.m`)

* **役割:** 解析結果の物理的検証。
* **主な機能:**
    * 高精細な **スティックフィギュアアニメーション** と平均姿勢のオーバーレイを生成する。
    * 動的な接地インジケータを備え、自動イベント検出の精度を検証する。

## 技術スタック

- **言語:** MATLAB
- **数理手法:** 特異値分解 (SVD)、主成分分析 (PCA)、デジタル信号処理（バターワース LPF）、スプライン補間。
- **応用分野:** ロボティクス、バイオメカニクス、医療 AI、モーションシミュレーション。

---

# Biomechanical Gait Analysis and Kinematic Synergy Extraction

## Overview

This repository contains a suite of MATLAB scripts developed for the graduate thesis: **"Analysis of Intra-limb and Inter-limb Coordination Structures in Hemiparetic Gait Based on Kinematic Synergies"**.

The project implements a full data processing pipeline—from raw motion capture data to high-level coordination metrics—to quantify how the human central nervous system simplifies complex movement control through "Kinematic Synergies."

## Analysis Workflow

The programs are designed in a modular "Bucket Relay" style, where each step processes the data and saves it for the next stage.

### [Step 1] Preprocessing (`preprocess.m`)

* **Role:** Raw data cleaning and coordinate transformation.
* **Key Features:**
    * **Noise Reduction:** Applies a 2nd-order Butterworth low-pass filter (20Hz) to raw marker data.
    * **Angle Calculation:** Computes **Elevation Angles ($\theta$)** for 7 segments (Trunk, Thighs, Shanks, Feet) relative to the vertical axis.
    * **Event Detection:** Employs a robust algorithm using **X-coordinate minima** of toe markers to detect "Lift-Off" and heel curvature for "Touch-Down."

### [Step 2] Segmentation & Normalization (`segment.m`)

* **Role:** Striding and time-standardization.
* **Key Features:**
    * **Phase Decomposition:** Subdivides each gait cycle into four distinct phases (DS1, SS1, DS2, SS2) based on ground contact flags.
    * **Spline Interpolation:** Normalizes time-series data to 200 points per stride (or 100 per phase) to increase resolution and enable statistical averaging.
    * **Outlier Removal:** Automatically excludes irregular strides (e.g., stumbles) based on duration thresholds.

### [Step 3] Planar Law Analysis (`planar_law.m`)

* **Role:** Quantifying the "Planar Law of Intersegmental Coordination."
* **Key Features:**
    * **Geometric Modeling:** Calculates the covariation plane formed by the thigh, shank, and foot elevation angles.
    * **Symmetry Metrics:** Computes the dot product of normal vectors between the paretic and non-paretic limbs to evaluate coordination asymmetry.

### [Step 4] Synergy Extraction via SVD (`SVD.m`)

* **Role:** Dimensionality reduction using Singular Value Decomposition.
* **Key Features:**
    * **Coordination Modeling:** Extracts **Spatial Modes** (inter-segmental weighting) and **Temporal Bases** (activation patterns).
    * **Consistency Logic:** Includes manual sign-flipping logic to ensure physical consistency (flexion/extension directions) across subjects and groups.

### [Step 5] Symmetry Analysis (`symmetry_analysis.m`)

* **Role:** Clinical evaluation of gait recovery.
* **Key Features:**
    * **Half-Cycle Comparison:** Maps paretic limb data onto the non-paretic synergy space by flipping joint indices for direct comparison.
    * **Trajectory Visualization:** Generates 2D/3D projection maps to visualize the deformation of coordination patterns in patients.

### [Visualizer] Motion Animation (`make_movie.m`)

* **Role:** Physical verification of the analysis.
* **Key Features:**
    * Generates high-fidelity **Stick-Figure Animations** and mean posture overlays.
    * Includes dynamic ground contact indicators to verify the accuracy of the automated event detection.

## Technical Stack

- **Language:** MATLAB
- **Mathematical Methods:** Singular Value Decomposition (SVD), Principal Component Analysis (PCA), Digital Signal Processing (Butterworth LPF), Spline Interpolation.
- **Application:** Robotics, Biomechanics, Medical AI, Motion Simulation.
