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
