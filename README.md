#  UWB-Based Vehicle Localization using TDoA & AoA

This MATLAB project implements and compares multiple 2D localization algorithms using **Ultra-Wideband (UWB)** signals — specifically **Time Difference of Arrival (TDoA)** and **Angle of Arrival (AoA)** — to estimate a vehicle’s position in a racetrack environment.

---

## Overview

- **Setup**: 10 UWB access points (APs), 4 vehicle-mounted tags, GPS RTK ground truth
- **Sampling Rate**: 10 Hz
- **Tracks**: 3 (Obstacle, Straight, Noisy)
- **Algorithms**:
  - TDoA-only with Hampel & Kalman filtering
  - AoA-only with EKF & NIS-based tuning
  - TDoA + AoA fusion with fallback to AoA
- **Performance Metrics**: MAE, RMSE, Q95, Availability, Reliability

---

##  Project Files

| Filename                        | Description |
|----------------------------------|-------------|
| `visualizeRawDataStandalone.m`  | 3D AoA ray animation, static AoA visualization, raw TDoA/AoA plots |
| `antennaGT_AoA_processing.m`    | Raw antenna positions, global AoA elevation/azimuth angles |
| `tdoa_only_multitrack.m`        | TDoA-only localization + Hampel + Kalman filter (3 tracks) |
| `aoa_ekf_multitrack_final.m`    | Final AoA-only EKF with NIS-based R tuning and full metrics |
| `tdoa_aoa_fusion_track1.m`      | TDoA+AoA fusion (Track 1), with fallback and NIS |
| `tdoa_aoa_fusion_track2.m`      | TDoA+AoA fusion (Track 2), full metrics and plots |
| `tdoa_aoa_fusion_track3.m`      | TDoA+AoA fusion (Track 3), Kalman filtering, NIS |
| `LNSM_Project_Data.mat`         | Full dataset: raw TDoA/AoA, AP positions, ground truth |
| `localization_Project_presentation 1.pptx` | Final project presentation slides |

---

## How to Run

1. Open MATLAB (R2021a or later recommended).
2. Add all `.m` files and `LNSM_Project_Data.mat` to your working directory.
3. Run any of the analysis scripts directly:

```matlab
run('aoa_ekf_multitrack_final.m')        % AOA-only EKF
run('tdoa_only_multitrack.m')            % TDoA-only all tracks
run('tdoa_aoa_fusion_track2.m')          % Fusion method on Track 2
run('visualizeRawDataStandalone.m')      % Raw AoA rays and animations
