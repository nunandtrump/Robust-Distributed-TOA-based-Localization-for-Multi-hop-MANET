# Robust-Distributed-TOA-based-Localization-for-Multi-hop-MANET
This repository contains MATLAB implementations of Time of Arrival (TOA) based localization algorithms. These algorithms are designed to estimate the position of a target using noisy measurements from anchor nodes, considering both measurement noise and anchor position uncertainties.

## Files Included

1. **TOA_iCHAN.m**
   - Implements the iCHAN algorithm for target position estimation.
   - Inputs: Noisy anchor positions (`LS_gauss`), noisy distance measurements (`r_measure`), maximum iterations (`Sim_num`), convergence threshold (`step`).
   - Outputs: Estimated target position (`OUT_prev`).

2. **TOA_local_CRLB.m**
   - Computes the Cramér-Rao Lower Bound (CRLB) for TOA-based localization.
   - Inputs: Anchor positions (`LS`), target position (`targ`), TOA measurement noise variance (`gauss_distance`), anchor position uncertainties (`anchor_errors`).
   - Outputs: CRLB value (`CRLB`), inverse of the Fisher Information Matrix (`FIM_inverse`).

3. **TOA_ICPC_static.m**
   - Implements a static algorithm for target position estimation using TOA measurements.
   - Inputs: Noisy anchor positions (`LS_gauss`), noisy distance measurements (`r_measure`), true target position (`targ`), anchor position uncertainties (`anchor_errors`), TOA measurement noise variance (`gauss_distance`).
   - Outputs: Estimated target position (`estimate_position`), Root Mean Square Error (`RMSE`).

4. **TOA_ICPC_dynamic.m**
   - Implements a dynamic algorithm for target position estimation, considering target movement.
   - Inputs: Noisy anchor positions (`LS_gauss`), noisy distance measurements (`r_measure`), true target position (`targ`), previous target position estimate (`targ_initial`), time elapsed (`time`), target velocity (`v_node`), anchor position uncertainties (`anchor_errors`), TOA measurement noise variance (`gauss_distance`).
   - Outputs: Estimated target position (`estimate_position`), Root Mean Square Error (`RMSE`).

5. **calculateCollinearity.m**
   - Determines if anchor nodes are collinear based on their error circles.
   - Inputs: Anchor positions (`anchor_node`), anchor position covariance matrices (`anchor_node_cov`), TOA measurement noise variance (`gauss_distance`).
   - Outputs: Collinearity indicator (`collinearity`), where `0` indicates non-collinear and `1` indicates collinear.
