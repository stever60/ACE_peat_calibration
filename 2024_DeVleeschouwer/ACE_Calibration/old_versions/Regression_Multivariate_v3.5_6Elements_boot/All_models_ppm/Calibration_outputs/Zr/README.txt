Element: Zr

Calibration responses are ln(mg kg^-1).
Predictions stored as *_Pred_ppm, *_L95_ppm, *_U95_ppm.
Additional uncertainty columns:
 *_L95_log_pc, *_U95_log_pc : % error in log-space.
 *_L95_pc_ppm, *_U95_pc_ppm : % error bounds in ppm.

Models: OLS, WLS, OLS(wt), WLS(wt), Bayes, RF, PLS(LOO), PLS(k).
Ranking: highest R², then lowest RMSEP, then lowest RMSE.
Plots show predicted concentrations (Y-axis) versus observed ICP-MS values (X-axis); horizontal error bars represent ICP analytical uncertainty, vertical error bars represent model prediction uncertainty (95% CI).
Model robustness was evaluated using a composite signal-to-noise framework incorporating prediction uncertainty, downcore smoothness, and calibration performance. Signal-to-noise ratios were calculated from prediction confidence intervals and profile roughness, scaled within each element, and combined with R² into a weighted robustness score. Models exhibiting unstable downcore behaviour were excluded from selection. For each element, the highest-ranked stable model was designated as the production model and used for subsequent interpretation.
        
