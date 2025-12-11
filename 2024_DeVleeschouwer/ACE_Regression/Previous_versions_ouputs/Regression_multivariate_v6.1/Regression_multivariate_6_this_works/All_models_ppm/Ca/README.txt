Element: Ca

DATA NOTE:
  - Calibration responses (e.g., Ca_ICP) are natural-log transformed (ln mg kg^-1).
  - All model fitting and errors (RMSE, RMSEP) are computed in log-space.
  - Predictions on XRF_new are back-transformed with exp() to mg kg^-1.

PREDICTION COLUMNS IN *_preds.csv:
  - For each model, columns follow the pattern:
      Model_Element_Pred_ppm : predicted concentration in mg kg^-1
      Model_Element_L95_ppm  : lower 95% prediction bound in mg kg^-1
      Model_Element_U95_ppm  : upper 95% prediction bound in mg kg^-1
    e.g. PLS_LOO_Ti_Pred_ppm, PLS_LOO_Ti_L95_ppm, PLS_LOO_Ti_U95_ppm.

MODELS FITTED:
  - OLS, WLS, OLS(wt), WLS(wt), Bayes, RF, PLS(LOO), PLS(k).

RANKING RULES:
  1) Highest R² ranked best.
  2) For similar R², lower RMSEP ranked higher.
  3) For similar R² and RMSEP, lower RMSE ranked higher.
