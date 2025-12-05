Element: Fe

FILES PRODUCED:
- Fe_Model_Ranking.csv : Ranked list of all regression models for this element.
- Fe_Diagnostics_Summary.csv : Normality, heteroscedasticity, cross-validation diagnostics per model.
- *_PredVsObs.pdf : Predicted vs observed plots by model.
- *_Residuals.pdf : Residual distribution plots by model.
- *_Influence.pdf : Cook's distance & leverage plots (LM models only).
- *_PLS_RMSEP.pdf : PLS RMSEP vs components (PLS models only).

RANKING LOGIC:
1. All PLS models (PLS_LOO, PLS_kfold) are ranked above OLS and WLS models,
   because PLS is better suited to multicollinearity and can cope better
   with heteroscedastic behaviour in geochemical data.
2. Among non-PLS models, if R2, RMSE and RMSEP are equal,
   WLS models are preferred over OLS_weighted, which are preferred over OLS,
   because they explicitly account for measurement error/variance.
3. Within each priority group, models are ranked by:
   - Highest R2 (better fit),
   - Lowest RMSE (better calibration error),
   - Lowest RMSEP (better cross-validated prediction error),
   - Lowest AIC (if available),
   - Lowest BIC (if available).

DIAGNOSTICS EXPLANATIONS:
- Shapiro–Wilk: Tests residual normality; p > 0.05 suggests residuals are approximately normal.
- Breusch–Pagan: Tests heteroscedasticity; p < 0.05 suggests non-constant residual variance.
- CV10 RMSE: 10-fold cross-validation error; lower values imply better predictive power.
- AIC/BIC: Penalised likelihood criteria; lower values indicate better trade-off between fit and complexity.

NOTES:
- Missing AIC/BIC (NA) means those criteria are not defined for that model type (e.g. random forest, PLS).
- Residual-based diagnostics and influence plots are only defined for LM-class models.
- All Site colours in PredVsObs plots follow ggpubr/ggsci 'jco' palette with Site-level ordering BI10, HER42PB, KER1, KER3, PB1.

Generated automatically by run_full_regressions().
