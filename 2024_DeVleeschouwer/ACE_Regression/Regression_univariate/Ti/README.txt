Element: Ti

FILES PRODUCED:
- Ti_Model_Ranking.csv : Summary of all regression models for this element,
  including R², RMSEP, RMSE, AIC and BIC, ranked as:
    1) R² (highest → lowest),
    2) RMSEP (lowest → highest),
    3) RMSE (lowest → highest).
- Ti_Diagnostics_Summary.csv :
  Shapiro–Wilk normality test for residuals,
  Breusch–Pagan heteroscedasticity test,
  and 10-fold cross-validation RMSE values for each model.
- *_PredVsObs.pdf : Predicted vs observed plots for each model.
- *_Residuals.pdf : Residual distribution plots for each model.
- *_Influence.pdf : Cook's distance & leverage plots (for LM models only).
- *_PLS_RMSEP.pdf : RMSEP vs number of components for PLS models.

RANKING LOGIC (PER ELEMENT):
1. Rank 1 always corresponds to the model with the highest R².
2. When models share similar R², those with lower RMSEP are ranked higher.
3. When R² and RMSEP are similar, models with lower RMSE are ranked higher.
4. AIC and BIC are reported for context only and are not used in ranking.

DIAGNOSTICS:
- Shapiro–Wilk: tests residual normality; p > 0.05 is usually desirable.
- Breusch–Pagan: tests for heteroscedasticity; p < 0.05 suggests
  non-constant error variance.
- CV10_RMSE: 10-fold cross-validation RMSE; lower values indicate
  better predictive performance.

NOTES:
- AIC/BIC may be NA where they are not defined (e.g., RF, PLS).
- Influence plots are only computed for classical lm() models.
- Site colours in PredVsObs plots use the 'jco' palette with Site
  ordered as: BI10, HER42PB, KER1, KER3, PB1.

Generated automatically by run_full_regressions().
