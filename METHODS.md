# Analysis methods

## Questions and population

The primary question tests the joint conditional association of Q11 (busy-period coverage) and Q12 (absence coverage) with Q6 (perceived sick-worker exclusion). Secondary questions compare the two Q6 staffing coefficients directly and assess their joint association with Q4 (perceived handwashing). The descriptive exploratory analysis tabulates leadership, absence coverage, and Q6 response combinations. Every observation represents one worker questionnaire. Workers receive equal weight.

The public file contains 579 workers from eight jurisdictions, surveyed in March 2018 through March 2019. No restaurant identifier is available. The analysis treats jurisdictions as independent and allows dependence among workers within a jurisdiction.

## Coding

Map the exact labels strongly disagree, disagree, neither agree nor disagree, agree, and strongly agree to 1 through 5. Blanks remain missing. Unknown labels stop processing.

The main outcomes combine responses 1 and 2, retaining 3, 4, and 5 separately. Staffing items remain numeric 1–5. Scores are arithmetic means requiring every constituent item for complete-case analysis:

| Score | Items |
| --- | --- |
| Leadership | Q21, Q23, Q24, Q25, Q26, Q27 |
| Employee Commitment for Q6 | Q1, Q2, Q3, Q4 |
| Employee Commitment for Q4 | Q1, Q2, Q3 |
| Resource Availability | Q7, Q8, Q9 |
| Management Commitment as stored | Q13, Q15, Q17 |

Numeric predictors are centered at 3. The management score uses stored values; prior reverse coding remains unresolved. D8 is excluded because its coding is unresolved. Jurisdiction 1 is the model reference.

## Main models and uncertainty

Cumulative proportional-odds logistic regression includes both staffing items, four applicable scores, and seven jurisdiction indicators. Higher outcome categories mean stronger agreement with the perceived practice. Full and staffing-omitted models use identical workers. Q6 uses 539 complete cases; Q4 uses 541.

The covariance uses individual scores for all model parameters, summed within jurisdiction, and the full observed information matrix. The only covariance multiplier is eight divided by seven. Scalar intervals and the direct coefficient difference use a t reference with seven degrees of freedom. Joint two-parameter comparisons use the Wald statistic divided by two with an F reference on two and seven degrees of freedom. These are approximate reference calculations with eight clusters. There is no binary significance decision rule.

Holm adjustment covers exactly the two secondary reference p-values. Confidence intervals remain pointwise. The Q6 staffing difference exponentiates to a ratio of odds ratios. Log-likelihood gains describe apparent fit on identical observations; no likelihood-ratio p-value, prediction validation, or causal estimate is reported.

## Diagnostics and sensitivities

- Fit separate binary regressions at every main outcome threshold, allowing every substantive slope to vary. Compare coefficients and predicted probabilities with the ordinal fit.
- Treat staffing as three categories, with neutral as reference; retain all five Q6 response levels; separately combine Q4 disagreement and neutrality into a three-level outcome.
- Omit Q13 from the management score and reverse Q13, Q15, or Q17 separately. Global reversal is an algebraic invariance check.
- Add food-service tenure (D1), current-restaurant tenure (D5), training at the current restaurant (D2), and survey language. Use five original tenure levels, with less than one year as reference; use no training and English as other references. Refit the unaugmented model on the same restricted sample.
- Omit each jurisdiction and retain point estimates. Compute the CV3 jackknife covariance for common staffing coefficients around the full estimate, multiplying summed deletion deviations by seven eighths. All eight valid deletions are required. Use t7 scalar intervals, without additional tests.
- Inspect staffing correlation, joint frequencies, and variance inflation factors conditional on jurisdiction. Do not select predictors from these diagnostics.

The main model remains the designated working summary. Sensitivities do not replace it based on favorable results. The Q4 lowest binary threshold contained eight disagreement responses and failed the separation/information checks; its estimates and intervals are unavailable.

## Missing-data sensitivity

Impute the 19 analysis items Q1–Q4, Q6–Q9, Q11–Q13, Q15, Q17, Q21, and Q23–Q27 using type-1 predictive mean matching with five observed donors. Each incomplete item is predicted from the other items and jurisdiction indicators. Regression-posterior draws use the flat coefficient/inverse-variance prior; donor predictions use fitted coefficients and recipient predictions use posterior draws. Original observed values are unchanged. Recalculate scores after imputation.

Use 50 independent chains of 50 iterations per root seed, 20261111 and 20261112, with distinct L'Ecuyer-CMRG streams. Inspect diagnostics before outcome fitting. Persistent drift or poor mixing triggers extension to 100 and then 200 iterations. For each outcome, include only workers whose outcome was originally observed: Q6 n=577 and Q4 n=576. Imputed outcomes can inform other imputations but do not enter their own outcome regression.

Average staffing log coefficients across completed data, then exponentiate. Also average the Q6 difference. Report no pooled sampling intervals or tests. If Monte Carlo error exceeds 0.01 or the absolute difference between seed runs exceeds 0.02 log-odds units, extend both runs to 100 and then 200 chains. Any unresolved required fit prevents aggregation for that outcome. Both completed 50-chain runs met these tolerances, and all 200 fits passed.

## Exploratory displays

Favorable leadership means a mean score at least 4. Absence coverage is inadequate at Q12=1 or 2, neutral at 3, and adequate at 4 or 5. Cross these definitions into six groups. Report group frequency among all 558 complete exploratory records and all four Q6 response percentages within each group. A group below 10 would be flagged without substantive interpretation. No such group occurred. Apply no hypothesis tests or confidence intervals to these descriptive comparisons.

## Numerical acceptance and execution

Require valid probabilities, positive-definite information, a maximum absolute score divided by sample size below 1e-6, and agreement with an independent fit. Investigate coefficient magnitude above 10, information condition number above 1e10, and adjacent threshold gaps below 1e-6. Preserve attempts and warnings.

A higher-precision restart for otherwise stable fits was added after two sensitivity fits narrowly missed the unchanged score tolerance. It uses the same observations, likelihood, variables, and acceptance criteria and retains the original attempt. The independent fit remains required. This numerical correction did not change the main estimates or substitute for the separated binary fit.

The public-package wrapper differs from the analysis-run wrapper only in its default input directory and the methods file included in the provenance manifest. The diagnostic acceptance helper accepts the reviewing person's name. The estimation, preparation, and imputation algorithms are unchanged.
