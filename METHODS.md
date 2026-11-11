# Analysis methods

## Questions and population

The primary question examines the joint conditional association of Q11 (busy-period coverage) and Q12 (absence coverage) with Q6 (perceived sick-worker exclusion). Secondary questions compare the two Q6 staffing coefficients directly and assess their joint association with Q4 (perceived handwashing). The descriptive exploratory analysis tabulates leadership, absence coverage, and Q6 response combinations. Every observation represents one worker questionnaire. Workers receive equal weight.

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
- Omit each jurisdiction and retain point estimates. Compute the CV3 jackknife covariance for common staffing coefficients around the full estimate, multiplying summed deletion deviations by seven eighths. All eight valid deletions are required. The original sensitivity uses t7 scalar intervals. The post-results extension below adds a joint reference calculation.
- Inspect staffing correlation, joint frequencies, and variance inflation factors conditional on jurisdiction. Do not select predictors from these diagnostics.

The main model remains the designated working summary. Sensitivities do not replace it based on favorable results. The Q4 lowest binary threshold contained eight disagreement responses and failed the separation/information checks; its estimates and intervals are unavailable.

## Missing-data sensitivity

Impute the 19 analysis items Q1–Q4, Q6–Q9, Q11–Q13, Q15, Q17, Q21, and Q23–Q27 using type-1 predictive mean matching with five observed donors. Each incomplete item is predicted from the other items and jurisdiction indicators. Regression-posterior draws use the flat coefficient/inverse-variance prior; donor predictions use fitted coefficients and recipient predictions use posterior draws. Original observed values are unchanged. Recalculate scores after imputation.

Use 50 independent chains of 50 iterations per root seed, 20261111 and 20261112, with distinct L'Ecuyer-CMRG streams. Inspect diagnostics before outcome fitting. Persistent drift or poor mixing triggers extension to 100 and then 200 iterations. For each outcome, include only workers whose outcome was originally observed: Q6 n=577 and Q4 n=576. Imputed outcomes can inform other imputations but do not enter their own outcome regression.

Average staffing log coefficients across completed data, then exponentiate. Also average the Q6 difference. Report no pooled sampling intervals or tests. If Monte Carlo error exceeds 0.01 or the absolute difference between seed runs exceeds 0.02 log-odds units, extend both runs to 100 and then 200 chains. Any unresolved required fit prevents aggregation for that outcome. Both completed 50-chain runs met these tolerances, and all 200 fits passed.

## Post-results extensions

The following procedures were specified after the original results had been inspected. They assess selected assumptions and preserve the original analyses and research questions.

### Joint covariance and jurisdiction contributions

Use the original full-model staffing coefficient vector with both the site-score covariance and full-estimate-centered CV3 covariance. Report each two-parameter Wald statistic and its uncalibrated F(2,7) reference probability. Do not create a new multiplicity family or significance decision. The cluster-jackknife construction follows MacKinnon, Nielsen, and Webb (2025); their binary-logit results do not establish calibration for this eight-jurisdiction ordinal application.

Decompose each covariance into two orthogonal directions and report how the coefficient projections and eigenvalues reconstruct the joint statistic. Decompose each covariance into its eight jurisdiction contributions. Report each jurisdiction's share of total variance and of variance along that covariance's narrowest direction. These are covariance shares, not additive contributions to the Wald statistic, effective degrees of freedom, or deletion criteria. All eight deletions must pass row, outcome, design, and model-identity checks. A failed deletion makes the entire outcome's jackknife extension unavailable.

### Ordinal model and staffing-score restrictions

Fit one partial proportional-odds model per outcome on the original complete cases. Employee Commitment and Resource Availability have separate coefficients at the three cumulative thresholds. Leadership, management, both staffing items, and jurisdiction retain common slopes. The model has 20 parameters. Fit the coherent four-category likelihood and require ordered cumulative probabilities at every observed employee/resource score pair, at the centered reference pair, and under the specified staffing comparisons. Compare the `ordinal::clm` fit with two feasible-start constrained fits. Reject boundary solutions and retain every failed attempt. Partial proportional odds and the crossing problem are described by Agresti (2010, Sections 3.6.1 and 3.6.4).

Separately fit a proportional-odds model with both staffing items represented by their original five categories, using neutral response 3 as reference. This 22-parameter model removes the equal-increment restriction on staffing scores. The numeric staffing model is nested within it. The earlier grouped three-category model is nonnested with the numeric model, but is nested within the five-category model. Report eight staffing-category odds ratios as point estimates. Do not perform an eight-parameter cluster Wald test when the score covariance has rank at most seven.

All new ordinal-model results are point estimates, without new confidence intervals or p-values. For probability summaries, fix each staffing comparison at neutral 3 versus agree 4, retaining each worker's other values. Within each jurisdiction, require at least five complete cases in both relevant joint staffing cells. A display also requires at least 30 eligible workers from at least three jurisdictions. These are reporting rules, not proof of continuous-covariate overlap. Average all four response-category probabilities over the same eligible workers for each model. The direct absence-versus-busy-period comparison concerns Q6 only and uses the intersection of the two eligible populations. Unsupported displays remain unavailable.

### Enriched missing-predictor sensitivity

Jointly complete the original 19 questionnaire items and D1, D2, and D5. Include observed jurisdiction and language. All categorical predictors use fixed reference-cell indicators, retaining documented columns even when a column is zero. This changes auxiliary information, predictor representation, and regularization together. Comparisons with the original numeric-predictor PMM cannot identify the separate contribution of each change. White, Royston, and Wood (2011) discuss auxiliary variables, categorical predictors, diagnostic checks, and departures from missing at random.

Use a proper normal/inverse-gamma working regression for matching. Scale each target to 0 through 1. The coefficient prior has intercept mean 0.5, other means zero, and conditional covariance equal to 100 times the residual variance times the identity matrix. The inverse-gamma variance prior has shape 2 and scale 0.04 under the density proportional to variance raised to minus shape minus one, multiplied by exp(-scale/variance). Its mean is 0.04. These are fixed regularization choices. They permit singular indicator designs and do not define an ordinal response likelihood. Donor matching uses posterior mean predictions for observed donors and a posterior coefficient draw for recipients. Copy an original observed response from the five closest donors. Observed responses remain unchanged.

Run an enriched baseline and four donor-preference scenarios. For originally missing Q11 and Q12 only, multiply each candidate donor's weight by the exponential of delta multiplied by (response minus 3), with delta zero at baseline and all four combinations of minus log(2) and plus log(2) for the two items. The four scenario names are `both_lower`, `both_higher`, `busy_lower` (busy lower and absence higher), and `absence_lower` (busy higher and absence lower). A positive log(2) tilt doubles the higher-versus-lower donor probability ratio for two represented adjacent categories relative to uniform selection. It is not an outcome odds ratio or an estimated nonresponse effect. The tilt cannot introduce a response absent from the local donor set and has no direct effect when all five donors share a category.

Report achieved donor-distribution changes and local donor support. The original file has only three missing Q11 responses and two missing Q12 responses. This is a narrow staffing-response sensitivity; it does not assess all missing culture responses or selection caused by outcome nonresponse. Proper conditional regression posteriors do not prove joint or substantive-model compatibility, as distinguished by Bartlett et al. (2015).

Each scenario starts with 50 independent chains of 50 iterations under each of the original two root seeds. Corresponding initial streams are deliberately reused across scenarios; scenarios are not independent replicates. Review every scenario's diagnostics before fitting the original full proportional-odds model to its originally observed-outcome population. Extend iterations within that scenario to 100 and then 200 for unresolved drift. Extend both roots to 100 and then 200 chains if its coefficient Monte Carlo error exceeds 0.01 or between-seed difference exceeds 0.02. Preserve existing streams and completed chains. A failed required fit prevents the scenario/outcome summary; never average a successful subset. Report mean log coefficients and exponentiate them. Do not pool sampling intervals or tests.

A separate bounded numerical retry applies only to a saved ordinal fit whose sole failure is the primary score tolerance. It preserves the completed dataset, observed-outcome sample, design, and likelihood. It refines the saved primary and alternate estimates and retains the original numerical acceptance limits, with a tighter target score divided by sample size below `1e-8`. The original fits remain unchanged. Final summaries replace only accepted retries and require every planned fit.

## Exploratory displays

Favorable leadership means a mean score at least 4. Absence coverage is inadequate at Q12=1 or 2, neutral at 3, and adequate at 4 or 5. Cross these definitions into six groups. Report group frequency among all 558 complete exploratory records and all four Q6 response percentages within each group. A group below 10 would be flagged without substantive interpretation. No such group occurred. Apply no hypothesis tests or confidence intervals to these descriptive comparisons.

## Numerical acceptance and execution

Require valid probabilities, positive-definite information, a maximum absolute score divided by sample size below 1e-6, and agreement with an independent fit. Investigate coefficient magnitude above 10, information condition number above 1e10, and adjacent threshold gaps below 1e-6. Preserve attempts and warnings.

A higher-precision restart for otherwise stable fits was added after two sensitivity fits narrowly missed the unchanged score tolerance. It uses the same observations, likelihood, variables, and acceptance criteria and retains the original attempt. The independent fit remains required. This numerical correction did not change the main estimates or substitute for the separated binary fit.

The public-package wrapper differs from the analysis-run wrapper only in its default input directory and the methods file included in the provenance manifest. The diagnostic acceptance helper accepts the reviewing person's name. The estimation, preparation, and imputation algorithms are unchanged.

The extension modules preserve the original modules. They bind source and input hashes, preserve failed attempts, and verify mathematical identities and software agreement on synthetic examples before study execution. Synthetic checks validate implementation within their tested scope. They do not validate missingness assumptions, sampling coverage, or causal identification.

## Method sources for the extensions

- Agresti A. *Analysis of Ordinal Categorical Data*. 2nd ed. Wiley; 2010. Sections 3.6.1 and 3.6.4.
- MacKinnon JG, Nielsen MØ, Webb MD. Cluster-robust jackknife and bootstrap inference for logistic regression models. *Econometric Reviews*. Published online July 1, 2025. https://doi.org/10.1080/07474938.2025.2515161
- White IR, Royston P, Wood AM. Multiple imputation using chained equations: Issues and guidance for practice. *Statistics in Medicine*. 2011;30(4):377–399. https://doi.org/10.1002/sim.4067
- Bartlett JW, Seaman SR, White IR, Carpenter JR. Multiple imputation of covariates by fully conditional specification: Accommodating the substantive model. *Statistical Methods in Medical Research*. 2015;24(4):462–487. https://doi.org/10.1177/0962280214521348
