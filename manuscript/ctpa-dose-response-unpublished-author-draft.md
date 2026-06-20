# Pulmonary Artery Size on CT-Pulmonary Angiogram and Risk of Mortality in an Emergency Department Cohort: A Dose-Response Analysis

> Status: Unpublished author manuscript draft converted from `manuscript/source/2023-10-10-dose-response-letter.docx` for repository indexing and review. This is not a peer-reviewed journal article, not an accepted manuscript, and not publisher-formatted text. Word comments and document package metadata were not copied into this Markdown version.

### Authors

Brian W. Locke1,2 MD (ORCID 0000-0002-3588-5238)

Brittany M Scarpato, MD MSCI1

Joseph Bledsoe, MD3,4

Daniel B Knox, MD2

Karen Conner, MD5

Gregory J Stoddard, MStat6

Meghan M Cirulis, MD MSCI1,2,7

C Gregory Elliott, MD2,7

Mark W Dodson, MD PhD2,7

1Division of Pulmonary and Critical Care Medicine, University of Utah, Salt Lake City, UT

2Department of Pulmonary and Critical Care Medicine, Intermountain Medical Center, Murray, UT

3Division of Emergency Medicine, Intermountain Medical Center, Murray, UT

4Department of Emergency Medicine, Stanford Medicine, Stanford, CA

5Division of Radiology, Intermountain Medical Center, Murray, UT

6Division of Epidemiology, University of Utah, Salt Lake City, UT

7Pulmonary Hypertension Care Center, Intermountain Medical Center, Murray, UT

Corresponding author: Brian W Locke (brian.locke@hsc.utah.edu)

Brian W Locke MD. Pulmonary and Critical Care Medicine Fellow. 26 North 1900 East. Wintrobe 701. Salt Lake City, UT. 84132

### Conflicts of Interest

This research was supported by the National Institutes of Health under Ruth L. Kirschstein National Research Service Award 5T32HL105321 from the NIH. (B.W.L.).

B.W.L. also receives funding for an unrelated project from an American Thoracic Society program supported by ResMed, Philips Respironics, and Fisher & Paykel Healthcare.

### Author contributions

Study concept design BWL, BMS, MMC, CGE, MWD

Data acquisition/analysis BMS, JB, DBK, KC, GJS

Manuscript drafting/revision/editing BWL, BMS, CGE, MWD

Literature Review BWL, BMS, CGE, MWD

Statistics BWL, GJS

Guarantor of integrity of the entire study BWL, BMS, MWD

Keywords: Pulmonary Artery Size, Pulmonary Hypertension, Computed Tomography, Mortality Risk

### Word count

Abstract: ***

Manuscript: ~1550

### Abbreviations

***

## Abstract

Purpose: Pulmonary artery (PA) size can be assess by computed tomography pulmonary angiogram (CTPA) as either the PA diameter (PAd) or the ratio with the ascending aorta (PA:AA). When dichotomized, both metrics are associated with mortality risk in the emergency population. However, the excess mortality risk conferred by various degrees of elevation is not known and may differ by age.

Methods: PAd and PA:AA were measured on 889 CTPAs obtained in two community emergency departments in Salt Lake City that were negative for acute pulmonary embolism. Cox-regression modeling with restricted cubic splines was used to model the influence of

Results: Both PAd and PA:AA were associated with larger increases in the hazard of death with greater degrees of enlargement. Enlarged PAd was more common in older patients, but not associated with increased mortality risk until greater degree of enlargement occurred. Elevated PA:AA was more common in younger patients, but excess mortality risk was not conferred until greater enlargement as present. Mortality risk prediction improved if age-adjusted norms were used.

Conclusion: Although verification in an independent dataset is needed, this data suggests that age-adjustment of PA size metrics should be considered .

## Introduction

An enlarged pulmonary artery (PA) on computed tomography of the chest has been associated with adverse outcomes, including mortality, in a variety of situations. PA size is quantified by one of two metrics: the PA diameter (PAd) or a the ratio of the PA diameter to ascending aorta (PA:AA). While differing thresholds have been proposed, the most commonly used thresholds differentiating normal from enlarged are PAd larger than 29 mm in men or 27mm in women and a PA:AA above 0.90. These thresholds are derived from the 90% coverage range in the purportedly healthy Framingham cohort in an analysis by Truong et al.1

There are two major limitations to this approach. First, dichotomizing a continuous measure discards useful information about the degree of elevation2 and potentially causes misleading results near the threshold value.3 One method for addressing this limitation is modeling the influence of the covariate of interest as a continuous measure with a strategy that does not assume a particular functional form (shape of the relationship), such as restricted cubic splines4. Second, the threshold that defines abnormal (uncommon) in the general population is not necessarily the same threshold where the findings are pathologic or have prognostic importance. For example, being very tall may be uncommon but minimally important, while a cholesterol level within the ‘normal range’ may still indicate an increased risk of cardiovascular disease.

We previously reported that enlarged PA by predicted an increased risk of death among patients without pulmonary embolism that received CT pulmonary angiogram (CTPA) in an Emergency department setting5. In that analysis, older patients much more commonly had elevated PAd while younger patients had enlarged PA:AA. Thus, we performed this secondary analysis to evaluate the “dose-response” relationship of enlarged PA size and mortality risk, with additional stratification by age to evaluate whether these frequent elevations were equally correlated with risk of death.

## Methods

This is a secondary analysis of a previously reported cohort study, reported in full elsewhere.5 Briefly, 3160 consecutive CTPAs that were negative for acute pulmonary embolism were obtained in Emergency Departments of Intermountain Medical Center or LDS Hospital (Salt Lake City, 4300 feet elevation) between May 22, 2009 and June 30, 2010. A random sample of 1000 CTPAs were reviewed. Patients were excluded for no follow-up within the health system, complete absence of comorbidity data, or repeat CTPAs.

PAd (at the main PA bifurcation) and PA:AA (AA at the same level as the PA) were manually measured. Demographics, weight, and medical diagnoses rendered prior to the CTPA were extracted from the electronic medical record. Vital status was assessed through October 1, 2021 using health records and local death registries.

The independent association of age (by decade), gender, and medical comorbidities with PA size (both PAd and PA:AA) was modelled using Poisson regression to estimate the prevalence ratio by age using the traditional thresholds from Truong et al1, after correcting for comorbidities. PA diameter was sex-normalized by subtracting 1mm from male values and adding 1mm to female values. Linear regression (mean), and quantile regression at 50th percentile (median) and 90th percentile were used to evaluate the influence of age on the frequency of various degrees of enlargement. Hazard of death was estimated using Cox regression with age (stratified by decade), gender, and PA size represented with a 4-knot restricted cubic splines. This method models the influence of PA size as a continuous variable but allows for the shape of the relationship to be non-linear. Knots were placed at the 5th, 35th, 65th, and 95th percentile.6 Hazard Ratios (HR) are referenced to the healthy population median value (PA:AA 0.77; PAd 25.1mm)1. A post-hoc analysis of the interaction between age and PA-size was performed by generating separate 3-knot restricted cubic splines (10th, 50th, 10th percentiles) for young (<50y), middle (50-65y), and older (65+) patients. Lastly, a model using 4 knot restricted cubic spline of the difference of the PAd and PA:AA from the age-adjusted median was compared to the same modeling using absolute PAd and PA:AA to evaluate the predictive value of age-normalization.

No imputation for missing data was performed. Statistical analyses were performed using Stata version 18 (StataCorp, College Station, TX). The study was approved under a waiver of informed consent by the Intermountain Medical Center Review Board (IRB # ).

## Results

After exclusions (n=10 repeat scans, n=82 with no follow-up, n=23 completely absent comorbidity data), 889 patients were included in the analytic sample. The median age was 52 [interquartile range (IQR) 37, 68] years, 63.9% (n=568) were women, and the median BMI was 29.3 [IQR 24.9,36.2]. The median PAd was 25.3 [IQR 22.5, 28.3] and median PA:AA was 0.85 [IQR 0.76, 0.93]. During a median follow-up of 10.7 years, n=263 (29.6%) patients died (median 3.5 years).

Association of age with enlarged PAd and PA:AA

Older patients were more likely to have an enlarged PAd, while younger patients were more likely to have an elevated PA:AA (table 1, figure 2). In regressions accounting for age and gender, but not comorbidities, the PAd increased with age for the mean (0.56 mm per decade, 95CI 0.42-0.70), median (0.68 mm per decade, 95CI 0.45-0.90), and 90th percentile size (1.0 mm per decade, 95CI 0.58-1.4). Conversely, the PA:AA decreased by . the mean diameter increased by ___ per year, the median diameter increased by ____ per year, and the 90th percentile diameter increased by ____. The ascending aorta generally increased in size faster than the PA with age (Figure 3), so the mean PA:AA decreased by 0.019 per decade (95CI -0.015, -0.023), median ratio decreased by 0.021 per decade (95CI -0.015, -0.027), and 90th percentile ratio decreased by 0.011 per decade (95CI -0.0051, -0.017)

After adjusting for gender, BMI, and comorbidities, patients below age 50 had no different likelihood of enlarged PAd (Relative Risk [RR] 0.82, 95CI 0.61-1.10) but a markedly increased risk of enlarged PA:AA (1.78, 1.36-2.31) compared to patients age 50-65y. Conversely, patients older than 65 had an increased adjusted risk of enlarged PAd (1.58, 1.21-2.05) but no increased risk of enlarged PA:AA (0.97, 0.70-1.34). The independent association with the prevalence of enlarged PA by decadd of age is graphically shown

Association of PA enlargement with risk of death

Cox-regression using restricted cubic splines showed the estimate shape of the relationship between PA size and hazard of death (Figure 1). For gender-normalized PA diameter, a step up in risk occurs around the traditional threshold (28mm, equivalent to 27mm in women and 29mm in men) with little excess risk beyond. Conversely, increases in PA:AA are associated with roughly log-linear increases in hazard of death with increasing abnormality.

Significance of age group

For patients under 50 years old, no excess risk of mortality was seen until PA:AA increased over roughly 1.0 (figure 5), but mortality risk increased at the proposed threshold of PAd. For patients age 50-65, risk began to increase near the usual thresholds for both PA:AA and PAd. For adults over age 65, PAd did not discriminate individuals who are at lower verses higher risk of mortality as well as PA:AA, or as well as PAd in other age strata.

## Discussion

There are three main findings of this secondary analysis of a large cohort of Emergency Department-obtained CTPAs. First, the degree of pulmonary artery enlargement by either PAd or PA:AA leads to a dose-dependent increased risk of mortality. Thus, simply dichotomizing PA size as enlarged or not inefficiently uses information and may give a misleading estimate of the prognostic significance when patients have size values near the threshold.

Second the prevalence of PA enlargement depends on age, with an opposite effect for PAd and PA:AA. PAd enlargement is more common in older adults while an elevated PA:AA is more common in younger adults. The ratio generally decreases with age despite PAd increasing with age because the AA increases even faster. This effect persists after controlling for the effect of accumulating comorbidities with age.

Third, the commonly seen enlargement of PAd in older adults does not seem to be associated with an increased risk of death. Similarly, the commonly seen elevated PA:AA seen in younger adults does not seem to be associated with an excess risk of death until elevation is markedly higher.

These findings suggest that PAd and PA:AA ought to be corrected for age regardless whether the criterion of abnormality is frequency of occurrence or prognostic significance. Roughly 1 in 3 of patients in their 20s or 30s were classified as having elevated PA:AA ratios by the traditional cutoff. However, we did not find evidence that these common, mild elevation in PA:AA are associated with an increased hazard of mortality. It is likely that this age-associated misclassification dilutes the predictive power of both PA size metrics. Thus, previously neutral results of association (e.g. in the general population), might be revisited to assess if associations aresignificant with better classification of whether the PA is enlarged for the patients age, or not

There are a few immediate clinical and research implications of these findings. First, methods such as interval likelihood ratios (for clinical use)3 and modeling strategies that do not assume a linear or dichotomous form (for research use)9 could help avoid inefficient and misleading use of the prognostic information contained in PA size measurements. Second age adjustment, particularly for PA:AA, should be considered.

One principal limitation to the current analysis is that the age stratification is a post-hoc hypotheses that was formulated after seeing the excess frequency of abnormal classifications in each age group. Thus, validation in an independent dataset is needed. However, it should be noted that a similar effect of age was found in the purportedly healthy Framingham cohort in an analysis by Truong et al.1, but this was not integrated into the proposed threshold for classifying abnormality.

Lastly, the mechanisms of the association of PA size with age or mortality are not clear. We evaluated an emergency department population, and thus it would be helpful to replicate this analysis in other circumstances to evaluate if this effect is specific to the emergency department.

In summary, we find that the excess risk conferred by enlarged PAd and elevated PA:AA occurs in a dose-response relationship, but that excess risk varies by age group. In particular, mild elevations in PA:AA are common among patients presenting to the ED, but do not appear to be associated with increased risk of death.

Tables and Figures;

Table 1

*Poisson regression controlling for age category, gender, body mass index, and comorbidities (pulmonary disease, congestive heart disease, pulmonary vascular disease, hypertension, diabetes, peripheral vascular disease, kidney disease). CI = confidence interval). Enlarged PA diameter: 27+mm in women, 29+mm in men. Enlarged PA:AA is 0.9+

Figure 1

The dose response relationship between PA size and hazard of death, as modeled by restricted cubic splines with 4 knots.

Figure 2: Unadjusted PA Size by Age Strata

Figure 3: Size of PA, AA, and Ratio by Age

Figure 4: Adjusted Prevalence of PA enlargement using traditional threshold

Figure 5 Increased Hazard of Mortality by Degree of PA enlargement and age across 3 strata of age

Age strata-referenced hazard of death by size of PA enlargement. For PAd, a 25% (yellow), 50% (orange), and 100% (red) excess risk of death occurs with diameter of at lower size in Age<50 and Age 50-65 as compared to in age > 65. Conversely, a much larger ratio must be seen to confer the same excess risk of death in Age < 50, as compared to 50-65 or 65+.

## References

1. Truong QA, Massaro JM, Rogers IS, et al. Reference values for normal pulmonary artery dimensions by noncontrast cardiac computed tomography: the Framingham Heart Study. Circulation: Cardiovascular Imaging. 2012;5(1):147-154.

2. Altman DG, Royston P. The cost of dichotomising continuous variables. BMJ. 2006;332(7549):1080. doi:10.1136/bmj.332.7549.1080

3. Kohn MA. Key concepts in clinical epidemiology: reporting on the accuracy of continuous tests. Journal of Clinical Epidemiology. 2022;141:157-160. doi:10.1016/j.jclinepi.2021.07.012

4. Desquilbet L, Mariotti F. Dose-response analyses using restricted cubic spline functions in public health research. Statistics in Medicine. 2010;29(9):1037-1057. doi:10.1002/sim.3841

5. Scarpato BM, Locke BW, Bledsoe J, et al. The association between pulmonary artery enlargement and mortality in an Emergency Department population undergoing computed tomography pulmonary angiography. Pulmonary Circulation. 2023;13(2):e12225.

6. Harrell FE. Regression Modeling Strategies: With Applications to Linear Models, Logistic Regression, and Survival Analysis. Vol 608. Springer; 2001.

7. Newman JH. Pulmonary hypertension. American journal of respiratory and critical care medicine. 2005;172(9):1072-1077.

8. Mahapatra S, Nishimura RA, Sorajja P, Cha S, McGoon MD. Relationship of pulmonary arterial capacitance and mortality in idiopathic pulmonary arterial hypertension. Journal of the American College of Cardiology. 2006;47(4):799-803.

9. Collins GS, Ogundimu EO, Cook JA, Manach YL, Altman DG. Quantifying the impact of different approaches for handling continuous predictors on the performance of a prognostic model. Statistics in medicine. 2016;35(23):4124-4135.

## Extracted Tables

> Tables were converted from the DOCX as plain Markdown for discoverability; verify formatting against the source DOCX before reuse.

### Table 1

|  | Total | <50 years | 50-65 years | 65+ years |
| --- | --- | --- | --- | --- |
|  | N=889 | N=425 | N=217 | N=247 |
| PA diameter | PA diameter | PA diameter | PA diameter | PA diameter |
| Median PA diameter (IQR) | 25.3 (22.5,28.3) | 24.3 (21.7,27.0) | 25.8 (23.3,28.7) | 26.7 (23.8,29.9) |
| Enlarged PAd | 28.7% (255) | 19.5% (83) | 29.0% (63) | 44.1% (109) |
| Adjusted Relative Risk* (95% CI) |  | 0.82 (0.61 – 1.10) | 1 (Referent) | 1.58 (1.22-2.05) |
| PA:AA Ratio | PA:AA Ratio | PA:AA Ratio | PA:AA Ratio | PA:AA Ratio |
| Median PA:AA Ratio (IQR) | 0.85 (0.76,0.93) | 0.88 (0.81,0.97) | 0.80 (0.74,0.90) | 0.80 (0.72,0.90) |
| Increased PA:AA | 34.3% (305) | 44.0% (187) | 25.8% (56) | 25.1% (62) |
| Adjusted Relative Risk* (95% CI) |  | 1.78 (1.37 – 2.31) | 1 (Referent) | 0.97 (0.70 – 1.34) |

### Table 2

|  |  |
| --- | --- |

### Table 3

|  |  |
| --- | --- |

### Table 4

|  |
| --- |

### Table 5

|  |  |
| --- | --- |
