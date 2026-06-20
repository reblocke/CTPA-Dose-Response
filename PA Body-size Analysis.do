version 17.0
args input_root output_root
if "`input_root'" == "" local input_root "data/private"
if "`output_root'" == "" local output_root "outputs/stata"

capture mkdir "outputs"
capture mkdir "`output_root'"
local output_date = subinstr(strtrim(c(current_date)), " ", "-", .)
local output_dir "`output_root'/`output_date'"
capture mkdir "`output_dir'"
capture mkdir "`output_dir'/Logs"
local a1=substr(c(current_time),1,2)
local a2=substr(c(current_time),4,2)
local a3=substr(c(current_time),7,2)
local b = "PA Body-size Analysis.do" // do file name
capture copy "`b'" "`output_dir'/Logs/(`a1'_`a2'_`a3')`b'", replace
capture log close _all
log using "`output_dir'/Logs/(`a1'_`a2'_`a3')PA Body-size Analysis.log", replace text

clear
set scheme cleanplots //cleanplots white_tableau white_w3d //3 options I usually use

local cleaned_data "`input_root'/cleaned_noempi.dta"
capture confirm file "`cleaned_data'"
if _rc {
	di as error "Required restricted input not found: `cleaned_data'"
	exit 601
}
use "`cleaned_data'", clear

/* Analysis */ 

/* 
Clinical question / point:

Framing: when evaluating PA diameter (or related metrics) w.r.t. their prognostic
implications, body size is a nuissance variable (something that confounds the desired
exposure->outcome relationship).

However, using ratios such as the PA to AA ratio, as a prognostic factor is 
potentially misleading, because there are multiple ways to achieve the same 
ratio which each may have different prognostic implications - for example, 
a normal PA:AA ratio might result from a small PA and small AA, a normal PA and normal AA, 
or a large PA and a large AA - each of which may have different implciations. 

We previously reported that both PAd and PA:AA were associated with mortality risk 
in a cohort of patients who'd receive PA imaging in an emergency room setting. 
In this investigation, we seek re-analyze this data to evaluate: 



1. Does resolving PA:AA collapse potentially useful prognostic information from the 
PA and AA?

2. Does correcting for AA size truly correct for body size, or does it capture 
other information about comorbidities?  Does correcting for the body size nuisance variable by correct for
 height, BMI, or BSA more effectively distill the prognostic information related 
 to PA size most effectively? 
 
We evaluate these associations with a series of Cox-proportional hazard regression
models. The model performance is evaluated using the integrated Brier score. 


* BSA (Mosteller), height_cm and weight_kg: gen bsa = sqrt(height_cm*weight_kg/3600)
(Done in the prep document)


Background: 



Ng et al. proposed the PA:AA ratio as a CT sign of pulmonary arterial hypertension 
- https://pubmed.ncbi.nlm.nih.gov/10524808/ - the ratio correlated with mean 
PA pressure and, unlike absolute PA diameter, showed no association with body surface area (BSA)

A widely cited ATS review stated the aorta serves as an "internal normalization," 
effectively correcting for BSA, sex, age, cardiac‑cycle phase, and technical 
factors—i.e., normalization beyond body size alone. 
https://pmc.ncbi.nlm.nih.gov/articles/PMC4298979/

New population analyses questioned the interpretability of PA:AA when aortic 
caliber varies, arguing that dependence on AA complicates biological meaning and
 risk signals. This reinforces that the ratio is not a body‑size index per se, 
 but a relative vascular‑caliber measure influenced by aortic remodeling. 
 https://pubmed.ncbi.nlm.nih.gov/38637237/
  
 
Some work on generating age, sex, and body-surface-area adjusted expectations of the 
PA and AA size has been done - though z-scores that are directly applicable to
our population aren't available - so we don't utilize these even though it would be nice


*/ 




/*

 Analytical Considerations: 
 
 While this is a study about *prediction* not *cause and effect* - we do care about
 controlling for covariates because we want to capture how much additional 
 information our prognostic factor(s) capture beyond what would be intuitively 
 obviously
 
 -- for example, if PA size just captures age (ie. gets bigger with age, so that
 folks with larger PA die sooner because they are older) that wouldn't be helpful
 to a clinician - so we want to control for age (because the clinician will known
 that before they'd know anything about the PA)
 
 

 
 Unlike (some) prior investigations, we do find that PA:AA correlates (negatively) with age 

*/




/* ---------

goal 1: describe the population

Relevant covariates: 

PA size (probably want to separate on this by tertile, and a summary paragraph)

- mpad
- age
- male 
- height
- weight
- bmi
- death
- time of death

(also comorbidities, for supplement)
---------- */ 

/* 
Table 1
*/ 


table1_mc, by (mpad_tertile) ///
vars( ///
mpad conts %4.1f \ ///
age conts %4.1f \ ///
male bin %4.1f \ /// 
height conts %4.1f \ ///
weight conts %4.1f \ ///
bmi conts %4.1f \ /// 
bsa_mosteller conts %4.2f \ ///
ascendingaorta conts %4.1f \ ///
mpaaa conts %4.2f \ ///
death bin %4.1f \ /// 
time_of_death conts %4.1f \ ///
time_of_censoring conts %4.1f \ ///
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol total(before) ///
saving("`output_dir'/Table 1 Baseline chars by PA tertile.xlsx", replace)



/* 
Supplemental Figure - comorbidities by exposure status
*/ 
table1_mc, by (mpad_tertile) ///
vars( ///
pulmdisease bin %4.1f \ /// 
chf bin %4.1f \ /// 
diabetes3 cat %4.1f \ /// 
hypertension bin %4.1f \ ///
pulmonarycircdisorders bin %4.1f \ /// 
peripheralvascdisorders bin %4.1f \ /// 
renalfailure bin %4.1f \ /// 
anyemergency bin %4.1f \ /// 
anyadmission bin %4.1f \ /// 
death bin %4.1f \ /// 
time_of_death conts %4.1f \ /// 
time_of_censoring conts %4.1f \ /// 
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol total(before) ///
saving("`output_dir'/Table s1 comorbs and utilization by PA tertile.xlsx", replace)


/* How much do ascendingaorta body-size co-vary?*/ 

/* Why do we care? 
If ascending aorta is purely capturing information about person-size (defined by any of the metrics: BSA, BMI, Height)
then, you'd expect that a PA-diameter measure normalized to size would be equally predictive of mortality.

(Similarly, if ascending aorta just captures information about size, age, and sex - then you'd expect normalizing
by those metrics together would be equally predictve) 

However, if Ascending Aorta is not strongly explained by those factors, then there's a potential that it's capturing
other information, that may be more predictive. 
*/

* Simple correlations: 
// how much does variation height explain v variations in ascending aorta?
spearman ascendingaorta height // just assuming ranks correlate; monotone association, strictly increasing/decreasing relationships, not necessarily linear
corr ascendingaorta height //linear association on the original scale
regress ascendingaorta c.height //r2 is the proportion of variance explained 

spearman ascendingaorta weight
ktau    ascendingaorta weight
corr ascendingaorta weight
regress ascendingaorta c.weight

spearman ascendingaorta bsa
ktau    ascendingaorta bsa
corr ascendingaorta bsa
regress ascendingaorta c.bsa

spearman ascendingaorta bmi
ktau    ascendingaorta bmi
corr ascendingaorta bmi
regress ascendingaorta c.bmi


* Partial Pearson controlling for age/sex
//how much does height+weight, BMI, and BSA explain variations in ascending aorta sized
//and whats the unique contribution of each to explain the variability? 
//Partial correlation (r2): Among the variance in Y not explained by the other predictors, the fraction that "X" explains is ___
//Semi-partial correlation (r2): The raw increase in the outcome variation explained that you get when adding X to a model that already contains the other predictors

regress ascendingaorta age i.male
pcorr ascendingaorta age i.male
//age and sex explain 37% of the variation in ascending-aorta size  (adjusted R2)


pcorr ascendingaorta bsa age i.male
pcorr ascendingaorta bmi age i.male

regress ascendingaorta height weight age i.male
pcorr ascendingaorta height weight age i.male
//adding height and weight increases that only to 44% of the variation - adjusted R2 is correct parameter (adjusts for improvement by adding more predictors)

//implication: most of the variation in ascending aorta size is not explained by height, weight, age, and sex (or, it due to non-linear relationships of those)

//thus, if PA:AA is a better predictor of mortality than PA:body size metrics, then it would suggest it is capturing other relevant comorbidity information. 

/* Flexibly modeled */ 

preserve 

* Flexible age for AA modeling too
capture drop age_rc*
mkspline2 age_rc = age, cubic nknots(3)

* Competing body-size parameterizations
regress ascendingaorta c.height c.weight i.male c.(age_rc*), vce(robust)
scalar R2_hw = e(r2)

regress ascendingaorta c.bsa i.male c.(age_rc*), vce(robust)
scalar R2_bsa = e(r2)

regress ascendingaorta c.bmi i.male c.(age_rc*), vce(robust)
scalar R2_bmi = e(r2)

* Partial R² for body-size (example for BSA)
regress ascendingaorta i.male c.(age_rc*), vce(robust)
scalar R2_reduced = e(r2)
regress ascendingaorta c.bsa i.male c.(age_rc*), vce(robust)
scalar R2_full = e(r2)
scalar R2_partial_bsa = (R2_full - R2_reduced) / (1 - R2_reduced)
di "R2(height+weight) = " %4.3f R2_hw " | R2(BSA) = " %4.3f R2_bsa " | R2(BMI) = " %4.3f R2_bmi
di "Partial R2 of BSA given age spline + sex = " %4.3f R2_partial_bsa

capture drop h_rc* w_rc*
mkspline2 h_rc = height, cubic nknots(3)
mkspline2 w_rc = weight, cubic nknots(3)

regress ascendingaorta c.(h_rc* w_rc*) i.male c.(age_rc*), vce(robust)
scalar R2_hw_spline = e(r2)
di "R2(height,weight) spline = " %4.3f R2_hw_spline

restore

















/* part two: adjusted associations */ 

/* ===== Age+sex–adjusted comparison of PA normalizations (raw ratios) ===== */
preserve

* Analytic sample (assumes variables already exist & are cleaned)
drop if missing(mpad, ascendingaorta, age, male, height, weight, bmi, bsa, death, lastfollowupyear)

* Survival setup
stset lastfollowupyear, failure(death==1)

* Age spline (Harrell 10/50/90) and basis
quietly _pctile age, p(10 50 90)
local k_age `r(r1)' `r(r2)' `r(r3)'
capture drop age_rc*
mkspline2 age_rc = age, cubic knots(`k_age')

* ---------- Raw ratio predictors ----------
capture drop r_aa r_h r_w r_bmi r_bsa r_hw pa_hat_hw
gen double r_aa  = mpad/ascendingaorta
gen double r_h   = mpad/height
gen double r_w   = mpad/weight
gen double r_bmi = mpad/bmi
gen double r_bsa = mpad/bsa

* Both-size normalizer: PA divided by raw-units linear prediction from H+W
quietly regress mpad c.height c.weight
predict double pa_hat_hw, xb
replace pa_hat_hw = . if pa_hat_hw<=0
gen double r_hw = mpad/pa_hat_hw

* ---------- Fit each model; collect AIC/BIC/C into locals ----------
local rvars  "r_aa r_h r_w r_bmi r_bsa r_hw"
local names  `" "PA:AA" "PA:Height" "PA:Weight" "PA:BMI" "PA:BSA" "PA:(Ht,Wt)" "'

local i = 0
foreach R of local rvars {
    local ++i
    capture drop rc*
    mkspline2 rc = `R', cubic nknots(4)

    * Baseline: age-spline + sex
    quietly stcox c.(rc*) c.(age_rc*) i.male, vce(robust)
    * If you instead want ht/wt in the baseline, use:
    * quietly stcox c.(rc*) c.(age_rc*) i.male c.height c.weight, vce(robust)

    quietly estat ic
    local A`i' = r(S)[1,5]
    local B`i' = r(S)[1,6]

    local C`i' = .
    capture noisily estat concordance
    if (_rc==0) local C`i' = r(C)
}

* ---------- Build a small results dataset and print ----------
clear
set obs 6
gen str16 model = ""
gen double AIC = . 
gen double BIC = .
gen double C   = .

local j = 0
foreach lbl of local names {
    local ++j
    replace model = `"`lbl'"' in `j'
    replace AIC   = `A`j''     in `j'
    replace BIC   = `B`j''     in `j'
    replace C     = `C`j''     in `j'
}

format AIC %9.3f BIC %9.3f C %6.4f
list, abbreviate(16) noobs sepby(model)

restore














preserve
/* ==== Adjusted PA diameter curve (rcs-4); xblc-safe ==== */

stset lastfollowupyear, failure(death==1)

* Age spline
capture drop age_rc*
mkspline2 age_rc = age, cubic nknots(3)

* Integerize PA diameter to tenths (e.g., 25.1 mm -> 251)
capture drop mpad10 rc*
gen long mpad10 = round(mpad*10)

* 4-knot spline of the predictor
mkspline2 rc = mpad10, cubic nknots(4)

* Full adjustment
stcox c.(rc*) c.(age_rc*) i.male c.height c.weight, vce(robust)
gen byte esample = e(sample)

* Range (p1–p99, clamped to ~19–39 mm)
quietly summarize mpad if esample, detail
local lo = ceil( max(r(p1)*10, 190) )
local hi = floor(min(r(p99)*10, 390) )

* Build at() from observed integer values ONLY
levelsof mpad10 if esample & inrange(mpad10, `lo', `hi'), local(atlist)

* Reference 25.1 mm => 251
xblc rc*, covname(mpad10) at(`atlist') reference(251) eform ///
     generate(pa or lb ub)

* Rescale x back to mm
replace pa = pa/10

twoway ///
  (line lb ub pa, sort lc(black black) lp(longdash longdash)) ///
  (line or pa,   sort lc(black) lp(solid)), ///
  yscale(log extend) ///
  ylabel(1 2 4, angle(horiz) format(%2.1fc) labsize(large)) ///
  xlabel(19(3)39, labsize(large)) xmtick(19(3)39) ///
  ytitle("Hazard Ratio of Mortality", size(large)) ///
  xtitle("Main Pulmonary Artery Diameter (mm)", size(large)) ///
  title("Adjusted Mortality Risk by PA Diameter", size(vlarge)) ///
  note("Model: rcs(PA, 4 on ×10) + rcs(age, 3) + sex + height + weight; robust SEs", size(small)) ///
  yline(1, lp(shortdash) lc(gs10)) ///
  xline(28, lp(shortdash_dot) lc(gs10)) ///
  text(0.85 33.8 "95% CI (black dashes)", size(medlarge)) ///
  text(2.5  24   "Reference = 25.1 mm", size(medlarge)) ///
  scheme(cleanplots)
  
restore

preserve
/* ==== Adjusted Ascending Aorta curve (rcs-4); xblc-safe ==== */

* Survival setup
stset lastfollowupyear, failure(death==1)

* Age spline (3 knots)
capture drop age_rc*
mkspline2 age_rc = age, cubic nknots(3)

* Integerize AA to tenths of a mm (e.g., 30.1 mm -> 301)
capture drop aa10 rc*
gen long aa10 = round(ascendingaorta*10)

* 4-knot spline on the predictor
mkspline2 rc = aa10, cubic nknots(4)

* Full adjustment (identical to heatmap set, minus the 2D PA/AA terms)
stcox c.(rc*) c.(age_rc*) i.male c.height c.weight, vce(robust)
gen byte esample = e(sample)

* Evaluation range: p1–p99, clamped to ~25–40 mm
quietly summarize ascendingaorta if esample, detail
local lo = ceil( max(r(p1)*10, 250) )
local hi = floor(min(r(p99)*10, 400) )

* Build at() from actually observed integer values ONLY
levelsof aa10 if esample & inrange(aa10, `lo', `hi'), local(atlist)

* Reference: AA = 30.1 mm -> 301  (use 301; change if you prefer another ref)
xblc rc*, covname(aa10) at(`atlist') reference(301) eform ///
     generate(aa or lb ub)

* Rescale x back to mm
replace aa = aa/10

* Plot (style matched to your PA and ratio figures)
twoway ///
  (line lb ub aa, sort lc(black black) lp(longdash longdash)) ///
  (line or aa,   sort lc(black)       lp(solid)), ///
  yscale(log extend) ///
  ylabel(1 2 4, angle(horiz) format(%2.1fc) labsize(large)) ///
  xlabel(25(3)40, labsize(large)) xmtick(25(3)40) ///
  ytitle("Hazard Ratio of Mortality", size(large)) ///
  xtitle("Ascending Aorta Diameter (mm)", size(large)) ///
  title("Adjusted Mortality Risk by Ascending Aorta", size(vlarge)) ///
  note("Model: rcs(AA,4 on ×10) + rcs(age,3) + sex + height + weight; robust SEs", size(small)) ///
  yline(1, lp(shortdash) lc(gs10)) ///
  xline(30, lp(shortdash_dot) lc(gs10)) ///
  text(0.85 38.5 "95% CI (black dashes)", size(medlarge)) ///
  text(2.5  30.1 "Reference = 30.1 mm", size(medlarge)) ///
  scheme(cleanplots)

restore


preserve
/* ==== Adjusted PA:AA ratio curve (rcs-4); safe for xblc ==== */

* Survival setup
stset lastfollowupyear, failure(death==1)

* Age spline
capture drop age_rc*
mkspline2 age_rc = age, cubic nknots(3)

* Integerize PA:AA to avoid decimal issues (e.g., 0.77 -> 77)
capture drop paaa100 rc*
gen long paaa100 = round(mpaaa*100)

* 4-knot spline of the predictor
mkspline2 rc = paaa100, cubic nknots(4)

* Full adjustment
stcox c.(rc*) c.(age_rc*) i.male c.height c.weight, vce(robust)
gen byte esample = e(sample)

* Range (p1–p99, clamped)
quietly summarize mpaaa if esample, detail
local lo = ceil( max(r(p1)*100,  60) )
local hi = floor(min(r(p99)*100, 130) )	

* Build at() from observed integer values ONLY
levelsof paaa100 if esample & inrange(paaa100, `lo', `hi'), local(atlist)

* Reference 0.77 => 77
xblc rc*, covname(paaa100) at(`atlist') reference(77) eform ///
     generate(pa or lb ub)

* Rescale x back to ratio
replace pa = pa/100

twoway ///
  (line lb ub pa, sort lc(black black) lp(longdash longdash)) ///
  (line or pa,   sort lc(black) lp(solid)), ///
  yscale(log extend) ///
  ylabel(1 2 4, angle(horiz) format(%2.1fc) labsize(large)) ///
  xlabel(0.6(0.1)1.3, labsize(large)) xmtick(0.6(0.1)1.3) ///
  ytitle("Hazard Ratio of Mortality", size(large)) ///
  xtitle("Pulmonary Artery : Ascending Aorta Ratio", size(large)) ///
  title("Adjusted Mortality Risk by PA:AA", size(vlarge)) ///
  note("Model: rcs(PA:AA, 4 on ×100) + rcs(age, 3) + sex + height + weight; robust SEs", size(small)) ///
  yline(1, lp(shortdash) lc(gs10)) ///
  xline(0.9, lp(shortdash_dot) lc(gs10)) ///
  text(0.85 1.06 "95% CI (black dashes)", size(medlarge)) ///
  text(2.5  0.77 "Reference = 0.77", size(medlarge)) ///
  scheme(cleanplots)

restore






/* Comparisons */ 
/* ==== 2.1 Nested models + information/discrimination ==== */
/* Models:
   M0     : rcs(age,3) + sex + height + weight
   M_PAd  : M0 + rcs(PAd,4) on ×10
   M_Ratio: M0 + rcs(PA:AA,4) on ×100
   M_Both : M0 + rcs(PAd,4) + rcs(PA:AA,4)
*/

preserve 
/* ==== 2.1 Nested models + information/discrimination (no programs) ==== */

* Survival setup
stset lastfollowupyear, failure(death==1)

* Age spline (3 knots)
capture drop age_rc*
mkspline2 age_rc = age, cubic nknots(3)

* Integerized predictors (avoid decimal issues)
capture drop mpad10 paaa100
gen long mpad10  = round(mpad*10)         // 25.1 -> 251
gen long paaa100 = round(mpaaa*100)       // 0.77 -> 77

* Result table
matrix Tab = J(4,3,.)
matrix colnames Tab = AIC BIC C
matrix rownames Tab = M0 M_PAd M_Ratio M_Both

* ---------- M0 ----------
stcox c.(age_rc*) i.male c.height c.weight, vce(robust)

tempname S aic bic cc
quietly estat ic
matrix `S' = r(S)
scalar `aic' = `S'[1,5]
scalar `bic' = `S'[1,6]
capture noisily estat concordance
scalar `cc' = .
if (_rc==0) scalar `cc' = r(C)

matrix Tab[1,1] = scalar(`aic')
matrix Tab[1,2] = scalar(`bic')
matrix Tab[1,3] = scalar(`cc')

* ---------- M_PAd (+ rcs(PAd,4)) ----------
capture drop rcP*
mkspline2 rcP = mpad10, cubic nknots(4)
stcox c.(age_rc*) i.male c.height c.weight c.(rcP*), vce(robust)

quietly estat ic
matrix `S' = r(S)
scalar `aic' = `S'[1,5]
scalar `bic' = `S'[1,6]
capture noisily estat concordance
scalar `cc' = .
if (_rc==0) scalar `cc' = r(C)

matrix Tab[2,1] = scalar(`aic')
matrix Tab[2,2] = scalar(`bic')
matrix Tab[2,3] = scalar(`cc')
capture drop rcP*

* ---------- M_Ratio (+ rcs(PA:AA,4)) ----------
capture drop rcR*
mkspline2 rcR = paaa100, cubic nknots(4)
stcox c.(age_rc*) i.male c.height c.weight c.(rcR*), vce(robust)

quietly estat ic
matrix `S' = r(S)
scalar `aic' = `S'[1,5]
scalar `bic' = `S'[1,6]
capture noisily estat concordance
scalar `cc' = .
if (_rc==0) scalar `cc' = r(C)

matrix Tab[3,1] = scalar(`aic')
matrix Tab[3,2] = scalar(`bic')
matrix Tab[3,3] = scalar(`cc')
capture drop rcR*

* ---------- M_Both (+ both sets of splines) ----------
capture drop rcP* rcR*
mkspline2 rcP = mpad10,  cubic nknots(4)
mkspline2 rcR = paaa100, cubic nknots(4)
stcox c.(age_rc*) i.male c.height c.weight c.(rcP* rcR*), vce(robust)

quietly estat ic
matrix `S' = r(S)
scalar `aic' = `S'[1,5]
scalar `bic' = `S'[1,6]
capture noisily estat concordance
scalar `cc' = .
if (_rc==0) scalar `cc' = r(C)

matrix Tab[4,1] = scalar(`aic')
matrix Tab[4,2] = scalar(`bic')
matrix Tab[4,3] = scalar(`cc')
capture drop rcP* rcR*

* ---------- Report ----------
di as txt "Information / discrimination for nested models:"
matlist Tab, format(%9.4f)

restore

/* 
	•	Adding PA diameter (PAd) to the baseline model gives the biggest, most defensible improvement.
AIC drop ≈19 (decisive), BIC drop ≈5 (positive), C ↑ by ≈0.012.
	•	Adding PA:AA instead of PAd helps less and is not justified by BIC.
AIC improves, but BIC gets worse than baseline, and C ↑ by only ≈0.008.
	•	Adding PA:AA on top of PAd gives negligible discrimination gain and is over‑penalized.
ΔC ≈ +0.001 vs PAd‑only; AIC and especially BIC say don't add it.
	•	If you must choose one 1‑D summary, the evidence here favors PAd over PA:AA.
If you want to use both vascular calibers, use the 2‑D spline surface rather than their ratio.

⸻

3) Manuscript‑ready sentences
	•	"Relative to the adjustment model (age spline, sex, height, weight), adding main PA diameter improved AIC by 19 and BIC by 5 (Harrell's C +0.012), indicating materially better fit and discrimination. Adding the PA:AA ratio improved AIC by 11 but worsened BIC (+2.9) with a smaller gain in C (+0.008), suggesting limited added value under stricter parsimony. Including both PAd and PA:AA increased C by only 0.001 over PAd alone while substantially worsening BIC (+17.6), indicating over‑penalized complexity without meaningful discrimination gain. Overall, PAd is the preferred 1‑D predictor; the ratio adds little once PAd is modeled."
	•	"Akaike weights favored the PAd‑only model (≈0.84) over `Both' (≈0.14) and `Ratio' (≈0.02), implying ≈56:1 support for PAd vs PA:AA among these 1‑D formulations."

⸻

4) What these metrics mean (one‑liners)
	•	AIC/BIC compare relative in‑sample predictive information with a complexity penalty; lower is better. BIC penalizes more strongly; differences >10 are "very strong" evidence.
	•	Harrell's C is discrimination; small increments (≈0.01) are common and should be interpreted with uncertainty (ideally via bootstrap or cross‑validation).
*/ 









/* Part three: 

Guidance from a variety of sources suggests that collapsing a predictor down to 
a ratio is potentially problematic because there are an infinite number of potential 
unique values to get each ratio - and different versions that resolve to the same
ratio may have very different prognostic implications.

Thus, we try to flexibly model both PA-diameter and AA diameter and evaluate the 
pattern of risk to illustrate how at the same PA:AA, different values may have different prognostic 
implications. 
*/ 

/* ===== Heatmap: observed cells only | age spline; sex averaged | height/weight linear at means ===== */

* 0) Snapshot current dataset
tempfile ORIG
save `ORIG', replace

* 1) Survival setup and spline knots on estimation sample
stset lastfollowupyear, failure(death==1)

capture drop mpad_rounded aa_rounded
gen mpad_rounded = round(mpad, 0.1)
gen aa_rounded   = round(ascendingaorta, 0.1)

quietly _pctile mpad_rounded if !missing(mpad_rounded), p(10 50 90)
local k_mpad `r(r1)' `r(r2)' `r(r3)'
quietly summarize mpad_rounded, meanonly
local mp_min = r(min)
local mp_max = r(max)

quietly _pctile aa_rounded if !missing(aa_rounded), p(10 50 90)
local k_aa `r(r1)' `r(r2)' `r(r3)'
quietly summarize aa_rounded, meanonly
local aa_min = r(min)
local aa_max = r(max)

quietly _pctile age if !missing(age), p(10 50 90)
local k_age `r(r1)' `r(r2)' `r(r3)'
quietly summarize age, meanonly
local age_min = r(min)
local age_max = r(max)

* Build spline bases on estimation data & fit model
capture drop mpad_rc* aa_rc* age_rc*
mkspline2 mpad_rc = mpad_rounded, cubic knots(`k_mpad')
mkspline2 aa_rc   = aa_rounded,   cubic knots(`k_aa')
mkspline2 age_rc  = age,          cubic knots(`k_age')

* ► Add linear height and weight
stcox c.(mpad_rc*)##c.(aa_rc*) c.(age_rc*) i.male c.height c.weight, vce(robust)
estimates store M1
estat concordance
gen byte esample = e(sample)

* Means to hold fixed in the grid (computed on the estimation sample)
quietly summarize age    if esample, meanonly
scalar AGEbar = r(mean)
quietly summarize male   if esample, meanonly
scalar PMALE  = r(mean)
quietly summarize height if esample, meanonly
scalar HGTbar = r(mean)
quietly summarize weight if esample, meanonly
scalar WGTbar = r(mean)

* Means of PA and AA in the estimation sample (for reference point)
quietly summarize mpad           if esample, meanonly
scalar MPAbar = r(mean)
quietly summarize ascendingaorta if esample, meanonly
scalar AAbar  = r(mean)
di as txt "Means used:  PA=" %5.2f MPAbar "  AA=" %5.2f AAbar "  (ratio=" %4.2f (MPAbar/AAbar) ")"

* 2) Observed-cell mask from estimation sample (NO hard range)
tempfile EST
save `EST', replace
use `EST', clear
keep if esample
gen mpad_cell = round(mpad, 1)
gen aa_cell   = round(ascendingaorta, 1)
contract mpad_cell aa_cell
rename mpad_cell mpad
rename aa_cell   ascendingaorta
rename _freq     Ncell

* dynamic plotting range from observed cells
quietly summarize mpad, meanonly
local xlo = floor(r(min))     // PA range (for Y axis)
local xhi = ceil(r(max))
quietly summarize ascendingaorta, meanonly
local ylo = floor(r(min))     // AA range (for X axis)
local yhi = ceil(r(max))

* 3) Use OBSERVED CELLS as prediction grid; duplicate for sex; add sentinels
expand 2
bysort mpad ascendingaorta: gen male = _n-1
gen age    = AGEbar
gen height = HGTbar
gen weight = WGTbar

local N0 = _N
set obs `=_N + 6'
* mpad min/max sentinels
replace mpad           = `mp_min' in `= `N0' + 1'
replace ascendingaorta = 30       in `= `N0' + 1'
replace age            = AGEbar   in `= `N0' + 1'
replace male           = 0        in `= `N0' + 1'
replace height         = HGTbar   in `= `N0' + 1'
replace weight         = WGTbar   in `= `N0' + 1'

replace mpad           = `mp_max' in `= `N0' + 2'
replace ascendingaorta = 30       in `= `N0' + 2'
replace age            = AGEbar   in `= `N0' + 2'
replace male           = 0        in `= `N0' + 2'
replace height         = HGTbar   in `= `N0' + 2'
replace weight         = WGTbar   in `= `N0' + 2'

* aa min/max sentinels
replace mpad           = 25       in `= `N0' + 3'
replace ascendingaorta = `aa_min' in `= `N0' + 3'
replace age            = AGEbar   in `= `N0' + 3'
replace male           = 0        in `= `N0' + 3'
replace height         = HGTbar   in `= `N0' + 3'
replace weight         = WGTbar   in `= `N0' + 3'

replace mpad           = 25       in `= `N0' + 4'
replace ascendingaorta = `aa_max' in `= `N0' + 4'
replace age            = AGEbar   in `= `N0' + 4'
replace male           = 0        in `= `N0' + 4'
replace height         = HGTbar   in `= `N0' + 4'
replace weight         = WGTbar   in `= `N0' + 4'

* age min/max sentinels
replace mpad           = 25        in `= `N0' + 5'
replace ascendingaorta = 30        in `= `N0' + 5'
replace age            = `age_min' in `= `N0' + 5'
replace male           = 0         in `= `N0' + 5'
replace height         = HGTbar    in `= `N0' + 5'
replace weight         = WGTbar    in `= `N0' + 5'

replace mpad           = 25        in `= `N0' + 6'
replace ascendingaorta = 30        in `= `N0' + 6'
replace age            = `age_max' in `= `N0' + 6'
replace male           = 0         in `= `N0' + 6'
replace height         = HGTbar    in `= `N0' + 6'
replace weight         = WGTbar    in `= `N0' + 6'

gen byte sentinel = (_n > `N0')

* 4) Predict on observed cells; average over sex; normalize
gen mpad_rounded = round(mpad, 0.1)
gen aa_rounded   = round(ascendingaorta, 0.1)
mkspline2 mpad_rc = mpad_rounded, cubic knots(`k_mpad')
mkspline2 aa_rc   = aa_rounded,   cubic knots(`k_aa')
mkspline2 age_rc  = age,          cubic knots(`k_age')

estimates restore M1
predict xb, xb
gen HR = exp(xb)

drop if sentinel==1   // sentinels served their purpose

gen w = cond(male==1, PMALE, 1-PMALE)
bysort mpad ascendingaorta: egen HRw  = total(HR*w)
bysort mpad ascendingaorta: egen wsum = total(w)
gen HR_avg = HRw/wsum

* reference = nearest available observed cell to the JOINT MEANS (MPAbar, AAbar)
gen double d2 = (mpad - MPAbar)^2 + (ascendingaorta - AAbar)^2
bysort mpad ascendingaorta: keep if _n==1     // ensure unique grid cells
sort d2
scalar HRREF        = HR_avg[1]
local  ref_m_used   = mpad[1]
local  ref_a_used   = ascendingaorta[1]
gen double normalized_HR = HR_avg/HRREF

* 5) Bin + label legend categories (explicit labels; drop code 0 if present)
capture drop hr_bin
local cuts "0.55 0.67 0.83 0.91 1 1.10 1.25 1.50 2 3 4 6 9 12 16.5"
egen hr_bin = cut(normalized_HR), at(`cuts') icodes

label define hr_lab 1 "≤0.67" 2 "0.67–0.83" 3 "0.83–0.91" 4 "0.91–1.00" ///
                    5 "1.00–1.10" 6 "1.10–1.25" 7 "1.25–1.50" 8 "1.50–2.00" ///
                    9 "2.00–3.00" 10 "3.00–4.00" 11 "4.00–6.00" 12 "6.00–9.00" ///
                    13 "9.00–12.00" 14 "≥12.00", replace
label values hr_bin hr_lab

* Build legend labels ONLY for bins that occur (exclude code 0)
levelsof hr_bin if hr_bin>0 & !missing(hr_bin), local(usedbins)

local ramp_labels
foreach k of local usedbins {
    * Ask for the label via the variable that carries it (hr_bin), not the label name
    local lab : label (hr_bin) `k'
    local ramp_labels `ramp_labels' `k' "`lab'"
}
local nbins : word count `usedbins'

* Ratio guide lines and reference marker
local overlays
foreach r in 0.7 0.8 0.9 1.0 1.1 {
    local overlays `overlays' || ///
        function y = `r'*x, range(`ylo' `yhi') ///
        lcolor(black%35) lpattern(shortdash) lwidth(thin) legend(off)
}
local overlays `overlays' || ///
    scatter mpad ascendingaorta if mpad==`ref_m_used' & ascendingaorta==`ref_a_used', ///
    msymbol(plus) mcolor(black) msize(medlarge) legend(off)

* 6) Heatmap: X=AA, Y=PA; categorical legend; split title across two lines
heatplot hr_bin mpad ascendingaorta, discrete ///
    color(RdYlBu, reverse n(`nbins')) ///
    ramp(right subtitle("HR (vs ref)") label(`ramp_labels')) ///
    title("Predicted Hazard Ratio of Death (vs mean)") ///
    subtitle("AA and PA + interaction (age at mean; sex averaged; ht/wt at means)", size(medsmall)) ///
    xtitle("Ascending Aorta Diameter (mm)") ///
    ytitle("Main PA Diameter (mm)") ///
    xlabel(`ylo'(5)`yhi') ylabel(`xlo'(5)`xhi') ///
    scheme(white_w3d) ///
    addplot(`overlays')

	
	
* Verify the marker location and HR normalization
di "Ref cell used:  PA=" `ref_m_used' "  AA=" `ref_a_used'
list normalized_HR if mpad==`ref_m_used' & ascendingaorta==`ref_a_used'
* should be 1

* Compare mean ratio to the guide lines
di "Mean PA:AA ratio = " %5.3f (MPAbar/AAbar)



* 7) Restore your original dataset
use `ORIG', clear
/* ===== End block ===== */




/* ==== 3.1 Iso‑ratio HR curves from the 2D surface (fixed age knots) ==== */
/* Fits 2D model: rcs(PAd,3) × rcs(AA,3) + rcs(age,3) + sex + height + weight
   Predicts along PA = r * AA, r ∈ {0.7,0.8,0.9,1.0,1.1}
   Normalizes to HR at (MPAbar, AAbar) nearest grid point. */

preserve

* --- Fit 2D model on analytic sample ---
stset lastfollowupyear, failure(death==1)

* Data-driven knots from the estimation data
quietly _pctile mpad,            p(10 50 90)
local k_mpad `r(r1)' `r(r2)' `r(r3)'
quietly _pctile ascendingaorta,  p(10 50 90)
local k_aa   `r(r1)' `r(r2)' `r(r3)'
quietly _pctile age,             p(10 50 90)
local k_age  `r(r1)' `r(r2)' `r(r3)'

* Build spline bases & fit (use the same knots)
capture drop age_rc* mpad_rc* aa_rc*
mkspline2 age_rc  = age,           cubic knots(`k_age')
mkspline2 mpad_rc = mpad,          cubic knots(`k_mpad')
mkspline2 aa_rc   = ascendingaorta, cubic knots(`k_aa')

stcox c.(mpad_rc*)##c.(aa_rc*) c.(age_rc*) i.male c.height c.weight, vce(robust)
estimates store M2D
gen byte __esmp = e(sample)

* Holds at estimation-sample means
quietly summarize age    if __esmp, meanonly
scalar AGEbar = r(mean)
quietly summarize height if __esmp, meanonly
scalar HGTbar = r(mean)
quietly summarize weight if __esmp, meanonly
scalar WGTbar = r(mean)
quietly summarize mpad           if __esmp, meanonly
scalar MPAbar = r(mean)
quietly summarize ascendingaorta if __esmp, meanonly
scalar AAbar  = r(mean)

* AA range for the iso-ratio traces (p10–p90 of AA)
quietly _pctile ascendingaorta if __esmp, p(10 90)
local AAlo = floor(r(r1))
local AAhi = ceil(r(r2))

* --- Build iso‑ratio grid ---
tempfile iso
tempname pst
postfile `pst' double ascendingaorta mpad ratio using `iso', replace
local ratios "0.7 0.8 0.9 1.0 1.1"
forvalues A = `AAlo'/`AAhi' {
    foreach r of local ratios {
        post `pst' (`A') (`= `r'*`A' ') (`r')
    }
}
postclose `pst'
use `iso', clear

gen double ratio_key = round(ratio, .001)   // 0.700, 0.800, ... robust to fp error

* Set display covariates
gen age    = AGEbar
gen height = HGTbar
gen weight = WGTbar
gen male   = 0     // set to 1 for male; or duplicate/average if desired

* Recreate spline bases on the grid **using the same knots** (no data-driven knots)
mkspline2 age_rc  = age,           cubic knots(`k_age')
mkspline2 mpad_rc = mpad,          cubic knots(`k_mpad')
mkspline2 aa_rc   = ascendingaorta, cubic knots(`k_aa')

* Predict and normalize
estimates restore M2D
predict double xb, xb
gen double HR = exp(xb)

gen double d2 = (mpad-MPAbar)^2 + (ascendingaorta-AAbar)^2
sort d2
scalar HRREF = HR[1]
gen double HRnorm = HR/HRREF

* Plot
twoway ///
    (line HRnorm ascendingaorta if ratio_key==0.70, sort) ///
    (line HRnorm ascendingaorta if ratio_key==0.80, sort) ///
    (line HRnorm ascendingaorta if ratio_key==0.90, sort) ///
    (line HRnorm ascendingaorta if ratio_key==1.00, sort) ///
    (line HRnorm ascendingaorta if ratio_key==1.10, sort), ///
    yscale(log) ylabel(0.67 1 1.5 2 3 4, angle(horiz)) ///
    legend(order(1 "r=0.7" 2 "r=0.8" 3 "r=0.9" 4 "r=1.0" 5 "r=1.1")) ///
    xtitle("Ascending Aorta (mm)") ///
    ytitle("HR (normalized to joint means)") ///
    title("Predicted HR along iso-ratio lines  (PA = r × AA)") ///
    scheme(cleanplots)

/* So every curve shows relative hazard versus an "average" patient at the same held covariates. This removes the baseline hazard and any multiplicative constant, making curves comparable and visually anchored at 1.0. */ 
	
* Clean
capture drop __esmp age_rc* mpad_rc* aa_rc*

restore







/* ==== 3.3 Model comparison: ratio vs 2D surface ==== */
/* Both adjusted for rcs(age,3) + sex + height + weight (robust). */
preserve 
stset lastfollowupyear, failure(death==1)

* Age spline
capture drop age_rc*
mkspline2 age_rc = age, cubic nknots(3)

* --- M_Ratio: rcs(PA:AA,4) on ×100 ---
capture drop paaa100 rcR*
gen long paaa100 = round(mpaaa*100)
mkspline2 rcR = paaa100, cubic nknots(4)
stcox c.(age_rc*) i.male c.height c.weight c.(rcR*), vce(robust)

tempname S aR bR cR
quietly estat ic
matrix `S' = r(S)
scalar `aR' = `S'[1,5]
scalar `bR' = `S'[1,6]
capture noisily estat concordance
scalar `cR' = .
if (_rc==0) scalar `cR' = r(C)
capture drop rcR*

* --- M_2D: rcs(PA,3) × rcs(AA,3) ---
quietly _pctile mpad, p(10 50 90)
local k_mpad `r(r1)' `r(r2)' `r(r3)'
quietly _pctile ascendingaorta, p(10 50 90)
local k_aa   `r(r1)' `r(r2)' `r(r3)'

capture drop mpad_rc* aa_rc*
mkspline2 mpad_rc = mpad,          cubic knots(`k_mpad')
mkspline2 aa_rc   = ascendingaorta, cubic knots(`k_aa')

stcox c.(mpad_rc*)##c.(aa_rc*) c.(age_rc*) i.male c.height c.weight, vce(robust)

tempname S2 a2 b2 c2
quietly estat ic
matrix `S2' = r(S)
scalar `a2' = `S2'[1,5]
scalar `b2' = `S2'[1,6]
capture noisily estat concordance
scalar `c2' = .
if (_rc==0) scalar `c2' = r(C)

* Report table
matrix M = ( scalar(`aR') , scalar(`bR') , scalar(`cR') \ ///
             scalar(`a2') , scalar(`b2') , scalar(`c2') )
matrix rownames M = M_Ratio M_2D
matrix colnames M = AIC BIC C
di as txt "Model comparison: Ratio (1D) vs 2D surface"
matlist M, format(%9.4f)

* Clean
capture drop mpad_rc* aa_rc* age_rc* paaa100

restore

