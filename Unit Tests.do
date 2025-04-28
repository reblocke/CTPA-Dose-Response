do define_calc_ibs.do

clear
set seed 123456   // for reproducibility
set obs 500       // 500 subjects

gen id = _n

// Generate event times from Exponential(rate=0.1)
gen time = -ln(runiform()) / 0.1

// Administrative censoring at 12 years
gen death = 1
replace death = 0 if time > 12
replace time  = 12 if time > 12

rename time lastfollowupyear

// Define st
stset lastfollowupyear, failure(death==1)

// Optionally, fit a no-covariate Cox model 
// (just to illustrate):
// stcox
// Then run your calc_ibs program, e.g.:
// calc_ibs " " , graph
//
// Because there's really no primary predictor, you might call:
// calc_ibs "age", graph
// if you wanted to see how the code runs, but in truth there's no real predictor 
// in the data. You can just leave it blank or pick a dummy variable.

stcox
calc_ibs " " , graph





clear
set seed 654321
set obs 500       // 500 subjects

gen id = _n

// Assign 50% male vs female
gen male = (_n <= 250)

// For male=0: hazard=0.1; for male=1: hazard=0.2
// So if male=1, rate=0.2; if male=0, rate=0.1
gen basehaz = cond(male==1, 0.2, 0.1)

// Simulate Exponential with that subject-specific rate
gen time = -ln(runiform()) / basehaz

// Administrative censoring at 12 years
gen death = 1
replace death = 0 if time > 12
replace time  = 12 if time > 12

rename time lastfollowupyear

// stset
stset lastfollowupyear, failure(death==1)

// Now fit the correct Cox model (one predictor: male):
// stcox male
//
// Then call your calc_ibs program, e.g.:
// calc_ibs "male" , graph
//
// The IBS you get should hover near ~4.7, subject to simulation error. 

stcox male
predict lp, xb
predict bs, basesurv
gen s = bs^exp(lp)
stcurve, survival
calc_ibs male








capture program drop calc_ibs
program define calc_ibs
    version 17.0
    syntax varlist(min=1)
    local modelspec "`varlist'"
    local sttime "`e(timevar)'"
    local ststatus "`e(failvar)'"
    
    // Call the external Python script with the proper arguments
    python script calc_ibs.py, args("`sttime'" "`ststatus'" "`modelspec'")
end

calc_ibs mpad age male

/*
program define calc_ibs, rclass
	syntax anything(name=modelspec), [graph]

	********************************************************************************
	* Calculate integrated and normalized Brier score for a given Cox model
	********************************************************************************

	preserve

	* Tokenize model spec: first variable = primary predictor
	tokenize `modelspec'
	local primary_predictor `1'
	macro shift
	local secondary_predictors `*'

	di as text "Primary Predictor: `primary_predictor'"
	di as text "Secondary Predictors: `secondary_predictors'"

	* Time points (0–12 years, yearly)
	local times 1 2 3 4 5 6 7 8 9 10 11 12

	* Fixed IPCW adjustment
	local ipcw_adjustment c.age i.male //i.obesity - has some missing values

	********************************************************************************
	* Compute Brier scores for each time point
	********************************************************************************
	foreach t of local times {
		stbrier `primary_predictor' `secondary_predictors', ///
			ipcw(`ipcw_adjustment') ///
			btime(`t') ///   		
			gen(br`t')
	}

	* Collapse to mean Brier score per time
	collapse (mean) br*

	* Add a row for time=0 with Brier=0
	set obs `=_N + 1'
	foreach t of local times {
		local varname = "br`t'"
		replace `varname' = 0 in l
	}
	gen order = _n
	replace order = 0 in l
	sort order

	* Reshape long
	gen id = _n
	reshape long br, i(id) j(_t)
	sort _t

	* Trapezoidal rule to calculate integrated Brier score
	gen t0 = _t[_n-1]
	gen br0 = br[_n-1]
	gen interval = _t - t0
	gen area = 0.5 * (br + br0) * interval
	replace area = . if _n == 1

	summarize area, meanonly
	scalar IBS = r(sum)
	di as result "Integrated Brier Score (0..12): " IBS
	return scalar IBS = IBS

	gen normalized_area = area / 12
	replace normalized_area = . if _n == 1

	summarize normalized_area, meanonly
	scalar normalized_IBS = r(sum)
	di as result "Normalized IBS (0..12): " normalized_IBS
	return scalar normalized_IBS = normalized_IBS

	* Optional: Plot Brier score curve
	if "`graph'" != "" {
		line br _t if _t <= 12, sort ytitle("Brier Score") ///
			xtitle("Time (years)") ///
			title("Brier Score by Time") ///
			scheme(white_w3d)
	}

	restore
end
*/ 



******************************************************
* 1) Example: Setting up survival in Stata
******************************************************
*    You mentioned:
*      stset lastfollowupyear, failure(death == 1)
*      stcox mpad c.age i.male
*
*    We won't actually run stcox here; just showing that
*    you already have a time variable (lastfollowupyear),
*    an event variable (death==1), and covariates (mpad, age, male).
******************************************************

* If not already done, ensure the data is survival-time ready:
stset lastfollowupyear, failure(death == 1)

***************************************************************
* 2) Python block: pass data from Stata to Python, fit a Cox
*    model via scikit-survival, and calculate IBS up to 12 years
***************************************************************

python:
import sfi  # sfi = Stata Function Interface
import numpy as np
import pandas as pd

# scikit-survival
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.metrics import integrated_brier_score

#
# 2a) Extract data from Stata
#

# Let's read the relevant variables into Python.
# We will create a pandas DataFrame for convenience.

# Pull variables from current Stata dataset
time_var    = 'lastfollowupyear'
event_var   = 'death'
covariates  = ['mpad', 'age', 'male']  # the same as stcox mpad c.age i.male (male is factor-coded in Stata, but let's keep it simple)

# Get all data as columns
time_data  = sfi.Data.get(time_var)
event_data = sfi.Data.get(event_var)

# Covariates come as a list of lists
cov_data = {var: sfi.Data.get(var) for var in covariates}

# Convert to pandas DataFrame
df = pd.DataFrame({
    time_var:  time_data,
    event_var: event_data
})
for var in covariates:
    df[var] = cov_data[var]

# Drop any missing if needed (common in real data)
df = df.dropna()

# Convert the event/time into the structured array format required by sksurv
# sksurv expects a structured array with dtype=[('event', bool), ('time', float)] 
structured_y = np.array([
    (bool(e), t) for e, t in zip(df[event_var], df[time_var])
    ],
    dtype=[('event', '?'), ('time', '<f8')]
)

# Covariates matrix
X = df[covariates]

#
# 2b) Fit a Cox model using scikit-survival
#
cox_model = CoxPHSurvivalAnalysis()
cox_model.fit(X, structured_y)

#
# 2c) Compute the Brier score over time, then integrate up to year = 12
#
# We need a grid of times at which to evaluate the predicted survival. 
# If your data is in *years*, you can go from 0 to 12 in small steps:
times_grid = np.linspace(0, 12, 100)  # 100 points up to year 12

# Predict the survival function for each subject in X
pred_surv_functions = cox_model.predict_survival_function(X)

# `pred_surv_functions` is a list of callables; evaluate them at times_grid
# producing a 2D array: predictions[sample, time_index]
pred_surv = np.row_stack([
    fn(times_grid) for fn in pred_surv_functions
])

# Now compute the IBS. For an "apparent" estimate, we can pass the same dataset
# as both "training" and "testing" in integrated_brier_score. 
# Real usage might do cross-validation for a more honest IBS.
from sksurv.metrics import brier_score

ibs_value = integrated_brier_score(structured_y, structured_y, pred_surv, times_grid)

# Print or return IBS so we can capture it in Stata
print("Integrated Brier Score (IBS) up to year 12: ", ibs_value)

end
