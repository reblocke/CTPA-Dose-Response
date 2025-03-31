import sys
import sfi
import numpy as np
import pandas as pd
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.metrics import integrated_brier_score

# Retrieve macros passed as command-line arguments:
# Expecting: sttime, ststatus, modelspec
sttime_var = sys.argv[1]
ststatus_var = sys.argv[2]
modelspec = sys.argv[3]
covariates = modelspec.split()

# Retrieve data from Stata
time_data  = sfi.Data.get(sttime_var)
event_data = sfi.Data.get(ststatus_var)

# Debug prints:
print("DEBUG: Type of time_data[0]:", type(time_data[0]), "Value:", time_data[0])
print("DEBUG: Type of event_data[0]:", type(event_data[0]), "Value:", event_data[0])

def flatten(val):
    if isinstance(val, list):
        return val[0]
    return val

# Build data for predictors
cov_data = {}
for var in covariates:
    cov_data[var] = sfi.Data.get(var)

df = pd.DataFrame({sttime_var: time_data, ststatus_var: event_data})
for var in covariates:
    df[var] = cov_data[var]
df = df.dropna()

# Flatten the values if needed and construct the structured array
structured_y = np.array([(bool(flatten(e)), float(flatten(t))) 
                          for e, t in zip(df[ststatus_var], df[sttime_var])],
                          dtype=[('event', '?'), ('time', '<f8')])

X = df[covariates]

# Fit Cox model
cox_model = CoxPHSurvivalAnalysis()
cox_model.fit(X, structured_y)

# Create grid of times from 0 to 12 years
times_grid = np.linspace(0, 12, 100)

# Predict survival functions and compute integrated Brier score
pred_surv_functions = cox_model.predict_survival_function(X)
pred_surv = np.vstack([fn(times_grid) for fn in pred_surv_functions])
ibs_value = integrated_brier_score(structured_y, structured_y, pred_surv, times_grid)

print("Integrated Brier Score (IBS) up to year 12 = ", ibs_value)