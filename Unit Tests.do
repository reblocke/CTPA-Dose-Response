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
