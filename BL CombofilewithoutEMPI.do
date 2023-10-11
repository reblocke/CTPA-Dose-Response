
*========================*
*----------begin analysis------*
cd "/Users/blocke/Box Sync/Residency Personal Files/Scholarly Work/Locke Research Projects/Pulm Artery Stuff w Scarps/Local Analysis"

//cd "C:\Users\reblo\Box\Residency Personal Files\Scholarly Work\Locke Research Projects\Pulm Artery Stuff w Scarps\Local Analysis"

capture mkdir "Results and Figures"
capture mkdir "Results and Figures/$S_DATE/" //make new folder for figure output if needed
capture mkdir "Results and Figures/$S_DATE/Logs/" //new folder for stata logs
local a1=substr(c(current_time),1,2)
local a2=substr(c(current_time),4,2)
local a3=substr(c(current_time),7,2)
local b = "BL CombofilewithoutEMPI.do" // do file name
copy "`b'" "Results and Figures/$S_DATE/Logs/(`a1'_`a2'_`a3')`b'"

use final_noempi, clear
count // n=990
count if lastfollowupday==0 // n=78 without any follow-up data

//ssc install schemepack, replace //if you hven't
//net install cleanplots, from("https://tdmize.github.io/data/cleanplots")
set scheme cleanplots //cleanplots white_tableau white_w3d //3 options I usually use

/* Data cleaning */ 

//ssc install missings, replace // if you don't have already
missings report
missings dropobs, force //drop all rows with no observations. None
missings dropvars, force //drop all columns with no observations
missings report
mdesc
codebook

label variable age "Age (years)"

recode age min/30=0 30/40=1 40/50=2 50/60=3 60/70=4 70/max=5, gen(age_decade)
label define age_dec_lab 0 "<30 years" 1 "30-40 years" 2 "40-50 years" 3 "50-60 years" 4 "60-70 years" 5 "70+ years"
label variable age_decade "Age (by decade)"
label values age_decade age_dec_lab 

label define male_lab 0 "Female" 1 "Male"
label variable male "Male Sex?"
label values male male_lab

bysort obesity: sum bmi, detail
label variable obesity "Obesity (ICD only)"

generate obesity_calc = . 
replace obesity_calc = 0 if !missing(obesity) & obesity == 0
replace obesity_calc = 1 if !missing(obesity) & obesity == 1
replace obesity_calc = 1 if !missing(bmi) & bmi >= 30
replace obesity_calc = 0 if !missing(bmi) & bmi < 30
label define obesity_calc_lab 0 "Not Obese" 1 "Obese"
label variable obesity_calc "Obesity"
label values obesity_calc obesity_calc_lab
tab obesity obesity_calc, missing
tab obesity_calc, missing
tab obesity, missing

label variable enlargedratio "PA:AA Ratio"
label define enlargedratio_lab 0 "Normal PA:AA (<0.9)" 1 "Increased PA:AA (0.9+)"
label values enlargedratio enlargedratio_lab

label variable enlargedpa "PA diameter"
label define enlargedpa_lab 0 "Normal PAd" 1 "Enlarged PAd"
label values enlargedpa enlargedpa_lab

label variable pulmdisease "Pulmonary Disease"
label variable chf "Congestive Heart Failure"
label variable diabetes3 "Diabetes"
label define diabetes_val 0 "No DM" 1 "Uncomplicated DM" 2 "DM with Complication(s)"
label values diabetes3 diabetes_val 
label variable hypertension "Hypertension"
label variable pulmonarycircdisorder "Pulm Circ. Disorder"
label variable peripheralvascdisorders "Periph Vasc. Disease"
label variable renalfailure "Kidney Disease"

gen known_assoc_comorb = 0
replace known_assoc_comorb = 1 if (pulmdisease == 1 | pulmonarycircdisorders == 1 | chf == 1) 
tab chf known_assoc_comorb, missing //sanity checks
tab pulmdisease known_assoc_comorb, missing
tab pulmonarycircdisorders known_assoc_comorb, missing
label variable known_assoc_comorb "Known Comorbidity w PA->Mortality Association?"
label define known_assoc_comorb_lab 0 "No Pulm Dz/CHF/Pulm Vasc Dz" 1 "Any of Pulm Dz/CHF/Pulm Vasc Dz"
label values known_assoc_comorb known_assoc_comorb_lab

// mean of each - splitting into complicated and uncomplicated. 
tab diabetes3
capture drop diabetes1
gen diabetes1=cond(diabetes3==1,1,0)
tab diabetes3 diabetes1, missing
capture drop diabetes2
gen diabetes2=cond(diabetes3==2,1,0)
tab diabetes3 diabetes2, missing

label variable diabetes1 "Uncomplicated Diabetes"
label variable diabetes2 "Diabetes with Complication(s)"

gen pa_confusion_matrix = .
replace pa_confusion_matrix = 0 if (enlargedpa == 0) & (enlargedratio == 0)
replace pa_confusion_matrix = 1 if (enlargedpa == 1) & (enlargedratio == 0)
replace pa_confusion_matrix = 2 if (enlargedpa == 0) & (enlargedratio == 1)
replace pa_confusion_matrix = 3 if (enlargedpa == 1) & (enlargedratio == 1)
label variable pa_confusion_matrix "PAd high, PA:AA high, neither, or both?"
label define pa_confusion_lab 0 "Neither" 1 "Only Enlarged PAd" 2 "Only High PA:AA" 3 "Both"
label values pa_confusion_matrix pa_confusion_lab
tab enlargedpa pa_confusion_matrix, col //sanity checks
tab enlargedratio pa_confusion_matrix 

label variable death "Death"
label define death_label 0 "Alive or Censored" 1 "Died"
label values death death_label

label variable anyemergency "Any Emergency Visits in Follow-up?"
label variable anyadmission "Any Hospital Admission in Follow-up?"

label variable lastfollowupyear "Follow-up or death(years)"

//mpaaa 0.9 is threshold by framingham: break out into more granular categories
recode mpaaa min/0.75=0 0.75/0.9=1 0.9/1.05=2 1.05/max=3, gen(mpaaa_cat)
label define mpaaa_cat_lab 0 "<0.75" 1 "0.75 to 0.9" 2 "0.9 to 1.05" 3 "1.05+"
//TODO:  Cutpoints for ratio PA were defined as normal ≤0.9; mild>0.9 to 1.0; moderate>1.0 to 1.1; and severe>1.1  Truong QA, Bhatia HS, Szymonifka J, et al. A four-tier classification system of pulmonary artery metrics on computed tomography for the diagnosis and prognosis of pulmonary hypertension. J Cardiovasc Comput Tomogr 2018;12(1):60–66.
label variable mpaaa_cat "PA:AA strata"
label values mpaaa_cat mpaaa_cat_lab 

gen sex_norm_mpad = .
replace sex_norm_mpad = mpad + 1 if male == 0
replace sex_norm_mpad = mpad - 1 if male == 1
sum sex_norm_mpad, detail
recode sex_norm_mpad min/25=0 25/28=1 28/31=2 31/max=3, gen(mpad_cat) 
label define mpad_cat_lab 0 "F:<24 mm, M:<26 mm" 1 "F:24-27 mm, M:26-29 mm" 2 "F:27-30 mm, M:29-32 mm" 3 "F:30+ mm, M32+ mm"
//TODO: Cutpoints for mPA were defined with ≤27 mm(F) and ≤29 mm(M) as the normal reference range; mild as >27 to <31 mm(F) and >29 to <31 mm(M); moderate≥31–34 mm; and severe>34 mm.  from -  Truong QA, Bhatia HS, Szymonifka J, et al. A four-tier classification system of pulmonary artery metrics on computed tomography for the diagnosis and prognosis of pulmonary hypertension. J Cardiovasc Comput Tomogr 2018;12(1):60–66.
label variable mpad_cat "PA_d strata"
label values mpad_cat mpad_cat_lab 
bysort mpad_cat: sum sex_norm_mpad, detail
bysort mpad_cat: sum mpad if male == 0, detail  //sanity check
bysort mpad_cat: sum mpad if male == 1, detail  //sanity check

gen sex_norm_aa = .
replace sex_norm_aa = ascendingaorta + 1.5 if male == 0 //men 3mm larger; more dispersion (not accounted for)
replace sex_norm_aa = ascendingaorta - 1.5 if male == 1
sum sex_norm_aa, detail
recode sex_norm_aa min/30=0 30/33.5=1 33.5/37=2 37/max=3, gen(aa_cat)
label define aa_cat_lab 0 "F:<28.5 mm, M:<31.5 mm" 1 "F:28.5-32 mm, M:31.5-35 mm" 2 "F:32-35.5 mm, M:35-38.5mm" 3 "F:35.5+mm, M:38.5+mm" //35 and 32 as thresholds (similar percentil as for MPAd) 33.5 as threshold for abn, roughly
label variable aa_cat "AA_d strata"
label values aa_cat aa_cat_lab 
bysort aa_cat: sum sex_norm_aa, detail
bysort aa_cat: sum ascendingaorta if male == 0, detail
bysort aa_cat: sum ascendingaorta if male == 1, detail 

recode sex_norm_aa min/33.5=0 33.5/max=1, gen(enlargedaa)
label define aa_abn_lab 0 "AA F:<32 mm M:<35mm" 1 "AA F:32+ mm M:35+ mm"
label variable enlargedaa "AA_d above ULN? (32mm in F, 35mm in M)"
label values enlargedaa aa_abn_lab
tab3way enlargedaa aa_cat male

gen pa_aa_confusion_matrix = .
replace pa_aa_confusion_matrix = 0 if (enlargedpa == 0) & (enlargedaa == 0)
replace pa_aa_confusion_matrix = 1 if (enlargedpa == 1) & (enlargedaa == 0)
replace pa_aa_confusion_matrix = 2 if (enlargedpa == 0) & (enlargedaa == 1)
replace pa_aa_confusion_matrix = 3 if (enlargedpa == 1) & (enlargedaa == 1)
label variable pa_aa_confusion_matrix "PAd high, AAd high, neither, or both?"
label define pa_aa_confusion_lab 0 "Neither" 1 "Only Enlarged PAd" 2 "Only High AAd" 3 "Both"
label values pa_aa_confusion_matrix pa_aa_confusion_lab
tab enlargedpa pa_aa_confusion_matrix //sanity checks
tab enlargedaa pa_aa_confusion_matrix 

recode age min/50=0 50/65=1 65/max=2, gen(age_bi_decade)
label define age_bi_dec_lab 0 "<50 years" 1 "50-65 years" 2 "65+ years"
label variable age_bi_decade "Age (by 15 yr)"
label values age_bi_decade age_bi_dec_lab 
bysort age_bi_decade: sum age, detail //sanity check

save cleaned_noempi, replace

/* Exploring Data */ 

clear
use cleaned_noempi


/* Univariate explorations of exposure, confounders, and outcomes */ 

// -- exposure
/*
//Evaluate distributions of these
bysort male: sum mpaaa, detail //not hugely different
tab mpaaa_cat, plot
tab age_decade mpaaa_cat, row

bysort male: sum mpad, detail //quite a bit different; men 2mm larger; *more dispersion*
tab mpad_cat, plot
tab age_decade mpad_cat
tab mpaaa_cat mpad_cat

bysort male: sum ascendingaorta, detail //quite a bit different; men 3mm larger; more dispersion
//35 and 32 as thresholds (similar percentil as for MPAd)
tab aa_cat, plot
tab age_decade aa_cat
tab mpaaa_cat aa_cat 
tab3way  mpad_cat aa_cat mpaaa_cat

// prevalence rates for PA_d and PA:AA by key confounders
tab age_decade enlargedratio, row // much more likely in <30 and 30-40 to be enlarged
tab obesity_calc enlargedratio, row // obese less likely to have enlarged ratio
poisson enlargedratio male obesity_calc ib3.age_decade, irr robust //after adjusting for age and sex, obesity is MORE likely to have elevated ratio. Sex not indep associated

tab age_decade enlargedpa, row // much more likely in 60-70 and 70+ to be enlarged
tab obesity_calc enlargedpa, row // obese more likely to have enlarged ratio
poisson enlargedpa male obesity_calc ib3.age_decade, irr robust //after adjusting for age, obesity is strongly associated (so is age). Sex not indep associated

//ssc install tab3way, replace //if you don't have already
tab3way enlargedpa enlargedratio age_decade //among <30 - many have increased ratio, normal PAd. Amongst 70+ - many have enlarged PAd, normal ratio

tab enlargedpa enlargedratio, row
tab enlargedpa enlargedratio, col

tab3way enlargedpa enlargedratio obesity_calc //those with only increased ratio (not PAd) are more common in non-obese.
tab3way enlargedpa enlargedratio male //those with only increased ratio (not PAd) are more common in women. May be confounding by age and bmi

bysort male: tab enlargedpa enlargedratio, row
diagt enlargedratio enlargedpa 
diagt enlargedratio enlargedpa if male == 0
diagt enlargedratio enlargedpa if male == 1
logistic enlargedratio enlargedpa male ib2.age_decade

tab3way enlargedpa enlargedratio known_assoc_comorb //those w NO CHF/PAH/Pulm less likely to have just PAd enlarged than just PA:AA. If CHF/PAH/or pulm - balanced.

*/ 

// -- confounders
/*
tab age_decade male, row // gender balance is relatively even by decade

tab age_decade obesity, row //less obese patients at extreme of age
bysort obesity: sum bmi, detail // roughly ~30% without obesity ICD have BMI over 30, roughly 15% with label have BMI < 30

tab age_decade obesity_calc, row //more obesity if you add info from recorded BMI
bysort obesity_calc: sum bmi, detail
tab male obesity_calc, chi2 //no diff by sex

tab age_decade pulmdisease, row //pulm dz less common in <30
tab male pulmdisease, row // no diff by sex
tab obesity_calc pulmdisease, row //more pulm dz in obese

tab age_decade chf, row //lots more chf in olds
tab male chf, row //similar ish
tab obesity_calc chf, row chi2 //more chf in obese

tab age_decade pulmonarycircdisorders, row //pretty stable by age
tab male pulmonarycircdisorders, row //stable by sex
tab obesity_calc pulmonarycircdisorders, row //more in obese

tab age_decade diabetes3, row //all go up with age
tab age_decade hypertension, row
tab age_decade peripheralvascdisorders, row
tab age_decade renalfailure, row
*/

// -- outcome 
/*
//death
tab age_decade death, row //lot more death in olds; not linear
tab male death, row chi2 //no diff
tab obesity_calc death, row chi2 //no diff 

tab3way death age_decade enlargedpa, colpct //in middle age decades, pa enlarged folks die more
tab3way death age_decade enlargedratio, colpct //in middle age decades, paratio high folks die more.. not true <40. Lot of 30-40 normal PA:AA deaths
tab3way death age_decade pa_confusion_matrix, colpct //INTERESTING - very few with PA:AA enlarged who are young died

//lastfollowupyear - //duration of follow-up
bysort male: sum lastfollowupyear, detail //females have longer follow-up
bysort obesity_calc: sum lastfollowupyear, detail //obese have longer follow-up
bysort age_decade: sum lastfollowupyear //younger have longer follow-up
bysort pa_confusion_matrix: sum lastfollowupyear //enlargedpa have shorter follow-up; pa:aa high = no diff
bysort death: sum lastfollowupyear, detail //people died at an avg of 4.2 +/- 3.46 years (wide spread)

hist lastfollowupyear if death == 0
hist lastfollowupyear if death == 1

hist lastfollowupyear if death == 1, by(enlargedpa) xtitle("Time of death, by PA diameter status (years)")
graph export "Results and Figures/$S_DATE/Hist Time of Death PAd status.png", as(png) name("Graph") replace
hist lastfollowupyear if death == 1, by(enlargedratio) xtitle("Time of death, by PA:AA ratio status (years)")
graph export "Results and Figures/$S_DATE/Hist Time of Death PAAA status.png", as(png) name("Graph") replace

//mean time of death by group
bysort enlargedpa: summ lastfollowupyear if death == 1, detail //PAd enlarged: 3.53 [1.26 - 6.20] years
bysort enlargedratio: summ lastfollowupyear if death == 1, detail //Enlarged Ratio 3.28 [1.04 - 6.15] years
summ lastfollowupyear if death == 1, detail //Total: 3.55 [1.15 - 6.57] years

*/ 

// Misc counts within cohort
total PE_Infarction Iatrogenic_PE_Infarct Septic_PE Other_PE_or_Infarct Chronic_PulmHeartDz_OrPPH Chronic_PE Chronic_PulmHeartDz_Unspec ChronicPulmDz_NOS 
tab icdcodepresent
tab pulmonarycircdisorders

/* ---------------------------------
//Exploratory analysis to see what correlates with risk of only have PAd enlargement or only having PA:AA
--------------------------------*/
/*
mlogit pa_confusion_matrix male obesity_calc ib3.age_decade, rrr 
/* Multinomial logistic regression shows: 
Obesity: increases risk of PAd only, both
Sex: not independently associated with anything
Age: old -> mostly indepenently associated with likelihood of just PAd enlarged
	young -> strong association with just PA:AA 
	v old -> increased risk of both
*/ 

mlogit pa_confusion_matrix male obesity_calc ib3.age_decade hypertension diabetes3 chf pulmdisease pulmonarycircdisorders peripheralvascdisorders renalfailure, rrr 
/* Multinomial logistic regression shows:  (with everything)
CHF -> increased likelihood of just PAd enlarged
Periph Vasc -> increased likelihood of just PA:AA
CHF and Peripheral vasc borderline increase likelihood of both
After controlling for everything else.
*/ 

//Exploratory analysis to see what correlates with risk of only have PAd enlargement or only having AAd enlargement
mlogit pa_aa_confusion_matrix male obesity_calc ib3.age_decade, rrr 
/* Multinomial logistic regression shows: 
Obesity: increases risk of PA and AA, but PAd more. 
Sex: increased likelihood of AA enlargment
Age: old -> increases range of each independently, and both
*/ 

mlogit pa_aa_confusion_matrix male obesity_calc ib3.age_decade hypertension diabetes3 chf pulmdisease pulmonarycircdisorders peripheralvascdisorders renalfailure, rrr 
/* Multinomial logistic regression shows:  (with everything)
gender less strongly associated after controlling for the rest; age and obesity findings similar
CHF -> strongly associated with both enlarged but not either in isolation
HTN, DM, pulm_circ -> no indep assoc
Pulm -> low risk of AA only.
Periph Vasc, Renal -> increased likelihood of just PA (?weird)
*/ 


mlogit pa_confusion_matrix male obesity_calc ib3.age_decade hypertension diabetes3 chf pulmdisease pulmonarycircdisorders peripheralvascdisorders renalfailure, rrr 
/* Multinomial logistic regression shows:  (with everything)
CHF -> increased likelihood of just PAd enlarged
Periph Vasc -> increased likelihood of just PA:AA
CHF and Peripheral vasc borderline increase likelihood of both
After controlling for everything else.
*/ 

//Exploratory analysis to see what correlates with risk of only have PAd enlargement or only having AAd enlargement
mlogit pa_aa_confusion_matrix male obesity_calc ib3.age_decade, rrr 
/* Multinomial logistic regression shows: 
Obesity: increases risk of PA and AA, but PAd more. 
Sex: increased likelihood of AA enlargment
Age: old -> increases range of each independently, and both
*/ 

mlogit pa_aa_confusion_matrix male obesity_calc ib3.age_decade hypertension diabetes3 chf pulmdisease pulmonarycircdisorders peripheralvascdisorders renalfailure, rrr 
/* Multinomial logistic regression shows:  (with everything)
gender less strongly associated after controlling for the rest; age and obesity findings similar
CHF -> strongly associated with both enlarged but not either in isolation
HTN, DM, pulm_circ -> no indep assoc
Pulm -> low risk of AA only.
Periph Vasc, Renal -> increased likelihood of just PA (?weird)

*/ 


mlogit pa_confusion_matrix male obesity_calc ib3.age_decade hypertension diabetes3 chf pulmdisease pulmonarycircdisorders peripheralvascdisorders renalfailure, rrr 
/* Multinomial logistic regression shows:  (with everything)
CHF -> increased likelihood of just PAd enlarged
Periph Vasc -> increased likelihood of just PA:AA
CHF and Peripheral vasc borderline increase likelihood of both
After controlling for everything else.
*/ 

//Exploratory analysis to see what correlates with risk of only have PAd enlargement or only having AAd enlargement
mlogit pa_aa_confusion_matrix male obesity_calc ib3.age_decade, rrr 
/* Multinomial logistic regression shows: 
Obesity: increases risk of PA and AA, but PAd more. 
Sex: increased likelihood of AA enlargment
Age: old -> increases range of each independently, and both
*/ 

mlogit pa_aa_confusion_matrix male obesity_calc ib3.age_decade hypertension diabetes3 chf pulmdisease pulmonarycircdisorders peripheralvascdisorders renalfailure, rrr 
/* Multinomial logistic regression shows:  (with everything)
gender less strongly associated after controlling for the rest; age and obesity findings similar
CHF -> strongly associated with both enlarged but not either in isolation
HTN, DM, pulm_circ -> no indep assoc
Pulm -> low risk of AA only.
Periph Vasc, Renal -> increased likelihood of just PA (?weird)

*/ 


mlogit pa_confusion_matrix male obesity_calc ib3.age_decade hypertension diabetes3 chf pulmdisease pulmonarycircdisorders peripheralvascdisorders renalfailure, rrr 
/* Multinomial logistic regression shows:  (with everything)
CHF -> increased likelihood of just PAd enlarged
Periph Vasc -> increased likelihood of just PA:AA
CHF and Peripheral vasc borderline increase likelihood of both
After controlling for everything else.
*/ 

//Exploratory analysis to see what correlates with risk of only have PAd enlargement or only having AAd enlargement
mlogit pa_aa_confusion_matrix male obesity_calc ib3.age_decade, rrr 
/* Multinomial logistic regression shows: 
Obesity: increases risk of PA and AA, but PAd more. 
Sex: increased likelihood of AA enlargment
Age: old -> increases range of each independently, and both
*/ 

mlogit pa_aa_confusion_matrix male obesity_calc ib3.age_decade hypertension diabetes3 chf pulmdisease pulmonarycircdisorders peripheralvascdisorders renalfailure, rrr 
/* Multinomial logistic regression shows:  (with everything)
gender less strongly associated after controlling for the rest; age and obesity findings similar
CHF -> strongly associated with both enlarged but not either in isolation
HTN, DM, pulm_circ -> no indep assoc
Pulm -> low risk of AA only.
Periph Vasc, Renal -> increased likelihood of just PA (?weird)

*/ 
*/

//TODO: association with TRJet or other echo parameters would be an interesting sanity check. 
*** some data that adding RV parameter improves discriminatory ability: 
//Outcome: ED use or subsequent admission outcomes interesting.


//*ORIGINAL ANALYSIS FROM HERE*
/// Entire cohort /// 
/*
*** vs primary outcome *** 
ttest age, by(enlargedpa) // there is a statistically significant difference in patients with enlarged measurements, age 
ttest age, by(enlargedratio)
sum mpad, detail // gives percentiles for whole cohort
sum mpad if male==0, detail
sum mpad if male==1, detail
sum mpad if lastfollowupyear!=0, detail
tab enlargedpa enlargedratio //
sum mpad if enlargedpa==0
sum mpad if enlargedpa==1
sum mpaaa if enlargedratio==1 & male==0
sum mpaaa if enlargedratio==1 & male==1
// if want to do by gender,can add & male==0, &male==1
mean mpaaa if enlargedratio==0
mean mpaaa if enlargedratio==1

*=============* table 1 analysis for everybody, even those without any follow up 
ttest age, by(enlargedpa)
tab enlargedpa male, expect
tab enlargedpa male, row chi2
tab enlargedpa pulmdisease, expect // all expected >5, appropriate for chi square analysis 
tab enlargedpa pulmdisease, row chi2 // dependent variable first ... just not sure how we should present this
tab enlargedpa chf, expect
tab enlargedpa chf, row chi2 
tab enlargedpa diabetes3, expect
tab enlargedpa diabetes3, row chi2
tab enlargedpa hypertension, expect
tab enlargedpa hypertension, row chi2
tab enlargedpa obesity, expect
tab enlargedpa obesity, row chi2
tab enlargedpa pulmonarycircdisorders, expect
tab enlargedpa pulmonarycircdisorders, row chi2
tab enlargedpa peripheralvascdisorders, expect
tab enlargedpa peripheralvascdisorders, row chi2
tab enlargedpa renalfailure, expect
tab enlargedpa renalfailure, row chi2

//Brian addition: easier way to do this
//ssc install table1_mc, replace // if you haven't

*========================*
* table 1 (patient characteristics) on people with any follow-up 
*=========================*
tab age 
sum age,detail
sum age // without as much detail 
tab lastfollowupyear 
sum lastfollowupyear, detail
sum lastfollowupyear if lastfollowupyear!=0, detail // if there was follow-up, shows us the average last time there was follow up  
tab firstfollowupyear
sum firstfollowupyear, detail
sum firstfollowupyear if firstfollowupyear!=0, detail // removes those with 0

tab enlargedpa pulmdisease, expect // all expected >5, appropriate for chi square analysis 
tab enlargedpa pulmdisease, row chi2 // dependent variable first ... just not sure how we should present this
tab enlargedpa chf, expect
tab enlargedpa chf, row chi2 
tab enlargedpa diabetes3, expect
tab enlargedpa diabetes3, row chi2
tab enlargedpa hypertension, expect
tab enlargedpa hypertension, row chi2
tab enlargedpa obesity, expect
tab enlargedpa obesity, row chi2
tab enlargedpa pulmonarycircdisorders, expect
tab enlargedpa pulmonarycircdisorders, row chi2
tab enlargedpa peripheralvascdisorders, expect
tab enlargedpa peripheralvascdisorders, row chi2
// ??????How do we deterine which variables impact size of pulmonary artery

// if balanced between the groups do we even really have to adjust for them?

tab enlargedratio male, expect
tab enlargedratio male, row chi2
tab enlargedratio pulmdisease, expect // all expected >5, appropriate for chi square analysis 
tab enlargedratio pulmdisease, row chi2 // dependent variable first ... just not sure how we should present this
tab enlargedratio chf, expect
tab enlargedratio chf, row chi2 
tab enlargedratio diabetes3, expect
tab enlargedratio diabetes3, row chi2
tab enlargedratio hypertension, expect
tab enlargedratio hypertension, row chi2
tab enlargedratio obesity, expect
tab enlargedratio obesity, row chi2
tab enlargedratio pulmonarycircdisorders, expect
tab enlargedratio pulmonarycircdisorders, row chi2
tab enlargedratio peripheralvascdisorders, expect
tab enlargedratio peripheralvascdisorders, row chi2


*========================*
// analysis
*=========================*
tab death enlargedratio, expect // all expected >5, appropriate for chi square analysis 
tab death enlargedratio, col chi2
tab death enlargedpa, expect 
tab death enlargedpa, col chi2 

*/ 

/* DROP THOSE WITH NEW FOLLOW-UP AND NO COVARIATE INFORMATION*/ 

// only want those with follow-up for cox 
drop if lastfollowupday==0
count // count = 912, 78 dropped 
drop if missing(pulmdisease) //it's the same 23 that are missing all of the covariates.
save finalwithfollowupandcovariates_noempi, replace


/* 
Part of Table 1 for 2023 CHEST / Letter
*/ 
//consider age as a decade or group

table1_mc, by (age_bi_decade) ///
vars( ///
mpad conts %4.1f \ ///
enlargedpa cat %4.1f \ /// 
mpaaa conts %4.2f \ ///
enlargedratio cat %4.1f \ /// 
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol total(before) ///
saving("Results and Figures/$S_DATE/PA enlargement by Age.xlsx", replace)

//mpad_cat cat %4.1f \ ///
//mpaaa_cat cat %4.1f \ ///


table1_mc, by(enlargedpa) ///
vars( ///
age conts %4.0f \ /// 
male bin %4.1f \ /// 
obesity_calc bin %4.1f \ /// 
pulmdisease bin %4.1f \ /// 
chf bin %4.1f \ /// 
diabetes3 cat %4.1f \ /// 
hypertension bin %4.1f \ ///
pulmonarycircdisorders bin %4.1f \ /// 
peripheralvascdisorders bin %4.1f \ /// 
renalfailure bin %4.1f \ /// 
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol total(before) ///
saving("Results and Figures/$S_DATE/Table 1a - main demographics PAd.xlsx", replace)

table1_mc, by(enlargedratio) ///
vars( ///
age conts %4.0f \ /// 
male bin %4.1f \ /// 
obesity_calc bin %4.1f \ /// 
pulmdisease bin %4.1f \ /// 
chf bin %4.1f \ /// 
diabetes3 cat %4.1f \ /// 
hypertension bin %4.1f \ ///
pulmonarycircdisorders bin %4.1f \ /// 
peripheralvascdisorders bin %4.1f \ /// 
renalfailure bin %4.1f \ /// 
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol total(before) ///
saving("Results and Figures/$S_DATE/Table 1b - main demographics PAAA.xlsx", replace)
//Took out: 


table1_mc, by(death) ///
vars( ///
age conts %4.0f \ /// 
age_bi_decade cat %4.0f \ /// 
male bin %4.1f \ /// 
bmi conts %4.1f \ ///
pulmdisease bin %4.1f \ /// 
chf bin %4.1f \ /// 
diabetes3 cat %4.1f \ /// 
hypertension bin %4.1f \ ///
mpad conts %4.1f \ ///
enlargedpa cat %4.1f \ /// 
mpaaa conts %4.2f \ ///
enlargedratio cat %4.1f \ /// 
anyemergency bin %4.1f \ /// 
anyadmission bin %4.1f \ /// 
lastfollowupyear conts %4.1f \ /// 
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol total(before) ///
saving("Results and Figures/$S_DATE/Table 3 - main demographics by death.xlsx", replace)
//Took out: 
//obesity bin %4.1f \ /// 
//pulmonarycircdisorders bin %4.1f \ /// 
//peripheralvascdisorders bin %4.1f \ /// 
//renalfailure bin %4.1f \ /// 

/* 
Exploratory analyses: tables of univariate associations with the confusion matrices
*/ 
/*
table1_mc, by(pa_confusion_matrix) ///
vars( ///
age conts %4.0f \ /// 
age_decade cat %4.0f \ /// 
male bin %4.1f \ /// 
obesity_calc bin %4.1f \ /// 
bmi contn %4.1f \ /// 
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
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol total(before) ///
saving("Results and Figures/$S_DATE/All demographcis by PAAA and PAd.xlsx", replace)

table1_mc, by(pa_aa_confusion_matrix) ///
vars( ///
age conts %4.0f \ /// 
age_decade cat %4.0f \ /// 
male bin %4.1f \ /// 
obesity_calc bin %4.1f \ /// 
bmi contn %4.1f \ /// 
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
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol total(before) ///
saving("Results and Figures/$S_DATE/All demographcis by AAd and PAd.xlsx", replace)


table1_mc, by(enlargedaa) ///
vars( ///
age conts %4.0f \ /// 
age_decade cat %4.0f \ /// 
male bin %4.1f \ /// 
obesity_calc bin %4.1f \ /// 
bmi contn %4.1f \ /// 
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
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol total(before) ///
saving("Results and Figures/$S_DATE/All demographcis by AA_d.xlsx", replace)

*/ 



*=======================*
//Tables for CHEST Abstract / Letter
*=======================*
//TODO: coding a separate male mpad and female mpad variable would be better than the sex-norm
//TODO: I think the corrected IRR's might be better. 

//TODO: consider whether these are more convincing by major age group or decade?


table1_mc, by(enlargedpa) ///
vars( ///
age_bi_decade cat %4.0f \ /// 
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol ///
saving("Results and Figures/$S_DATE/CHEST-Age by PA_d simple.xlsx", replace)

table1_mc, by(enlargedratio) ///
vars( ///
age_bi_decade cat %4.0f \ /// 
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol ///
saving("Results and Figures/$S_DATE/CHEST-Age by PA_AA simple.xlsx", replace)

/* Analysis using the categories from Truong and PAd made up */ 

/*
table1_mc, by(age_decade) ///
vars( ///
sex_norm_mpad conts %4.1f \ ///
mpaaa conts %4.2f \ ///
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol ///
saving("Results and Figures/$S_DATE/CHEST-PAd and PAAA by decade.xlsx", replace)

table1_mc, by(mpad_cat) ///
vars( ///
age_decade cat %4.0f \ /// 
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol ///
saving("Results and Figures/$S_DATE/CHEST-Age by PA_d.xlsx", replace)

table1_mc, by(mpaaa_cat) ///
vars( ///
age_decade cat %4.0f \ /// 
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol ///
saving("Results and Figures/$S_DATE/CHEST-Age by PA_AA.xlsx", replace)

table1_mc, by(mpad_cat) ///
vars( ///
age conts %4.0f \ /// 
age_bi_decade cat %4.0f \ /// 
male bin %4.1f \ /// 
obesity bin %4.1f \ /// 
bmi contn %4.1f \ /// 
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
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol total(before) ///
saving("Results and Figures/$S_DATE/CHEST-All demographcis by PA_d.xlsx", replace)

table1_mc, by(mpaaa_cat) ///
vars( ///
age conts %4.0f \ /// 
age_bi_decade cat %4.0f \ /// 
male bin %4.1f \ /// 
obesity bin %4.1f \ /// 
bmi contn %4.1f \ /// 
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
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol total(before) ///
saving("Results and Figures/$S_DATE/CHEST-All demographcis by PA_AA.xlsx", replace)


tab age_bi_decade mpaaa_cat, row
tab age_bi_decade mpad_cat, row
tab mpaaa_cat mpad_cat

bysort enlargedpa: sum age, detail 
bysort enlargedratio: sum age, detail 


bysort mpad_cat: sum age, detail //even moreso at the extremes
bysort mpaaa_cat: sum age, detail 

summarize age, detail //overall median age 51

//Age to split the group?
bysort deceased: summarize age, detail
//alive: median 39.5
//deceased: median 69 --> half of events occur in above half in below
//unknown: 50
*/



/* 
REGRESSIONS EVALUATING WHICH CHARACTERSTICS WERE INDEPENDENTLY ASSOCIATED WITH LARGE PA
*/ 
// [ ] these IRR tables would be good to include in the Quarto doc. 

// PREVALENCE ESTIMATES

//Subdivided by decade
poisson enlargedpa ib3.age_decade i.male i.obesity_calc i.hypertension ib0.diabetes3 i.chf i.pulmdisease i.pulmonarycircdisorders i.peripheralvascdisorders i.renalfailure, irr robust
estimates store irr_pa_decade

poisson enlargedratio ib3.age_decade i.male i.obesity_calc i.hypertension ib0.diabetes3 i.chf i.pulmdisease i.pulmonarycircdisorders i.peripheralvascdisorders i.renalfailure, irr robust
estimates store irr_ratio_decade

//I like this one better for poster
coefplot irr_pa_decade, bylabel("PA diameter") subtitle(, size(large)) || ///
 irr_ratio_decade, bylabel("PA:AA Ratio") subtitle(, size(large)) ||, ///
 keep(0.age_decade 1.age_decade 2.age_decade 3.age_decade 4.age_decade 5.age_decade) ///
 legend(off) /// 
 eform /// 
 xscale(log) /// 
 xline(1) ///
 xlabel(0.5 1 2 4, labsize(large)) ///
 xscale(extend) ///
 xtitle("Age-, Sex-, BMI-, Comorbidity-Adjusted Relative Risk of Enlargement" "(PAd over 25 mm in women, 27mm in men. PA:AA over 0.9)" , size(medlarge)) ///
 yscale(extend) ///
 ylabel(, labsize(large)) ///
 ciopts(recast(rcap) lwidth(thick)) ///
 mlabel(string(@b,"%9.2f") + " [" + string(@ll,"%9.2f") + " - " + string(@ul,"%9.2f") + "]") /// 
 mlabsize(medlarge) ///
 mlabposition(12) /// 
 mlabgap(*1) /// 
 baselevels /// 
 scheme(white_tableau)
graph export "Results and Figures/$S_DATE/Sensitivity Regression Effect of IPW.png", as(png) name("Graph") replace



//Subdivided by 15 year category 
poisson enlargedpa ib1.age_bi_decade i.male i.obesity_calc i.hypertension ib0.diabetes3 i.chf i.pulmdisease i.pulmonarycircdisorders i.peripheralvascdisorders i.renalfailure, irr robust
estimates store irr_pa_bidecade

poisson enlargedratio ib1.age_bi_decade i.male i.obesity_calc i.hypertension ib0.diabetes3 i.chf i.pulmdisease i.pulmonarycircdisorders i.peripheralvascdisorders i.renalfailure, irr robust
estimates store irr_ratio_bidecade

coefplot irr_pa_bidecade, bylabel("PA diameter > 25 mm (W), > 27mm (M)") || irr_ratio_bidecade, bylabel("PA:AA Ratio > 0.9") ||, keep(0.age_bi_decade 1.age_bi_decade 2.age_bi_decade) legend(pos(2) ring(0) size(2.5)) eform xscale(log) xline(1) xscale(extend) xlabel(0.75 1 1.5 2)  xtitle("Relative Risk of Enlargement, controlling for gender, BMI, comorbidities" , size(small)) yscale(extend) headings( 0.age_bi_decade = "{bf:Age}")ciopts(recast(rcap) lwidth(thick)) mlabel(string(@b,"%9.2f") + " [ " + string(@ll,"%9.2f") + " - " + string(@ul,"%9.2f") + " ] " + cond(@pval<.001, "***", cond(@pval<.01, "**", cond(@pval<.05, "*", "")))) mlabsize(medsmall) mlabposition(12) mlabgap(*1) baselevels scheme(white_tableau)

//AA 
poisson enlargedaa ib3.age_decade male obesity_calc hypertension diabetes1 diabetes2 chf pulmdisease pulmonarycircdisorders peripheralvascdisorders renalfailure, irr robust

//[ ] todo: Splines 
/*
poisson enlargedpa age male obesity_calc hypertension diabetes1 diabetes2 chf pulmdisease pulmonarycircdisorders peripheralvascdisorders renalfailure, irr robust
outreg2 using "Results and Figures/$S_DATE/IRRs by Ratio and Diam", excel replace eform dec(2) pdec(3) stat(coef ci pval) title(Independent Associations with Enlarged PA and PAAA) label 

poisson enlargedratio age male obesity_calc hypertension diabetes1 diabetes2 chf pulmdisease pulmonarycircdisorders peripheralvascdisorders renalfailure, irr robust
outreg2 using "Results and Figures/$S_DATE/IRRs by Ratio and Diam", excel append eform dec(2) pdec(3) stat(coef ci pval) title(Independent Associations with Enlarged PA and PAAA) label 
*/

//LINEAR REGRESSION WITH SIZE
//PAd
regress mpad ib3.age_decade i.male i.obesity_calc i.hypertension ib0.diabetes3 i.chf i.pulmdisease i.pulmonarycircdisorders i.peripheralvascdisorders i.renalfailure

regress mpad ib1.age_bi_decade i.male i.obesity_calc i.hypertension ib0.diabetes3 i.chf i.pulmdisease i.pulmonarycircdisorders i.peripheralvascdisorders i.renalfailure
estimates store mpad_bi_decade

//Ratio
regress mpaaa ib3.age_decade i.male i.obesity_calc i.hypertension ib0.diabetes3 i.chf i.pulmdisease i.pulmonarycircdisorders i.peripheralvascdisorders i.renalfailure

regress mpaaa ib1.age_bi_decade i.male i.obesity_calc i.hypertension ib0.diabetes3 i.chf i.pulmdisease i.pulmonarycircdisorders i.peripheralvascdisorders i.renalfailure
estimates store mpaaa_bi_decade

//Aorta
regress ascendingaorta ib3.age_decade i.male i.obesity_calc i.hypertension ib0.diabetes3 i.chf i.pulmdisease i.pulmonarycircdisorders i.peripheralvascdisorders i.renalfailure

regress ascendingaorta ib1.age_bi_decade i.male i.obesity_calc i.hypertension ib0.diabetes3 i.chf i.pulmdisease i.pulmonarycircdisorders i.peripheralvascdisorders i.renalfailure
estimates store aa_bi_decade

//output regression coefficients to a table.
outreg2 [mpad_bi_decade mpaaa_bi_decade aa_bi_decade ] using "Results and Figures/$S_DATE/MPAD MPAAA AA regressions", word replace keep(0.age_bi_decade 1.age_bi_decade 2.age_bi_decade) dec(3) stat(coef ci pval) title(Change in Size by Age) label 


/*-----------
Kaplan Meier Curves by threshold.
Note: if I could figure out how to get the decade adjustment working properly it would be nice to include some of these. 
--------------*/ 
stset lastfollowupyear, failure(death==1)
/*
//PA d
sts test enlargedpa, logrank 
//survival
sts graph, by(enlargedpa) risktable(0(2.5)10) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) risktable(, order(1 "Normal PAd" 2 "Enlarged PAd"))
graph export "Results and Figures/$S_DATE/Unadjusted KM PAd.png", as(png) name("Graph") replace

//hazard
sts graph, by(enlargedpa) hazard xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) 
graph export "Results and Figures/$S_DATE/Unadjusted Hazard by PAd.png", as(png) name("Graph") replace

//adjusted hazard
sts graph, by(enlargedpa) adjustfor(male age) hazard xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) 
graph export "Results and Figures/$S_DATE/Adjusted Age Sex Hazard by PAd.png", as(png) name("Graph") replace

//without comorbidities
sts test enlargedpa if (known_assoc_comorb == 0), logrank
sts graph if (known_assoc_comorb == 0), by(enlargedpa) risktable(0(2.5)10) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) risktable(, order(1 "Normal PAd" 2 "Enlarged PAd")) title ("K-M Survival Estimate in Patients without Pulm or Pulm Vasc Dz")
graph export "Results and Figures/$S_DATE/No Known Comorb Unadjusted KM PAd.png", as(png) name("Graph") replace

//with comorbidities
sts test enlargedpa if (known_assoc_comorb == 1), logrank
sts graph if (known_assoc_comorb == 1), by(enlargedpa) risktable(0(2.5)10) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) risktable(, order(1 "Normal PAd" 2 "Enlarged PAd")) title ("K-M Survival Estimate in Patients with Pulm and Pulm Vasc Dz")
graph export "Results and Figures/$S_DATE/Known Comorb Unadjusted KM PAd.png", as(png) name("Graph") replace

//PA;AA
sts test enlargedratio, logrank
sts graph, by(enlargedratio) risktable(0(2.5)10) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) risktable(, order(1 "Normal PA:AA" 2 "Increased PA:AA"))
graph export "Results and Figures/$S_DATE/Unadjusted KM PAAA.png", as(png) name("Graph") replace

//hazard
sts graph, by(enlargedratio) hazard xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) 
graph export "Results and Figures/$S_DATE/Unadjusted Hazard by PAAA.png", as(png) name("Graph") replace

//Adjusted hazard
sts graph, by(enlargedratio) adjustfor(male age) hazard xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) 
graph export "Results and Figures/$S_DATE/Adjusted Age Sex Hazard by PAAA.png", as(png) name("Graph") replace

//subgroup without comorbidities - note, swapped curves; age is confounding
sts test enlargedratio if (known_assoc_comorb == 0), logrank
sts graph if (known_assoc_comorb == 0), by(enlargedratio) risktable(0(2.5)10) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) risktable(, order(1 "Normal PA:AA" 2 "Increased PA:AA")) title ("K-M Survival Estimate in Patients without Pulm or Pulm Vasc Dz")
graph export "Results and Figures/$S_DATE/No Known Comorb Unadjusted KM PAAA.png", as(png) name("Graph") replace

//with comorbidities
sts test enlargedratio if (known_assoc_comorb == 1), logrank
sts graph if (known_assoc_comorb == 1), by(enlargedratio) risktable(0(2.5)10) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) risktable(, order(1 "Normal PA:AA" 2 "Increased PA:AA")) title ("K-M Survival Estimate in Patients with Pulm and Pulm Vasc Dz")
graph export "Results and Figures/$S_DATE/Known Comorb Unadjusted KM PAAA.png", as(png) name("Graph") replace

*/



*========================*
// Hazard ratio per each disease
*=========================*
/*
bysort enlargedpa: stcox pulmdisease age0
bysort pulmdisease: stcox enlargedpa age // 
stcox enlargedpa age0 if pulmdisease==1 //  HR 1.4 if adjusting forage
stcox enlargedpa if pulmdisease==1 // HR 1.8
stcox enlargedpa if chf==1 // 1.46 (0.99-2.14)
stcox enlargedpa if hypertension==1
stcox enlargedpa if obesity==1
stcox enlargedpa if peripheralvascdisorders==1
stcox enlargedpa if pulmonarycircdisorders==1 // 2.45
stcox enlargedpa if renalfailure==1 // 1.16 (0.61-2.19)
stcox enlargedpa if diabetes3>0 // 1.44
stcox enlargedratio if pulmdisease==1 // 1.3 - CI 0.95-1.78
stcox enlargedratio if chf==1 // 1.78 (1.2-2.64)
stcox enlargedratio if hypertension==1 // 1.35
stcox enlargedratio if diabetes3>0
stcox enlargedratio if peripheralvascdisorders==1
stcox enlargedratio if renalfailure==1
stcox enlargedratio if pulmonarycircdisorders==1

// multivariable regression
stcox enlargedpa age male obesity // limited model
stcox enlargedpa age0 male0 peripheralvascdisorders0 obesity0 pulmdisease0 diabetes10 diabetes20 renalfailure0 pulmonarycircdisorders0 chf0

// multivariable regression for ratio
stcox enlargedratio age male obesity // limited model
stcox enlargedratio age0 male0 peripheralvascdisorders0 obesity0 pulmdisease0 diabetes10 diabetes20 renalfailure0 pulmonarycircdisorders0 chf0
*/  


// DATA EXPLORATION, REGRESSIONS, AND FIGURE FROM INITIAL PAPER BMC CIRC
/*
// Unadjusted enlargedratio
stcox enlargedratio
estimates store unadjusted_ratio_full
stcox enlargedratio if (known_assoc_comorb == 1)
estimates store unadjusted_ratio_comorb
stcox enlargedratio if (known_assoc_comorb == 0)
estimates store unadjusted_ratio_no_comorb

// model 1 enlargedratio
stcox enlargedratio ib3.age_decade male
stcurve, survival at(enlargedratio = (0 1)) xlabel(0(2.5)12.5) ylabel(0(.10)1) note("Adjusted by Age, Sex") xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(ring(0) position(2) rows(2) order(1 "PA:AA normal" 2 "PA:AA increased")) title("Survival by PA:AA increased vs not")
graph export "Results and Figures/$S_DATE/AgeSex KM PAAA.png", as(png) name("Graph") replace

//subgroups
stphplot, by(enlargedratio) //this honestly doesn't look that bad?
stcoxkm, by(enlargedratio) pred1opts(symbol(none) lpattern(dash)) obs2opts(symbol(none) color(gray)) pred2opts(symbol(none) color(gray) lpattern(dash)) legend(ring(0) position(2) col(2)) //this honestly doesn't look all that bad
estat phtest, detail //really just the decades that violate proportional hazards
//The reson this violates proportional hazards is because low age individuals have almost no effect (few die) while mid range do (of people that die, more have these comorbidities), then eldely have less again
//this may be an issue of competing risks: how many people in the <50 age group die cardiovascular deaths? 
//-- still a notable finding, that many patients under 50 have nominally elevated PA:AA, but do not face elevated 10y mortality; unclear if this is because the timeframe si not long enough for cardiovascular related death, or the true normal range changes in younger folks. 

stcox enlargedratio male, strata(age_decade) // this is the right way to not assume proportional hazards, but can't graph it
estimates store ratio_model1_full
estat phtest, detail 

stcox enlargedratio male if (known_assoc_comorb == 1), strata(age_decade)
estimates store ratio_model1_comorb
stcox enlargedratio male if (known_assoc_comorb == 0), strata(age_decade)
estimates store ratio_model1_no_comorb

stcox enlargedratio##known_assoc_comorb male, strata(age_decade) // p = 0.076 for interaction; ugh

//model2 enlargedratio
stcox enlargedratio ib3.age_decade male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf
stcurve, survival at(enlargedratio = (0 1)) xlabel(0(2.5)12.5) ylabel(0(.10)1) note("Adjusted by Age, Sex, and Comorbidities") xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(ring(0) position(2) rows(2) order(1 "PA:AA normal" 2 "PA:AA increased")) title("Survival by PA:AA increased vs not")
graph export "Results and Figures/$S_DATE/AgeSexComorb KM PAAA.png", as(png) name("Graph") replace
estat phtest, detail //really just the decades that violate proportional hazards

stcox enlargedratio male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf, strata(age_decade) //used
estimates store ratio_model2_full
stphplot, by(enlargedratio) ylabel(, format(%9.1f)) //  if the lines are parallel and Cox is OK 
//Subgroups
stcox enlargedratio male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf if (known_assoc_comorb == 0), strata(age_decade)
estimates store ratio_model2_no_comorb
stcox enlargedratio male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf if (known_assoc_comorb == 1), strata(age_decade)
estimates store ratio_model2_comorb

stcox enlargedratio##known_assoc_comorb male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf, strata(age_decade) //p = 0.073 for interaction model 2

stcox i.mpaaa_cat##known_assoc_comorb male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf, strata(age_decade) 
testparm i.mpaaa_cat##known_assoc_comorb // p 0.029

//chf pulm pulm vasc individiually
stcox enlargedratio if (chf == 1)
estimates store ratio_unadjusted_chf
stcox enlargedratio if (pulmdisease == 1)
estimates store ratio_unadjusted_pulm
stcox enlargedratio if (pulmonarycircdisorders == 1)
estimates store ratio_unadjusted_pvasc
stcox enlargedratio male if (chf == 1), strata(age_decade)
estimates store ratio_model_1_chf
stcox enlargedratio male if (pulmdisease == 1), strata(age_decade)
estimates store ratio_model_1_pulm
stcox enlargedratio male if (pulmonarycircdisorders == 1), strata(age_decade)
estimates store ratio_model_1_pvasc
stcox enlargedratio male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf if (chf == 1), strata(age_decade)
estimates store ratio_model_2_chf
stcox enlargedratio male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf if (pulmdisease == 1), strata(age_decade)
estimates store ratio_model_2_pulm
stcox enlargedratio male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf if (pulmonarycircdisorders == 1), strata(age_decade)
estimates store ratio_model_2_pvasc

//overlap between comorbidities - [/] probably useful to include in the quarto  - meh
pvenn2 pulmdisease pulmonarycircdisorders chf, plabel("Pulmonary_Disease" "Pulmonary_Vascular_Disease" "Congestive_Heart_Failure") drawtotal(0) n(3000)
graph export "Results and Figures/$S_DATE/ComorbidityVennDiagram.png", as(png) name("Graph") replace

//output regression coefficients to a table.
outreg2 [unadjusted_ratio_full unadjusted_ratio_comorb unadjusted_ratio_no_comorb ratio_model1_full ratio_model1_comorb ratio_model1_no_comorb ratio_model2_full ratio_model2_comorb ratio_model2_no_comorb ] using "Results and Figures/$S_DATE/ratio_HRs_by_model", excel replace keep(enlargedratio) eform dec(3) stat(coef ci pval) title(HR of Mortality) label 

//regression coefficients for the individual comorbidities
outreg2 [ratio_unadjusted_chf ratio_unadjusted_pulm ratio_unadjusted_pvasc ratio_model_1_chf ratio_model_1_pulm ratio_model_1_pvasc ratio_model_2_chf ratio_model_2_pulm ratio_model_2_pvasc] using "Results and Figures/$S_DATE/indiv_comorb_ratio_HRs_by_model", excel replace keep(enlargedratio) eform dec(3) stat(coef ci pval) title(HR of Mortality) label 

//unadjustedpa
stcox enlargedpa
estimates store unadjusted_pa_full
stcox enlargedpa if (known_assoc_comorb == 1)
estimates store unadjusted_pa_comorb
stcox enlargedpa if (known_assoc_comorb == 0)
estimates store unadjusted_pa_no_comorb

// model1 enlargedpa
stcox enlargedpa ib3.age_decade male
stcurve, survival at(enlargedpa = (0 1)) xlabel(0(2.5)12.5) ylabel(0(.10)1) note("Adjusted by Age, Sex") xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(ring(0) position(2) rows(2) order(1 "PA_d normal" 2 "PA_d > ULN")) title("Survival by PA_D > ULN vs not")
estat phtest, detail //really just the decades that violate proportional hazards
stphplot, by(enlargedpa) 

stcox enlargedpa male, strata(age_decade)
estimates store pa_model1_full
estat phtest, detail 
stcox enlargedpa male if (known_assoc_comorb == 1), strata(age_decade)
estimates store pa_model1_comorb
stcox enlargedpa male if (known_assoc_comorb == 0), strata(age_decade)
estimates store pa_model1_no_comorb
stcox enlargedpa##known_assoc_comorb male, strata(age_decade) // no interaction

// model2 enlargedpa
stcox enlargedpa ib3.age_decade male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf
estat phtest, detail //really just the decades that violate proportional hazards
stcurve, survival at(enlargedpa = (0 1)) xlabel(0(2.5)12.5) ylabel(0(.10)1) note("Adjusted by Age, Sex, and Comorbidities") xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(ring(0) position(2) rows(2) order(1 "PA_d normal" 2 "PA_d > ULN")) title("Survival by PA_D > ULN vs not")

stcox enlargedpa male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf, strata(age_decade)
estimates store pa_model2_full
stcox enlargedpa male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf if (known_assoc_comorb == 0), strata(age_decade)
estimates store pa_model2_no_comorb
stcox enlargedpa male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf if (known_assoc_comorb == 1), strata(age_decade)
estimates store pa_model2_comorb
stcox enlargedpa##known_assoc_comorb male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf, strata(age_decade) //no interaction p = 0.23

//Individual comorbidity ratios
stcox enlargedpa if (chf == 1)
estimates store pa_unadjusted_chf
stcox enlargedpa if (pulmdisease == 1)
estimates store pa_unadjusted_pulm
stcox enlargedpa if (pulmonarycircdisorders == 1)
estimates store pa_unadjusted_pvasc
stcox enlargedpa male if (chf == 1), strata(age_decade)
estimates store pa_model_1_chf
stcox enlargedpa male if (pulmdisease == 1), strata(age_decade)
estimates store pa_model_1_pulm
stcox enlargedpa male if (pulmonarycircdisorders == 1), strata(age_decade)
estimates store pa_model_1_pvasc
stcox enlargedpa male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf if (chf == 1), strata(age_decade)
estimates store pa_model_2_chf
stcox enlargedpa male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf if (pulmdisease == 1), strata(age_decade)
estimates store pa_model_2_pulm
stcox enlargedpa male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf if (pulmonarycircdisorders == 1), strata(age_decade)
estimates store pa_model_2_pvasc


//output regression coefficients to a table.
outreg2 [unadjusted_pa_full unadjusted_pa_comorb unadjusted_pa_no_comorb pa_model1_full pa_model1_comorb pa_model1_no_comorb pa_model2_full pa_model2_comorb pa_model2_no_comorb] using "Results and Figures/$S_DATE/pa_HRs_by_model", excel replace keep(enlargedpa) eform dec(3) stat(coef ci pval) title(HR of Mortality) label 

outreg2 [pa_unadjusted_chf pa_unadjusted_pulm pa_unadjusted_pvasc pa_model_1_chf pa_model_1_pulm pa_model_1_pvasc pa_model_2_chf pa_model_2_pulm pa_model_2_pvasc] using "Results and Figures/$S_DATE/indiv_comorb_pa_HRs_by_model", excel replace keep(enlargedpa) eform dec(3) stat(coef ci pval) title(HR of Mortality) label 


// Visualizations: "Forest Plot" by subgroup
/* WITH SUBGROUPS */ 
coefplot (unadjusted_ratio_full, keep(enlargedratio) label("Unadjusted, Full Sample") offset(0.45)) (ratio_model1_full, keep(enlargedratio) label("Model 1, Full Sample") offset(0.4)) (ratio_model2_full, keep(enlargedratio) label("Model 2, Full Sample") offset(0.35)) (unadjusted_ratio_comorb, keep(enlargedratio) label("Unadjusted, Subgroup +Pulm/CHF/PH") offset(0.05)) (ratio_model1_comorb, keep(enlargedratio) label("Model 1, Subgroup +Pulm/CHF/PH") offset(0)) (ratio_model2_comorb, keep(enlargedratio) label("Model 2, Subgroup +Pulm/PH") offset(-0.05)) (unadjusted_ratio_no_comorb, keep(enlargedratio) label("Unadjusted, Subgroup -Pulm/CHF/PH") offset(-0.35)) (ratio_model1_no_comorb, keep(enlargedratio) label("Model 1, Subgroup -Pulm/CHF/PH") offset(-0.4)) (ratio_model2_no_comorb, keep(enlargedratio) label("Model 2, Subgroup -Pulm/CHF/PH") offset(-0.45)) || , drop(_cons) eform xscale(log) xline(1) baselevels xlabel(0.5 1 2 4) xscale(extend) xtitle("Adjusted Hazard Ratio of Mortality" , size(small)) ciopts(recast(rcap)) ylab(, labs(vsmall)) xsize(7) ysize(7) note("Unadjusted" "Model 1: Adjusted for Age, Sex" "Model 2: Adjusted for Age, Sex, Comorbidities") mlabel format(%9.2g) mlabposition(12) mlabgap(*1) legend(position(6) rows(3) size(vsmall))

coefplot (unadjusted_pa_full, keep(enlargedpa) label("Unadjusted, Full Sample") offset(0.45)) (pa_model1_full, keep(enlargedpa) label("Model 1, Full Sample") offset(0.4)) (pa_model2_full, keep(enlargedpa) label("Model 2, Full Sample") offset(0.35)) (unadjusted_pa_comorb, keep(enlargedpa) label("Unadjusted, Subgroup +Pulm/CHF/PH") offset(0.05)) (pa_model1_comorb, keep(enlargedpa) label("Model 1, Subgroup +Pulm/CHF/PH") offset(0)) (pa_model2_comorb, keep(enlargedpa) label("Model 2, Subgroup +Pulm/CHF/PH") offset(-0.05)) (unadjusted_pa_no_comorb, keep(enlargedpa) label("Unadjusted, Subgroup -Pulm/CHF/PH") offset(-0.35)) (pa_model1_no_comorb, keep(enlargedpa) label("Model 1, Subgroup -Pulm/PH") offset(-0.4)) (pa_model2_no_comorb, keep(enlargedpa) label("Model 2, Subgroup -Pulm/CHF/PH") offset(-0.45)) || , drop(_cons) eform xscale(log) xline(1) baselevels xlabel(0.5 1 2 4) xscale(extend) xtitle("Adjusted Hazard Ratio of Mortality" , size(small)) ciopts(recast(rcap)) ylab(, labs(vsmall)) xsize(7) ysize(7) note("Unadjusted" "Model 1: Adjusted for Age, Sex" "Model 2: Adjusted for Age, Sex, Comorbidities") mlabel format(%9.2g) mlabposition(12) mlabgap(*1) legend(position(6) rows(3) size(vsmall))

//combined
coefplot (unadjusted_ratio_full, keep(enlargedratio) label("Unadjusted, Full Sample") offset(0.35)) (ratio_model1_full, keep(enlargedratio) label("Model 1, Full Sample") offset(0.3)) (ratio_model2_full, keep(enlargedratio) label("Model 2, Full Sample") offset(0.25)) (unadjusted_ratio_comorb, keep(enlargedratio) label("Unadjusted, Subgroup +Pulm/CHF/PH") offset(0.05)) (ratio_model1_comorb, keep(enlargedratio) label("Model 1, Subgroup +Pulm/CHF/PH") offset(0)) (ratio_model2_comorb, keep(enlargedratio) label("Model 2, Subgroup +Pulm/CHF/PH") offset(-0.05)) (unadjusted_ratio_no_comorb, keep(enlargedratio) label("Unadjusted, Subgroup -Pulm/CHF/PH") offset(-0.35)) (ratio_model1_no_comorb, keep(enlargedratio) label("Model 1, Subgroup -Pulm/CHF/PH") offset(-0.4)) (ratio_model2_no_comorb, keep(enlargedratio) label("Model 2, Subgroup -Pulm/CHF/PH") offset(-0.45)), bylabel("PA:AA > 0.9") || (unadjusted_pa_full, keep(enlargedpa) label("Unadjusted, Full Sample") offset(0.35)) (pa_model1_full, keep(enlargedpa) label("Model 1, Full Sample") offset(0.3)) (pa_model2_full, keep(enlargedpa) label("Model 2, Full Sample") offset(0.25)) (unadjusted_pa_comorb, keep(enlargedpa) label("Unadjusted, Subgroup +Pulm/CHF/PH") offset(0.05)) (pa_model1_comorb, keep(enlargedpa) label("Model 1, Subgroup +Pulm/CHF/PH") offset(0)) (pa_model2_comorb, keep(enlargedpa) label("Model 2, Subgroup +Pulm/CHF/PH") offset(-0.05)) (unadjusted_pa_no_comorb, keep(enlargedpa) label("Unadjusted, Subgroup -Pulm/CHF/PH") offset(-0.35)) (pa_model1_no_comorb, keep(enlargedpa) label("Model 1, Subgroup -Pulm/CHF/PH") offset(-0.4)) (pa_model2_no_comorb, keep(enlargedpa) label("Model 2, Subgroup -Pulm/CHF/PH") offset(-0.45)), bylabel("PAd > ULN") || , drop(_cons) eform xscale(log) xline(1) baselevels xlabel(0.5 1 2 4) xscale(extend) xtitle("Hazard Ratio of Mortality" , size(small)) ciopts(recast(rcap)) ylab(, labs(vsmall)) xsize(7) ysize(7) mlabel format(%9.2g) mlabposition(12) mlabgap(*1) legend(size(vsmall) position(6) rows(3) note("Unadjusted" "Model 1: Adjusted for Age (Decade), Sex" "Model 2: Adjusted for Age (Decade), Sex, Obesity, DM, CKD, PVD, CHF, PVasc, Pulm Dz")) swapnames 
graph export "Results and Figures/$S_DATE/Combined HR Unadjusted Model1 Model2.png", as(png) name("Graph") replace

/* FINAL FIGURE FOR MANUSCRIPT: without subgroups */
coefplot (unadjusted_pa_full, keep(enlargedpa) label("Unadjusted")) (pa_model1_full, keep(enlargedpa) label("Model 1")) (pa_model2_full, keep(enlargedpa) label("Model 2")), bylabel("Enlarged Pulmonary Artery Diameter") || (unadjusted_ratio_full, keep(enlargedratio) label("Unadjusted")) (ratio_model1_full, keep(enlargedratio) label("Model 1 adjusts for age (by decade) and sex")) (ratio_model2_full, keep(enlargedratio) label("Model 2 adjusts for age (by decade), sex, and comorbidities")), bylabel("Increased Pulmonary Artery : Ascending Aorta Ratio") || , eform swapnames xscale(log) xline(1) xscale(extend) xtitle("Hazard Ratio of Mortality" , size(medsmall)) ciopts(recast(rcap) lwidth(thick)) mlabel(string(@b,"%9.2f") + " [ " + string(@ll,"%9.2f") + " - " + string(@ul,"%9.2f") + " ] " + cond(@pval<.001, "***", cond(@pval<.01, "**", cond(@pval<.05, "*", "")))) mlabsize(medsmall) mlabposition(12) mlabgap(*1) legend(size(medsmall) position(6) rows(3) note("* = p < 0.05; ** = p < 0.01; *** = p < 0.001", size(medsmall)))  scheme(cleanplots) ylab("", labs(medsmall)) msymbol(S)
graph export "Results and Figures/$S_DATE/No SG - Combined HR Unadjusted Model1 Model2.png", as(png) name("Graph") replace
//Comorbodities: Obesity, Diabetes, Renal Failure, Peripheral Vascular Disease, Heart Failure, Pulmonary Vascular Disease, Pulmonary Disease



//todo: add category labels and synchronize colors, remove legend. 

/*
coefplot (ratio_model1_full, keep(enlargedratio) label("Model 1, Full Sample") offset(0.4)) (ratio_model2_full, keep(enlargedratio) label("Model 2, Full Sample") offset(0.3)) (ratio_model1_comorb, keep(enlargedratio) label("Model 1, Subgroup +Pulm/CHF/PH") offset(0.1)) (ratio_model2_comorb, keep(enlargedratio) label("Model 2, Subgroup +Pulm/CHF/PH") offset(0)) (ratio_model1_no_comorb, keep(enlargedratio) label("Model 1, Subgroup -Pulm/CHF/PH") offset(-0.2)) (ratio_model2_no_comorb, keep(enlargedratio) label("Model 2, Subgroup -Pulm/CHF/PH") offset(-0.30)) || , drop(_cons) eform xscale(log) xline(1) baselevels xlabel(0.5 1 2 4) xscale(extend) xtitle("Adjusted Hazard Ratio of Mortality" , size(small)) ciopts(recast(rcap)) ylab(, labs(vsmall)) xsize(7) ysize(7) note("Model 1: Adjusted for Age, Sex" "Model 2: Adjusted for Age, Sex, Comorbidities") mlabel format(%9.2g) mlabposition(12) mlabgap(*1) legend(position(6) rows(3))

coefplot (pa_model1_full, keep(enlargedpa) label("Model 1, Full Sample") offset(0.4)) (pa_model2_full, keep(enlargedpa) label("Model 2, Full Sample") offset(0.3)) (pa_model1_comorb, keep(enlargedpa) label("Model 1, Subgroup +Pulm/CHF/PH") offset(0.1)) (pa_model2_comorb, keep(enlargedpa) label("Model 2, Subgroup +Pulm/CHF/PH") offset(0)) (pa_model1_no_comorb, keep(enlargedpa) label("Model 1, Subgroup -Pulm/CHF/PH") offset(-0.2)) (pa_model2_no_comorb, keep(enlargedpa) label("Model 2, Subgroup -Pulm/CHF/PH") offset(-0.30)), drop(_cons) eform xscale(log) xline(1) baselevels xlabel(0.5 1 2 4) xscale(extend) xtitle("Adjusted Hazard Ratio of Mortality" , size(small)) ciopts(recast(rcap)) ylab(, labs(vsmall)) xsize(7) ysize(7) note("Model 1: Adjusted for Age, Sex" "Model 2: Adjusted for Age, Sex, Comorbidities") mlabel format(%9.2g) mlabposition(12) mlabgap(*1) legend(position(6) rows(3))

//combined
coefplot (ratio_model1_full, keep(enlargedratio) label("Model 1, Full Sample") offset(0.4)) (ratio_model2_full, keep(enlargedratio) label("Model 2, Full Sample") offset(0.3)) (ratio_model1_comorb, keep(enlargedratio) label("Model 1, Subgroup +Pulm/CHF/PH") offset(0.1)) (ratio_model2_comorb, keep(enlargedratio) label("Model 2, Subgroup +Pulm/CHF/PH") offset(0)) (ratio_model1_no_comorb, keep(enlargedratio) label("Model 1, Subgroup -Pulm/CHF/PH") offset(-0.2)) (ratio_model2_no_comorb, keep(enlargedratio) label("Model 2, Subgroup -Pulm/CHF/PH") offset(-0.30)), bylabel("PA:AA > 0.9") || (pa_model1_full, keep(enlargedpa) label("Model 1, Full Sample") offset(0.4)) (pa_model2_full, keep(enlargedpa) label("Model 2, Full Sample") offset(0.3)) (pa_model1_comorb, keep(enlargedpa) label("Model 1, Subgroup +Pulm/CHF/PH") offset(0.1)) (pa_model2_comorb, keep(enlargedpa) label("Model 2, Subgroup +Pulm/CHF/PH") offset(0)) (pa_model1_no_comorb, keep(enlargedpa) label("Model 1, Subgroup -Pulm/CHF/PH") offset(-0.2)) (pa_model2_no_comorb, keep(enlargedpa) label("Model 2, Subgroup -Pulm/CHF/PH") offset(-0.30)), bylabel("PAd > ULN") || , drop(_cons) eform xscale(log) xline(1) baselevels xlabel(0.5 1 2 4) xscale(extend) xtitle("Adjusted Hazard Ratio of Mortality" , size(small)) ciopts(recast(rcap)) ylab(, labs(vsmall)) xsize(7) ysize(7) mlabel format(%9.2g) mlabposition(12) mlabgap(*1) legend(position(6) rows(3) note("Model 1: Adjusted for Age (Decade), Sex" "Model 2: Adjusted for Age (Decade), Sex, Obesity, DM, CKD, Peripheral Vascular Disease")) swapnames 
*/ 


*/


//Nonlinear effect of age using the unstratified proportional hazards regression

//NOTE THE FOLLOWING CODE CHUNK DOES NOT TEST WHAT ITS PURPORTED TO TEST
/*
stcox enlargedratio ib3.age_decade male
estimates store ratio_age_unstrat
stcox enlarged pa ib3.age_decade male
estimates store pa_age_unstrat

//Perhaps worth including this in the quarto doc
coefplot ratio_age_unstrat, drop(male enlargedratio) bylabel("HR By Decade (PA:AA Model 1)") || pa_age_unstrat, drop(male enlargedpa) bylabel("HR By Decade (PAd Model 1)")  ||, drop(_cons ) eform xscale(log) xline(1) baselevels xlabel(0.03 0.06 0.125 0.25 0.5 1 2 4) xscale(extend) xtitle("Adjusted Hazard Ratio of Mortality" , size(small)) ciopts(recast(rcap)) ylab(, labs(vsmall)) xsize(7) ysize(7)
graph export "Results and Figures/$S_DATE/HR of Age by Decade.png", as(png) name("Graph") replace
*/

  
/* ------------------
ANALYSIS OF THE PA CONFUSION AND AA CONFUSION MATRIX SURVIVAL RATES
--------------------*/ 
  
//TODO: figure out how to get the K-M graph adjustment to work
//Note: the best version of this graph would adjust for age_dace and male - but we need to mean center to accoont for the command's default of holding covariates at 0 by default - see stoddard ch 5-23 p 43

//Looking at all 4 possibilities...  of none, PAd only, PA:AA only, both
sts test pa_confusion_matrix, logrank
//this is unadjusted by age so I don't think all that useful.
sts graph, by(pa_confusion_matrix) xlabel(0(2.5)12.5) ci xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(3) rows(8)) 

//Adjustment doesn't seem to be working properly
sts graph, by(pa_confusion_matrix) adjustfor(age_decade male, atmeans) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(3) rows(8)) 

sts graph, by(pa_confusion_matrix) hazard xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(3) rows(8)) 

sts graph, by(pa_confusion_matrix) adjustfor(age_decade male, atmeans) hazard xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(3) rows(8)) 
stphplot, by(pa_confusion_matrix)

//this is adjusting for age_decade and sex [ ] include this in the quarto document. 
stcox ib0.pa_confusion_matrix ib3.age_decade male
stcurve, survival at(pa_confusion_matrix = (0 1 2 3)) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) note("Adjusted by Age (Decade) and Sex") legend(ring(0) position(2) rows(2) order(1 "Normal" 2 "PA_d > ULN only" 3 "PA:AA > 0.9 only" 4 "Both")) title("Survival by PA_D > ULN, PA:AA > 0.9, Neither, Or Both")
graph export "Results and Figures/$S_DATE/PA Confusion Matrix KM age male.png", as(png) name("Graph") replace

//fully adjusted version
stcox ib0.pa_confusion_matrix ib3.age_decade male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf
stcurve, survival at(pa_confusion_matrix = (0 1 2 3)) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") note("Adjusted by Age (Decade), Sex, and Comorbidities") legend(ring(0) position(2) rows(2) order(1 "Normal" 2 "PA_d > ULN only" 3 "PA:AA > 0.9 only" 4 "Both")) title("Survival by PA_D > ULN, PA:AA > 0.9, Neither, Or Both")
graph export "Results and Figures/$S_DATE/PA Confusion Matrix KM age male comorbs.png", as(png) name("Graph") replace

//Looking at all 4 possibilities...  of none, PAd only, AAd only, both
sts test pa_aa_confusion_matrix, logrank
//this is unadjusted by age so I don't think all that useful.
sts graph, by(pa_aa_confusion_matrix) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") title("Survival by PAd > ULN, AAd > ULN, Neither, Or Both (All Ages)") risktable(0(2.5)10, size(vsmall)) legend(position(6) rows(2)) legend(ring(0) position(5) rows(2) order(1 "Normal" 2 "PA_d > ULN only" 3 "AAd > ULN only" 4 "Both"))
graph export "Results and Figures/$S_DATE/PA AA Confusion Matrix KM.png", as(png) name("Graph") replace
//this is adjusting for age_decade and sex. [ ] include this one in the quarto document. 


stcox ib0.pa_aa_confusion_matrix ib3.age_decade male
stcurve, survival at(pa_aa_confusion_matrix = (0 1 2 3)) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) note("Adjusted by Age (Decade) and Sex") legend(ring(0) position(2) rows(2) order(1 "Normal" 2 "PA_d > ULN only" 3 "AAd > ULN only" 4 "Both")) title("Survival by PA_D > ULN, AA_d > ULN, Neither, Or Both")
graph export "Results and Figures/$S_DATE/PA AA Confusion Matrix KM age male.png", as(png) name("Graph") replace

//fully adjusted version
stcox ib0.pa_aa_confusion_matrix ib3.age_decade male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf
stcurve, survival at(pa_aa_confusion_matrix = (0 1 2 3)) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") note("Adjusted by Age (Decade), Sex, and Comorbidities") legend(ring(0) position(2) rows(2) order(1 "Neither" 2 "PA_d > ULN only" 3 "AA_d > ULN only" 4 "Both")) title("Survival by PA_D > ULN, AA_d > ULN, Neither, Or Both")
graph export "Results and Figures/$S_DATE/PA AA Confusion Matrix KM age male comorbs.png", as(png) name("Graph") replace

//order(1 "Normal" 2 "PA_d > ULN only" 3 "PA:AA > 0.9 only" 4 "Both")
//label define pa_aa_confusion_lab 0 "Neither" 1 "Only Enlarged PAd" 2 "Only High AAd" 3 "Both"


*=======================*
//Regressions for CHEST Abstract
*=======================*

//Categorical Exposure HR's -- do the cutoffs make sense across age ranges?

/* -------------------
ANALYSIS BY CATEGORY OF ENLARGEMENT 
-------------------*/ 

//PA_d - Analysis

//Adjusted by age (and implicitly by sex)
stcox ib1.mpad_cat##ib1.age_bi_decade //mpad_cat already accounts for sex.
testparm mpad_cat#age_bi_decade //no interaction, p = 0.8; not expecting much excitement here
sts graph, by(mpad_cat) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) risktable(0(2.5)10, size(small)) title("Survival (all ages)")
graph export "Results and Figures/$S_DATE/PAd all ages KM.png", as(png) name("Graph") replace
sts graph if age_bi_decade == 0, by(mpad_cat) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) risktable(0(2.5)10, size(small)) title("Survival in <50 year olds")
graph export "Results and Figures/$S_DATE/PAd age below 50 KM.png", as(png) name("Graph") replace
sts graph if age_bi_decade == 1, by(mpad_cat) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) risktable(0(2.5)10, size(small)) title("Survival in 50-65 year olds")
graph export "Results and Figures/$S_DATE/PAd age 50 75 KM.png", as(png) name("Graph") replace
sts graph if age_bi_decade == 2, by(mpad_cat) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) risktable(0(2.5)10, size(small)) title("Survival in 65+ year olds")
graph export "Results and Figures/$S_DATE/PAd age above 65 KM.png", as(png) name("Graph") replace

stcox ib1.mpad_cat, strata(age_bi_decade) //Note: MPAD cat already accounts for sex. [ although on further thought, it doesn't necessarily take care of any adjustments in our cohort that aren't explained by the 2mm difference...]
estimates store pasize_in_all
stcox ib1.mpad_cat if age_bi_decade == 0
estimates store pasize_in_young
stcox ib1.mpad_cat peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure if age_bi_decade == 0
stcox ib1.mpad_cat if age_bi_decade == 1
estimates store pasize_in_med
stcox ib1.mpad_cat peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure if age_bi_decade == 1
stcox ib1.mpad_cat if age_bi_decade == 2
estimates store pasize_in_old
stcox ib1.mpad_cat peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure if age_bi_decade == 2

//[ ] include this in the quarto document 
coefplot (pasize_in_all, drop(male)), bylabel("All Ages (stratified)") || (pasize_in_young, drop(male)), bylabel("Age < 50") || (pasize_in_med, drop(male)), bylabel("Age 50-65") || (pasize_in_old, drop(male)), bylabel("Age 65+") ||, drop(_cons) eform xscale(log) xline(1) baselevels xlabel(0.25 0.5 1 2 4) xscale(extend) xtitle("Adjusted Hazard Ratio of Mortality" , size(medium)) ciopts(recast(rcap)) ylab(, labs(medium)) xsize(5) ysize(10) ytitle("PA Diameter", size(medlarge)) mlabel format(%9.2g) mlabposition(12) mlabgap(*1) byopts(compact cols(1))
graph export "Results and Figures/$S_DATE/PAd HR at each age.png", as(png) name("Graph") replace


//PA:AA Analysis by severity category
stcox ib1.mpaaa_cat##ib1.age_bi_decade male
testparm mpaaa_cat#age_bi_decade

sts graph, by(mpaaa_cat) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) risktable(0(2.5)10, size(small)) title("Survival (all ages)")
graph export "Results and Figures/$S_DATE/PAAA All Ages.png", as(png) name("Graph") replace
sts graph if age_bi_decade == 0, by(mpaaa_cat) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) risktable(0(2.5)10, size(small)) title("Survival in <50 year olds")
graph export "Results and Figures/$S_DATE/PAAA age below 50 KM.png", as(png) name("Graph") replace
sts graph if age_bi_decade == 1, by(mpaaa_cat) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) risktable(0(2.5)10, size(small)) title("Survival in <50-65 year olds")
graph export "Results and Figures/$S_DATE/PAAA age 50 65 KM.png", as(png) name("Graph") replace
sts graph if age_bi_decade == 2, by(mpaaa_cat) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) risktable(0(2.5)10, size(small)) title("Survival in 65+ year olds")
graph export "Results and Figures/$S_DATE/PAAA age above 65 KM.png", as(png) name("Graph") replace

stcox ib1.mpaaa_cat male, strata(age_bi_decade)
estimates store ratio_in_all
stcox ib1.mpaaa_cat male if age_bi_decade == 0
estimates store ratio_in_young
stcox ib1.mpaaa_cat male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure if age_bi_decade == 0
stcox ib1.mpaaa_cat male if age_bi_decade == 1
estimates store ratio_in_med
stcox ib1.mpaaa_cat male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure if age_bi_decade == 1
stcox ib1.mpaaa_cat male if age_bi_decade == 2
estimates store ratio_in_old
stcox ib1.mpaaa_cat male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure if age_bi_decade == 2

//include this in the quarto analysis
coefplot (ratio_in_all, drop(male)), bylabel("All Ages (stratified)") || (ratio_in_young, drop(male)), bylabel("Age < 50") || (ratio_in_med, drop(male)), bylabel("Age 50-65") || (ratio_in_old, drop(male)), bylabel("Age 65+") ||, drop(_cons) eform xscale(log) xline(1) baselevels xlabel(0.25 0.5 1 2 4) xscale(extend) xtitle("Adjusted Hazard Ratio of Mortality" , size(medium)) ciopts(recast(rcap)) ylab(, labs(medium)) xsize(5) ysize(10) ytitle("PA:AA ratio", size(medlarge)) mlabel format(%9.2g) mlabposition(12) mlabgap(*1) byopts(compact cols(1))
graph export "Results and Figures/$S_DATE/PAAA HR at each age.png", as(png) name("Graph") replace

//AORTA SIZE - Analysis - 
stcox ib1.aa_cat##ib1.age_bi_decade //aa_cat already accounts for sex.
testparm aa_cat#age_bi_decade //no interaction, p = 0.67; not expecting much excitement here
sts graph if age_bi_decade == 0, by(aa_cat) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) risktable(0(2.5)10, size(small)) title("Survival in <50 year olds")
graph export "Results and Figures/$S_DATE/AA age below 50 KM.png", as(png) name("Graph") replace
sts graph if age_bi_decade == 1, by(aa_cat) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) risktable(0(2.5)10, size(small)) title("Survival in 50-65 year olds")
graph export "Results and Figures/$S_DATE/AA age 50 65 KM.png", as(png) name("Graph") replace
sts graph if age_bi_decade == 2, by(aa_cat) xlabel(0(2.5)12.5) xtitle("Follow-up (years)") ytitle("Probability of Survival") legend(position(6) rows(2)) risktable(0(2.5)10, size(small)) title("Survival in 65+ year olds")
graph export "Results and Figures/$S_DATE/AA age above 65 KM.png", as(png) name("Graph") replace

stcox ib1.aa_cat, strata(age_bi_decade) //Note: AA cat already accounts for sex.
estimates store aasize_in_all
stcox ib1.aa_cat if age_bi_decade == 0
estimates store aasize_in_young
stcox ib1.aa_cat peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure if age_bi_decade == 0
stcox ib1.aa_cat if age_bi_decade == 1
estimates store aasize_in_med
stcox ib1.aa_cat peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure if age_bi_decade == 1
stcox ib1.aa_cat if age_bi_decade == 2
estimates store aasize_in_old
stcox ib1.aa_cat peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure if age_bi_decade == 2

coefplot (aasize_in_all, drop(male)), bylabel("All Ages (stratified)") || (aasize_in_young, drop(male)), bylabel("Age < 50") || (aasize_in_med, drop(male)), bylabel("Age 50-65") || (aasize_in_old, drop(male)), bylabel("Age 65+") ||, drop(_cons) eform xscale(log) xline(1) baselevels xlabel(0.25 0.5 1 2 4) xscale(extend) xtitle("Adjusted Hazard Ratio of Mortality" , size(medium)) ciopts(recast(rcap)) ylab(, labs(medium)) xsize(5) ysize(10) ytitle("AA Diameter", size(medlarge)) mlabel format(%9.2g) mlabposition(12) mlabgap(*1) byopts(compact cols(1))
graph export "Results and Figures/$S_DATE/AA HR at each age.png", as(png) name("Graph") replace


//EVALUATION INTERACTIONS 

// By categorical enlargement
stcox ib3.age_decade ib1.mpaaa_cat male
stcox ib1.mpaaa_cat male, strata(age_decade)  //except fit a spline instead of MPAAA_cat
	
stcox ib3.age_decade ib1.mpad_cat male
stcox ib1.mpad_cat male, strata(age_decade) //except fit a spline instead of MPAD_cat

//Is there an interaction by decade? You lose a lot of power by dropping the decade resolution, or increasing the categories of enlarged ratio
stcox ib3.age_decade##ib0.enlargedratio male
testparm i.age_decade#enlargedratio  // p 0.042

stcox ib1.age_bi_decade##ib0.enlargedratio male
testparm i.age_bi_decade#enlargedratio  // p 0.05

stcox ib3.age_decade##enlargedratio male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf //still present w full models? 
testparm i.age_decade#enlargedratio // p 0.0582

stcox ib3.age_decade##enlargedpa male
testparm i.age_decade#enlargedpa // p 0.11 for interaction with enlargedpa

stcox ib3.age_decade##enlargedpa male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf 
testparm i.age_decade#enlargedpa // p 0.07


/* ----------------------
FOLLOW-UP / DEATH TIMES
-----------------------*/

hist lastfollowupyear if death == 1, by(enlargedpa age_bi_decade) xtitle("Time of death, years") 
graph export "Results and Figures/$S_DATE/Hist Time of Death PAd status.png", as(png) name("Graph") replace
hist lastfollowupyear if death == 1, by(enlargedratio age_bi_decade) xtitle("Time of deaths (years)")
graph export "Results and Figures/$S_DATE/Hist Time of Death PAAA status.png", as(png) name("Graph") replace

hist lastfollowupyear if death == 0, by(enlargedpa age_bi_decade) xtitle("Time of last follow-up, years") 
graph export "Results and Figures/$S_DATE/Hist Time of last follow-up PAd status.png", as(png) name("Graph") replace
hist lastfollowupyear if death == 0, by(enlargedratio age_bi_decade) xtitle("Time of deaths (years)")
graph export "Results and Figures/$S_DATE/Hist Time of last follow-up PAAA status.png", as(png) name("Graph") replace


/* -----------------
DISTRIBUTRION OF VALUES
-------------------*/
//[ ] include this in the quarto analysis 
range paaa_range 0.6 1.2

twoway kdensity mpaaa if age_bi_decade==0, recast(area) fcolor(green%33) lcolor(green) lwidth(*0.75) || ///
kdensity mpaaa if age_bi_decade==1, recast(area) fcolor(midblue%33) lcolor(midblue) lwidth(*0.75) || ///
kdensity mpaaa if age_bi_decade==2, recast(area) fcolor(orange%33) lcolor(orange) lwidth(*0.75) ||, ///
legend(subtitle("Age", size(large)) order(1 "Age < 50" 2 "Age 50-65" 3 "Age 65+") ring(0) size(large)) ///
xtitle("Pulmonary Artery to Ascending Aorta Ratio", size(large)) ///
ytitle("Relative Frequency", size(large)) ///
xline(0.9) xlabel(0.6(0.1)1.2, labsize(large)) ///
ylabel("") ///
title("Distribution of PA:AA by Age Strata", size(vlarge)) ///
text(3 0.95 "Mean: {&darr} by 0.02 per decade" "Median: {&darr} by 0.02 per decade" "90th %tile: {&darr} by 0.01 per decade", placement(e) size(medlarge)) ///
scheme(white_tableau)
graph export "Results and Figures/$S_DATE/Frequency PA_AA by Age Group.png", as(png) name("Graph") replace

//[ ] inclue this in the quarto analysis 
twoway kdensity mpad if age_bi_decade==0, recast(area) fcolor(green%33) lcolor(green) lwidth(*0.75) || ///
kdensity mpad if age_bi_decade==1, recast(area) fcolor(midblue%33) lcolor(midblue) lwidth(*0.75) || ///
kdensity mpad if age_bi_decade==2, recast(area) fcolor(orange%33) lcolor(orange) lwidth(*0.75) ||, ///
legend(subtitle("Age", size(large)) order(1 "Age < 50" 2 "Age 50-65" 3 "Age 65+") ring(0) size(large)) ///
xtitle("Sex-Normalized PA diameter (mm)", size(large)) ///
ytitle("Relative Frequency", size(large)) ///
xline(28) xlabel(19(3)37, labsize(large)) ///
ylabel("") ///
title("Distribution of Pulmonary Artery Diameter by Age Strata", size(vlarge)) ///
text(0.085 29.5 "Mean: {&uarr} by 0.57mm per decade" "Median: {&uarr} by 0.64mm per decade" "90th %tile: {&uarr} by 1.0mm per decade", placement(e) size(medlarge)) ///
scheme(white_tableau)
graph export "Results and Figures/$S_DATE/Frequency PA_d by Age Group.png", as(png) name("Graph") replace

//Not sure if its worth combining these... but it actually seems like the labels might throw people off. 
//graph save "temp_mpaaa_old.gph", replace
//graph combine temp_mpaaa_old.gph temp_mpaaa_mid.gph temp_mpaaa_young.gph, ///
//xcommon col(1)

//--------------------------
/* Expected Ranges by Age */ 
//--------------------------

//todo: improve the colors --> abandon the schemes and just color the scatter
//todo: make the lines thicker.

twoway lfitci mpad age, stdf ciplot(rline) || scatter mpad age ||, xlabel(,labsize(medlarge)) ylabel(,labsize(medlarge)) xtitle("Age (years)", size(4)) ytitle("PA Size (mm)", size(4)) scheme(white_w3d) legend(off) title("Pulmonary Artery Diameter") xsize(5) ysize(5)
graph save "temp_mpad_age.gph", replace

twoway lfitci ascendingaorta age, stdf ciplot(rline) || scatter ascendingaorta age ||, xlabel(,labsize(medlarge)) ylabel(,labsize(medlarge)) xtitle("Age (years)", size(4)) ytitle("Ascending Aorta (mm)", size(4)) scheme(white_tableau) legend(off) title("Ascending Aorta Diameter") xsize(5) ysize(5)
graph save "temp_aad_age.gph", replace

twoway lfitci mpaaa age, stdf ciplot(rline) || scatter mpaaa age ||, xlabel(,labsize(medlarge)) ylabel(,labsize(medlarge)) xtitle("Age years", size(4)) ytitle("PA:AA Ratio", size(4)) scheme(white_jet) legend(off) title("PA to AA Ratio") xsize(5) ysize(5)
graph save "temp_ratio_age.gph", replace

graph combine temp_mpad_age.gph temp_aad_age.gph temp_ratio_age.gph, ///
 col(3) ///
 imargin(0 0 0 0) graphregion(margin(l=22 r=22)) ///
 note("95% Prediction Intervals Indicated by Grey Lines", pos(6)) ///
 xsize(9) ysize(3)
graph export "Results and Figures/$S_DATE/PA AA PAAA by Age continuous.png", as(png) name("Graph") replace
 // title("Variation in PA, AA, and PA:AA by Age") ///

/* -----------------
QUANTILE REGRESSION

Why? you can adjust for the mortality rate on the back-end, or swap it to the front end
------------------*/

regress mpaaa c.age 
sqreg mpaaa c.age, q(0.5 0.9)
estat coefplot

regress mpad c.age i.male
sqreg mpad c.age i.male, q(0.5 0.9)
estat coefplot
 
//TODO: could output these into tables as a potentially interesting finding. 
 
 
 
 
/* ----------
Generate Z-scores & percentiles after (linear) age adjustment
*/ 
//TODO; ensure these are accurately modeled linearly 
//TODO: any advantage to using heteroskedastic model? 
//TODO: what is the size of the std deviation? -> instead of the percentile, could do "___cm larger than age (sex)_median " --> reprogram this using residuals. 

regress mpaaa c.age
predict mpaaa_age_z, rstandard
hist mpaaa_age_z, normal
twoway scatter mpaaa_age_z age //no heteroskedasticity
gen mpaaa_age_ptile = normal(mpaaa_age_z)
list mpaaa_age_ptile mpaaa_age_z
hist mpaaa_age_ptile 

preserve
gen mpaaa_age_ptile_rounded = round(mpaaa_age_ptile, 0.01)
mkspline2 rc = mpaaa_age_ptile_rounded, cubic nknots(4) displayknots            
assert float(mpaaa_age_ptile_rounded) == float(rc1) // VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
stcox rc*
estat ic
estat concordance
estat phtest, detail
levelsof mpaaa_age_ptile_rounded, local(levels)
xblc rc*, covname(mpaaa_age_ptile_rounded) at(`r(levels)') reference(0.5) eform generate(pa or lb ub) // options for reference: minimum, mean. Could do most common? 
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)), yscale(log extend) scheme(cleanplots) legend(off) xlabel(0(0.1)1) xmtick(0(0.05)1) ylabel( 1 2 4, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("Pulmonary Artery to Ascending Aorta Ratio Percentile (Age-Adjusted)") title("Risk of Mortality by PA:AA") yline(1, lp("shortdash") lc(gs10)) xline(0.9, lp("shortdash_dot") lc(gs10)) note("Lines added at traditional threshold of abnormal (0.9) and HR = 1." "Hazard Ratio Referenced to Population Mean 0.77)" "Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAAA Splines - Whole Cohort Age-Z-score.png", as(png) name("Graph") replace
restore


regress mpad c.age i.male
predict mpad_age_z, rstandard
hist mpad_age_z, normal
twoway scatter mpad_age_z age //no heteroskedasticity
gen mpad_age_ptile = normal(mpad_age_z)
list mpad_age_ptile mpad_age_z
hist mpad_age_ptile   //interesting mpad has a bit more of a fat tail than the mpaaa
 
 
preserve
gen mpad_age_ptile_rounded = round(mpad_age_ptile, 0.01)
mkspline2 rc = mpad_age_ptile_rounded, cubic nknots(4) displayknots            
assert float(mpad_age_ptile_rounded) == float(rc1) // VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
stcox rc*
estat ic
estat concordance
estat phtest, detail
levelsof mpad_age_ptile_rounded, local(levels)
xblc rc*, covname(mpad_age_ptile_rounded) at(`r(levels)') reference(0.5) eform generate(pa or lb ub) // options for reference: minimum, mean. Could do most common? 
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)), yscale(log extend) scheme(cleanplots) legend(off) xlabel(0(0.1)1) xmtick(0(0.05)1) ylabel( 1 2 4, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("Pulmonary Artery Diameter Percentile (Age & Sex-Adjusted)") title("Risk of Mortality by PAd") yline(1, lp("shortdash") lc(gs10)) xline(0.9, lp("shortdash_dot") lc(gs10)) note("Lines added at traditional threshold of abnormal (0.9) and HR = 1." "Hazard Ratio Referenced to Population Mean 0.77)" "Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAd Splines - Whole Cohort Age-Z-score.png", as(png) name("Graph") replace

restore
 
//could also do this by decile categorical.... because I wonder if the ends are not being pulled around as much as you'd think?
 
 
 
 
/* ------------------------
RESTRICTED CUBIC SPLINES ANALYSIS---------
----------------------*/

//Splines


//TODO: HOW FRAGILE IS THIS ANALYSIS TO CHANGES?  - continuous age seems not to make much difference (postrcspline to check )
//todo: [x] reference as minimum hazard - does not seem particularly more informative. 
//TODO: maybe try population mean?  PA_d 25.1; PA:AA 0.77 -- SEEMS GOOD. Go with this. 
//TODO: try to add density plot   [separate]

//TODO: explore difference with unrestricted cubic splines and different spline points [ doesn't seem to make much difference... where does 4 come from? Harrell's textbook. 4 or 5 defensible in this case


// Harrell modeling strategies: (https://hbiostat.org/rmsc/genreg.html#sec-relax.linear) 
// tails behavior badly - hence restriction
// n dependent on fit vs n, suggestion: n<30 -> 3, n 30-100 -> 4, n 100+ -> 5
// location of knots doesn't matter much 3: .1,.5,.9 // 4: .05 .35 .65 . 95 // 5: 0.05 .275 .5 .725 .95 
// ---> can use AIC to choose ==> in this case, AIC lower for lower knots. 

/*-----------------------------
--------PA:AA ANALYSIS---------
-----------------------------*/
/*
------- WHOLE COHORT, UNADJUSTED
*/ 
/*
preserve
gen mpaaa_rounded = round(mpaaa, 0.01)
mkspline2 rc = mpaaa_rounded, cubic nknots(4) displayknots            
assert float(mpaaa_rounded) == float(rc1) // VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
stcox rc* // MODEL WITH FULL SPLINE
estat ic
levelsof mpaaa_rounded, local(levels)
xblc rc*, covname(mpaaa_rounded) at(`r(levels)') reference(0.77) eform generate(pa or lb ub) // options for reference: minimum, mean. Could do most common? 
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,0.65,1.3), yscale(log extend) scheme(cleanplots) legend(off) xlabel(0.65(0.05)1.3) xmtick(0.65(0.05)1.3) ylabel(0.5 1 2 4, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("PA:AA") title("HR for mortality; entire cohort; Unadjusted") yline(1, lp("shortdash") lc(gs10)) xline(0.9, lp("shortdash_dot") lc(gs10)) note("Healthy Population Mean (PA:AA = 0.77, Truong 2012) taken as reference." "Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (0.9) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAAA Splines - Whole Cohort Unadj.png", as(png) name("Graph") replace
restore
*/

/*
------- WHOLE COHORT, JUST ADJUSTED FOR SEX AND AGE
MAIN FIGURE FOR CHEST ABSTRACT
*/ 



preserve
gen mpaaa_rounded = round(mpaaa, 0.01)
mkspline2 rc = mpaaa_rounded, cubic nknots(4) displayknots            
assert float(mpaaa_rounded) == float(rc1) // VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
stcox rc* male, strata(age_bi_decade) // MODEL WITH FULL SPLINE
//stcox c.age##(c.rc?) male // MODEL WITH FULL SPINE AND INTERACTION TERM WITH AGE? - [ ] how to visualize this? 
estat ic
estat concordance
estat phtest, detail
levelsof mpaaa_rounded if inrange(mpaaa_rounded,0.65,1.2), local(levels)
xblc rc*, covname(mpaaa_rounded) at(`r(levels)') reference(0.77) eform generate(pa or lb ub) // options for reference: minimum, mean. Could do most common? 
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) ///
 (line or pa, sort lc(black) lp(l)) if inrange(pa,0.65,1.2), /// 
 yscale(log extend) /// 
 scheme(cleanplots) /// 
 legend(off) /// 
 xlabel(0.7(0.1)1.2, labsize(large)) ///
 xmtick(0.7(0.1)1.2) ///
 ylabel( 1 2 4, angle(horiz) format(%2.1fc) labsize(large)) /// 
 ytitle("Hazard Ratio of Mortality", size(large)) /// 
 xtitle("Pulmonary Artery to Ascending Aorta Ratio", size(large)) /// 
 title("Age-adjusted Mortality Risk by PA:AA", size(vlarge)) /// 
 yline(1, lp("shortdash") lc(gs10)) /// 
 xline(0.9, lp("shortdash_dot") lc(gs10)) ///
 text(0.85 1.06 "95% confidence interval (black dash)", size(medlarge)) ///
 text(2.5 .77 "Referenced to" "healthy population" "PA:AA mean (0.77)", size(medlarge))
graph export "Results and Figures/$S_DATE/HR PAAA Splines - Whole Cohort Adj AgeSex.png", as(png) name("Graph") replace
restore

// note("Lines added at traditional threshold of abnormal (0.9) and HR = 1." "Hazard Ratio Referenced to Population Mean 0.77)" "Dashed line = 95% Confidence interval")


//[ ] consider combining axes
/*
graph save "temp_mpaaa_old.gph", replace
graph combine temp_mpaaa_old.gph temp_mpaaa_mid.gph temp_mpaaa_young.gph, ///
xcommon col(1) //size?
*/ 



/*
------- WHOLE COHORT, ADJUSTED FOR SEX, AGE, AND COMORBIDITIES
*/ 
/*
preserve
gen mpaaa_rounded = round(mpaaa, 0.01)
mkspline2 rc = mpaaa_rounded, cubic nknots(4) displayknots            
assert float(mpaaa_rounded) == float(rc1)
stcox rc* male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf, strata(age_decade)
estat ic
levelsof mpaaa_rounded, local(levels)
xblc rc*, covname(mpaaa_rounded) at(`r(levels)') reference(0.77) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,0.65,1.3), yscale(log extend) legend(off) scheme(cleanplots) xlabel(0.65(0.05)1.3) xmtick(0.65(0.05)1.3) ylabel(1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("PA:AA") title("HR for mortality; entire cohort; adjusted for age(decade), sex, and comorbidities") yline(1, lp("shortdash") lc(gs10)) xline(0.9, lp("shortdash_dot") lc(gs10)) note("Healthy Population mean (PA:AA = 0.77, Truong 2012) taken as reference." "Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (0.9) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAAA Splines - Whole Cohort Adj AgeSexComorb.png", as(png) name("Graph") replace
restore
*/


//Hmmm, not sure the subgroups here. 
/*
------- Age < 50 , JUST ADJUSTED FOR SEX AND AGE
*/ 
preserve
//keep if age_bi_decade == 0
gen mpaaa_rounded = round(mpaaa, 0.01) if age_bi_decade == 0
mkspline2 rc = mpaaa_rounded, cubic nknots(3) displayknots            
// VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
assert float(mpaaa_rounded) == float(rc1) if age_bi_decade == 0
// MODEL WITH FULL SPLINE
stcox rc* male ib1.age_decade if age_bi_decade == 0 // c.age //strata(age_decade)
levelsof mpaaa_rounded, local(levels) 
xblc rc*, covname(mpaaa_rounded) at(`r(levels)') reference(0.77) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(forest_green forest_green) lp(longdash longdash)) ///
 (line or pa, sort lc(forest_green) lp(l)) if inrange(pa,0.65,1.3), /// 
 yscale(log extend) ///
 scheme(cleanplots) ///
 legend(off) ///
 xlabel(0.65(0.05)1.3) ///
 xmtick(0.65(0.05)1.3) ///
 ylabel(1 2 4 8, angle(horiz) format(%2.1fc)) ///
 ytitle("Hazard Ratio of Mortality") ///
 xtitle("PA:AA") ///
 title("HR for mortality; Age < 50; adjusted for age(decade) and sex") ///
 yline(1, lp("shortdash") lc(gs10)) ///
 xline(0.9, lp("shortdash_dot") lc(gs10)) ///
 note("Sample mean (PA:AA = 0.85) taken as reference. Restricted cubic splines with 3 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (0.9) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAAA Splines - Age less 50 Adj AgeSex.png", as(png) name("Graph") replace
restore

/*
------- Age < 50, ADJUSTED FOR SEX, AGE, AND COMORBIDITIES
*/ 
/*
preserve
keep if age_bi_decade == 0
gen mpaaa_rounded = round(mpaaa, 0.01)
mkspline2 rc = mpaaa_rounded, cubic nknots(4) displayknots            
// VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
assert float(mpaaa_rounded) == float(rc1)
// MODEL WITH FULL SPLINE
stcox rc* male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf, strata(age_decade)
levelsof mpaaa_rounded, local(levels)
xblc rc*, covname(mpaaa_rounded) at(`r(levels)') reference(0.77) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,0.65,1.3), yscale(log extend) legend(off) scheme(cleanplots) xlabel(0.65(0.05)1.3) xmtick(0.65(0.05)1.3) ylabel(1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("PA:AA") title("HR for mortality; Age < 50; adjusted for age(decade), sex, and comorbidities") yline(1, lp("shortdash") lc(gs10)) xline(0.9, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA:AA = 0.85) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (0.9) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAAA Splines - Age less 50 Adj AgeSexComorb.png", as(png) name("Graph") replace
restore
*/

/*
------- Age 50-65 , JUST ADJUSTED FOR SEX AND AGE
*/ 
preserve
keep if age_bi_decade == 1
gen mpaaa_rounded = round(mpaaa, 0.01)
mkspline2 rc = mpaaa_rounded, cubic nknots(4) displayknots            
// VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
assert float(mpaaa_rounded) == float(rc1)
// MODEL WITH FULL SPLINE
stcox rc* male, strata(age_decade)
levelsof mpaaa_rounded, local(levels)
xblc rc*, covname(mpaaa_rounded) at(`r(levels)') reference(0.77) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,0.65,1.3), yscale(log extend) scheme(cleanplots) legend(off) xlabel(0.65(0.05)1.3) xmtick(0.65(0.05)1.3) ylabel(1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("PA:AA") title("HR for mortality; Age 50-65; adjusted for age(decade) and sex") yline(1, lp("shortdash") lc(gs10)) xline(0.9, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA:AA = 0.85) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (0.9) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAAA Splines - Age 50-65 Adj AgeSex.png", as(png) name("Graph") replace
restore

/*
------- Age 50-65, ADJUSTED FOR SEX, AGE, AND COMORBIDITIES
*/ 
/*
preserve
keep if age_bi_decade == 1
gen mpaaa_rounded = round(mpaaa, 0.01)
mkspline2 rc = mpaaa_rounded, cubic nknots(4) displayknots            
// VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
assert float(mpaaa_rounded) == float(rc1)
// MODEL WITH FULL SPLINE
stcox rc* male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf, strata(age_decade)
levelsof mpaaa_rounded, local(levels)
xblc rc*, covname(mpaaa_rounded) at(`r(levels)') reference(0.77) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,0.65,1.3), yscale(log extend) legend(off) scheme(cleanplots) xlabel(0.65(0.05)1.3) xmtick(0.65(0.05)1.3) ylabel(1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("PA:AA") title("HR for mortality; Age 50-65; adjusted for age(decade), sex, and comorbidities") yline(1, lp("shortdash") lc(gs10)) xline(0.9, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA:AA = 0.85) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (0.9) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAAA Splines - Age 50-65 Adj AgeSexComorb.png", as(png) name("Graph") replace
restore
*/

/*
------- Age 65+, JUST ADJUSTED FOR SEX AND AGE
*/ 
preserve
keep if age_bi_decade == 2
gen mpaaa_rounded = round(mpaaa, 0.01)
mkspline2 rc = mpaaa_rounded, cubic nknots(4) displayknots            
assert float(mpaaa_rounded) == float(rc1)
stcox rc* male, strata(age_decade)
levelsof mpaaa_rounded, local(levels)
xblc rc*, covname(mpaaa_rounded) at(`r(levels)') reference(0.77) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,0.65,1.3), yscale(log extend) scheme(cleanplots) legend(off) xlabel(0.65(0.05)1.3) xmtick(0.65(0.05)1.3) ylabel(1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("PA:AA") title("HR for mortality; Age 65+; adjusted for age(decade) and sex") yline(1, lp("shortdash") lc(gs10)) xline(0.9, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA:AA = 0.85) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (0.9) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAAA Splines - Age 65plus Adj AgeSex.png", as(png) name("Graph") replace

restore

/*
------- Age 65+, ADJUSTED FOR SEX, AGE, AND COMORBIDITIES
*/ 
/*
preserve
keep if age_bi_decade == 2
gen mpaaa_rounded = round(mpaaa, 0.01)
mkspline2 rc = mpaaa_rounded, cubic nknots(4) displayknots            
// VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
assert float(mpaaa_rounded) == float(rc1)
// MODEL WITH FULL SPLINE
stcox rc* male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf, strata(age_decade)
levelsof mpaaa_rounded, local(levels)
xblc rc*, covname(mpaaa_rounded) at(`r(levels)') reference(0.77) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,0.65,1.3), yscale(log extend) legend(off) scheme(cleanplots) xlabel(0.65(0.05)1.3) xmtick(0.65(0.05)1.3) ylabel(1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("PA:AA") title("HR for mortality; Age 65+; adjusted for age(decade), sex, and comorbidities") yline(1, lp("shortdash") lc(gs10)) xline(0.9, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA:AA = 0.85) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (0.9) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAAA Splines - Age 65plus Adj AgeSexComorb.png", as(png) name("Graph") replace
restore
*/


/* TESTRUN 
------- Age AGE <69, ADJUSTED FOR SEX, AGE, AND COMORBIDITIES
//is this because 69 is the age where half of deaths occured above and below?
*/ 
/*
preserve
keep if age <= 69
gen mpaaa_rounded = round(mpaaa, 0.01)
mkspline2 rc = mpaaa_rounded, cubic nknots(4) displayknots            
// VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
assert float(mpaaa_rounded) == float(rc1)
// MODEL WITH FULL SPLINE
stcox rc* male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf, strata(age_decade)
levelsof mpaaa_rounded, local(levels)
xblc rc*, covname(mpaaa_rounded) at(`r(levels)') reference(0.81) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,0.65,1.3), yscale(log extend) legend(off) scheme(cleanplots) xlabel(0.65(0.05)1.3) xmtick(0.65(0.05)1.3) ylabel(1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("PA:AA") title("HR for mortality; Age <69; adjusted for age(decade), sex, and comorbidities") yline(1, lp("shortdash") lc(gs10)) xline(0.9, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA:AA = 0.85) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (0.9) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAAA Splines - Age below 50 Adj AgeSexComorb.png", as(png) name("Graph") replace
restore

/* TESTRUN 
------- Age AGE 69+, ADJUSTED FOR SEX, AGE, AND COMORBIDITIES
*/ 
preserve
keep if age > 69
gen mpaaa_rounded = round(mpaaa, 0.01)
mkspline2 rc = mpaaa_rounded, cubic nknots(4) displayknots            
// VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
assert float(mpaaa_rounded) == float(rc1)
// MODEL WITH FULL SPLINE
stcox rc* male peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf, strata(age_decade)
levelsof mpaaa_rounded, local(levels)
xblc rc*, covname(mpaaa_rounded) at(`r(levels)') reference(0.77) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,0.65,1.3), yscale(log extend) legend(off) scheme(cleanplots) xlabel(0.65(0.05)1.3) xmtick(0.65(0.05)1.3) ylabel(1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("PA:AA") title("HR for mortality; Age 69+; adjusted for age(decade), sex, and comorbidities") yline(1, lp("shortdash") lc(gs10)) xline(0.9, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA:AA = 0.85) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (0.9) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAAA Splines - Age 50 plus Adj AgeSexComorb.png", as(png) name("Graph") replace
restore
*/



// TODO: COULD YOU RE-INDEX THIS SO THAT THE BELOW-THRESHOLD AVERAGE IS TAKEN AS HR 1.0 AND A HORIZONTAL LINE IS DRAWN FOR ABOVE THE INDEX TO DEMONSTRATE HOW MUCH DATA THE DICHOTOMIZATION IS LOSING?


/*-----------------------------
---------PA_d ANALYSIS---------
-----------------------------*/

/*
------- WHOLE COHORT, UNADJUSTED
*/ 
/*
summarize sex_norm_mpad, detail
preserve
mkspline2 rc = sex_norm_mpad, cubic nknots(4) displayknots            
assert float(sex_norm_mpad) == float(rc1)
stcox rc*
estat ic
levelsof sex_norm_mpad, local(levels)
xblc rc*, covname(sex_norm_mpad) at(`r(levels)') reference(25.1) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,17.5,41.5), yscale(log extend) scheme(cleanplots) legend(off) xlabel(17.5(3)41.5) xmtick(17.5(1)41.5) ylabel(0.5 1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("Sex-Normalized PA_d (mm)") title("HR for mortality; entire cohort; Unadjusted") yline(1, lp("shortdash") lc(gs10)) xline(28, lp("shortdash_dot") lc(gs10)) note("Healthy Population mean (PA_d = 25.1mm, Truong 2012) taken as reference." "Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (28mm) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAd Splines - Whole Cohort Unadj.png", as(png) name("Graph") replace
restore
*/
 
/*
------- WHOLE COHORT, JUST ADJUSTED FOR SEX AND AGE
MAIN FIGURE FOR CHEST ABSTRACT
*/ 


preserve 
summarize sex_norm_mpad, detail

mkspline2 rc = sex_norm_mpad, cubic nknots(4) displayknots            
assert float(sex_norm_mpad) == float(rc1)
stcox rc*, strata(age_decade)
estat ic
levelsof sex_norm_mpad if inrange(sex_norm_mpad, 19, 39), local(levels)
xblc rc*, covname(sex_norm_mpad) at(`r(levels)') reference(25.1) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) ///
 (line or pa, sort lc(black) lp(l)) if inrange(pa,19,39), ///
 yscale(log extend) ///
 scheme(cleanplots) ///
 legend(off) ///
 xlabel(19(3)37, labsize(large)) ///
 xmtick(19(3)37) ///
 ylabel(1 2 4, angle(horiz) format(%2.1fc) labsize(large)) ///
 ytitle("Hazard Ratio of Mortality", size(large)) ///
 xtitle("Sex-Normalized Pulmonary Artery Diameter (mm)", size(large)) ///
 title("Age-adjusted Mortality Risk by PA Diameter", size(vlarge)) ///
 yline(1, lp("shortdash") lc(gs10)) ///
 xline(28, lp("shortdash_dot") lc(gs10)) ///
 text(0.85 33.75 "95% confidence interval (black dash)", size(medlarge)) ///
 text(2.5 24 "Referenced to" "healthy population" "PAd mean (25.1 mm)", size(medlarge))
graph export "Results and Figures/$S_DATE/HR PAd Splines - Whole Cohort Adj AgeSex.png", as(png) name("Graph") replace
restore

///  note("Lines added at traditional threshold of abnormal (28mm) and HR = 1." "Hazard Ratio Referenced to Population Mean (25.1mm)" "Dashed line = 95% Confidence interval")


/*
------- WHOLE COHORT, ADJUSTED FOR SEX, AGE, AND COMORBIDITIES
*/ 
/*
summarize sex_norm_mpad, detail
preserve
mkspline2 rc = sex_norm_mpad, cubic nknots(4) displayknots            
assert float(sex_norm_mpad) == float(rc1)
stcox rc* peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf, strata(age_decade)
levelsof sex_norm_mpad, local(levels)
xblc rc*, covname(sex_norm_mpad) at(`r(levels)') reference(25.1) eform generate(pa or lb ub)
estat ic
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,17.5,41.5), yscale(log extend) scheme(cleanplots) legend(off) xlabel(17.5(3)41.5) xmtick(17.5(1)41.5) ylabel(0.5 1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("Sex-Adjusted PA_d (mm)") title("HR for mortality; entire cohort; adjusted for age(decade), sex, and comorbidities") yline(1, lp("shortdash") lc(gs10)) xline(28, lp("shortdash_dot") lc(gs10)) note("Healthy Population mean (PA_d = 25.1mm, Truong 2012) taken as reference." "Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (28mm) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAd Splines - Whole Cohort Adj AgeSexComorb.png", as(png) name("Graph") replace
restore
*/


/*
------- Age < 50 , JUST ADJUSTED FOR SEX AND AGE
*/ 
summarize sex_norm_mpad, detail
preserve
keep if age_bi_decade == 0
mkspline2 rc = sex_norm_mpad, cubic nknots(4) displayknots            
// VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
assert float(sex_norm_mpad) == float(rc1)
// MODEL WITH FULL SPLINE
stcox rc* , strata(age_decade)
estat ic
levelsof sex_norm_mpad, local(levels)
xblc rc*, covname(sex_norm_mpad) at(`r(levels)') reference(25.5) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,17.5,41.5), yscale(log extend) scheme(cleanplots) legend(off) xlabel(17.5(3)41.5) xmtick(17.5(1)41.5) ylabel(0.5 1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("Sex-Adjusted PA_d (mm)") title("HR for mortality; Age < 50; adjusted for age(decade) and sex") yline(1, lp("shortdash") lc(gs10)) xline(28, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA_d = 25.5mm) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (28mm) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAd Splines - Age less 50 Adj AgeSex.png", as(png) name("Graph") replace
restore

/*
------- Age < 50, ADJUSTED FOR SEX, AGE, AND COMORBIDITIES
*/ 
/*
summarize sex_norm_mpad, detail
preserve
keep if age_bi_decade == 0
mkspline2 rc = sex_norm_mpad, cubic nknots(4) displayknots            
// VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
assert float(sex_norm_mpad) == float(rc1)
// MODEL WITH FULL SPLINE
stcox rc* peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf, strata(age_decade)
levelsof sex_norm_mpad, local(levels)
xblc rc*, covname(sex_norm_mpad) at(`r(levels)') reference(25.5) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,17.5,41.5), yscale(log extend) scheme(cleanplots) legend(off) xlabel(17.5(3)41.5) xmtick(17.5(1)41.5) ylabel(0.5 1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("Sex-Adjusted PA_d (mm)") title("HR for mortality; Age < 50; adjusted for age(decade), sex, and comorbidities") yline(1, lp("shortdash") lc(gs10)) xline(28, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA_d = 25.5mm) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (28mm) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAd Splines - Age less 50 Adj AgeSexComorb.png", as(png) name("Graph") replace
restore
*/

/*
------- Age 50-65 , JUST ADJUSTED FOR SEX AND AGE
*/ 
summarize sex_norm_mpad, detail
preserve
keep if age_bi_decade == 1
mkspline2 rc = sex_norm_mpad, cubic nknots(4) displayknots            
// VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
assert float(sex_norm_mpad) == float(rc1)
// MODEL WITH FULL SPLINE
stcox rc*, strata(age_decade)
levelsof sex_norm_mpad, local(levels)
xblc rc*, covname(sex_norm_mpad) at(`r(levels)') reference(25.5) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,17.5,41.5), yscale(log extend) scheme(cleanplots) legend(off) xlabel(17.5(3)41.5) xmtick(17.5(1)41.5) ylabel(0.5 1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("Sex-Adjusted PA_d (mm)") title("HR for mortality; Age 50-65; adjusted for age(decade) and sex") yline(1, lp("shortdash") lc(gs10)) xline(28, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA_d = 25.5mm) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (28mm) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAd Splines - Age 50-65 Adj AgeSex.png", as(png) name("Graph") replace
restore

/*
------- Age 50-65, ADJUSTED FOR SEX, AGE, AND COMORBIDITIES
*/ 
/*
summarize sex_norm_mpad, detail
preserve
keep if age_bi_decade == 1
mkspline2 rc = sex_norm_mpad, cubic nknots(4) displayknots            
// VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
assert float(sex_norm_mpad) == float(rc1)
// MODEL WITH FULL SPLINE
stcox rc* peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf, strata(age_decade)
levelsof sex_norm_mpad, local(levels)
xblc rc*, covname(sex_norm_mpad) at(`r(levels)') reference(25.5) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,17.5,41.5), yscale(log extend) scheme(cleanplots) legend(off) xlabel(17.5(3)41.5) xmtick(17.5(1)41.5) ylabel(0.5 1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("Sex-Adjusted PA_d (mm)") title("HR for mortality; Age 50-65; adjusted for age(decade), sex, and comorbidities") yline(1, lp("shortdash") lc(gs10)) xline(28, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA_d = 25.5mm) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (28mm) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAd Splines - Age 50-65 Adj AgeSexComorb.png", as(png) name("Graph") replace
restore
*/

/*
------- Age 65+ , JUST ADJUSTED FOR SEX AND AGE
*/ 
summarize sex_norm_mpad, detail
preserve
keep if age_bi_decade == 2
mkspline2 rc = sex_norm_mpad, cubic nknots(4) displayknots            
// VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
assert float(sex_norm_mpad) == float(rc1)
// MODEL WITH FULL SPLINE
stcox rc*, strata(age_decade)
levelsof sex_norm_mpad, local(levels)
xblc rc*, covname(sex_norm_mpad) at(`r(levels)') reference(25.5) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,17.5,41.5), yscale(log extend) scheme(cleanplots) legend(off) xlabel(17.5(3)41.5) xmtick(17.5(1)41.5) ylabel(0.5 1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("Sex-Adjusted PA_d (mm)") title("HR for mortality; Age 65+; adjusted for age(decade) and sex") yline(1, lp("shortdash") lc(gs10)) xline(28, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA_d = 25.5mm) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (28mm) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAd Splines - Age 65 plus Adj AgeSex.png", as(png) name("Graph") replace
restore

/*
------- Age 65+, ADJUSTED FOR SEX, AGE, AND COMORBIDITIES
*/ 
/*
summarize sex_norm_mpad, detail
preserve
keep if age_bi_decade == 2
mkspline2 rc = sex_norm_mpad, cubic nknots(4) displayknots            
// VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
assert float(sex_norm_mpad) == float(rc1)
// MODEL WITH FULL SPLINE
stcox rc* peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf, strata(age_decade)
levelsof sex_norm_mpad, local(levels)
xblc rc*, covname(sex_norm_mpad) at(`r(levels)') reference(25.5) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,17.5,41.5), yscale(log extend) scheme(cleanplots) legend(off) xlabel(17.5(3)41.5) xmtick(17.5(1)41.5) ylabel(0.5 1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("Sex-Adjusted PA_d (mm)") title("HR for mortality; Age 65+; adjusted for age(decade), sex, and comorbidities") yline(1, lp("shortdash") lc(gs10)) xline(28, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA_d = 25.5mm) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (28mm) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAd Splines - Age 65 plus Adj AgeSexComorb.png", as(png) name("Graph") replace
restore
*/

/*TESTRUN
------- Age <69, ADJUSTED FOR SEX, AGE, AND COMORBIDITIES
*/ 
/*
summarize sex_norm_mpad, detail
preserve
keep if age <= 69
mkspline2 rc = sex_norm_mpad, cubic nknots(4) displayknots            
// VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
assert float(sex_norm_mpad) == float(rc1)
// MODEL WITH FULL SPLINE
stcox rc*, strata(age_decade)
levelsof sex_norm_mpad, local(levels)
xblc rc*, covname(sex_norm_mpad) at(`r(levels)') reference(25.5) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,17.5,41.5), yscale(log extend) scheme(cleanplots) legend(off) xlabel(17.5(3)41.5) xmtick(17.5(1)41.5) ylabel(0.5 1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("Sex-Adjusted PA_d (mm)") title("HR for mortality; Age <69; adjusted for age(decade) and sex") yline(1, lp("shortdash") lc(gs10)) xline(28, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA_d = 25.5mm) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (28mm) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAd Splines - Age below 69 Adj AgeSexComorb.png", as(png) name("Graph") replace
restore 

/*TESTRUN
------- Age <69, ADJUSTED FOR SEX, AGE, AND COMORBIDITIES
*/ 
/*
summarize sex_norm_mpad, detail
preserve
keep if age <= 69
mkspline2 rc = sex_norm_mpad, cubic nknots(4) displayknots            
// VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
assert float(sex_norm_mpad) == float(rc1)
// MODEL WITH FULL SPLINE
stcox rc* peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf, strata(age_decade)
levelsof sex_norm_mpad, local(levels)
xblc rc*, covname(sex_norm_mpad) at(`r(levels)') reference(25.5) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,17.5,41.5), yscale(log extend) scheme(cleanplots) legend(off) xlabel(17.5(3)41.5) xmtick(17.5(1)41.5) ylabel(0.5 1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("Sex-Adjusted PA_d (mm)") title("HR for mortality; Age <69; adjusted for age(decade), sex, and comorbidities") yline(1, lp("shortdash") lc(gs10)) xline(28, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA_d = 25.5mm) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (28mm) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAd Splines - Age below 69 Adj AgeSexComorb.png", as(png) name("Graph") replace
restore 


/*TESTRUN
------- Age 69 plus, ADJUSTED FOR SEX, AGE, AND COMORBIDITIES
*/ 

summarize sex_norm_mpad, detail
preserve
keep if age > 69
mkspline2 rc = sex_norm_mpad, cubic nknots(4) displayknots            
// VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
assert float(sex_norm_mpad) == float(rc1)
// MODEL WITH FULL SPLINE
stcox rc* peripheralvascdisorders obesity_calc diabetes1 diabetes2 renalfailure pulmdisease pulmonarycircdisorders chf, strata(age_decade)
levelsof sex_norm_mpad, local(levels)
xblc rc*, covname(sex_norm_mpad) at(`r(levels)') reference(25.5) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,17.5,41.5), yscale(log extend) scheme(cleanplots) legend(off) xlabel(17.5(3)41.5) xmtick(17.5(1)41.5) ylabel(0.5 1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("Sex-Adjusted PA_d (mm)") title("HR for mortality; Age 69+; adjusted for age(decade), sex, and comorbidities") yline(1, lp("shortdash") lc(gs10)) xline(28, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA_d = 25.5mm) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (28mm) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAd Splines - Age 69 plus Adj AgeSexComorb.png", as(png) name("Graph") replace
restore
*/ 
*/ 

/* -------------------
Visualizing the effect of age on the HR by PA size
---------------------*/ 

//Using age (3 knot spline) and PA size (3 knot spline) 
//Not sure this is really capturing the relationship we want.. because it's not allowing a different shape of the relationship a different times

preserve
mkspline2 rcage = age, cubic nknots(3) displayknots
assert float(age) == float(rcage1) 
mkspline2 rcpad = sex_norm_mpad, cubic nknots(3) displayknots            
assert float(sex_norm_mpad) == float(rcpad1)
stcox rcage* rcpad* 
predict age_pad_spline_model, xb 
heatplot age_pad_spline_model age sex_norm_mpad, aspectratio(0.7) color(YlOrRd) cut(1(0.5)5.5) xlabel(15(05)40, angle(vertical) labsize(3))   xbwidth(2) xtitle("Sex-Adjusted PA_d (mm)", size(4)) ybwidth(5) ytitle("Age (years)", size(4)) ylabel(20(10)100 ,labsize(3)) ramp(right format(%3.2f) space(18) subtitle("Predicted" "HR of" "Mortality", size(small) justification(center)) label(1(0.5)5.5, labsize(3))) p(lcolor(black%10) lwidth(*0.15)) clip     
graph export "Results and Figures/$S_DATE/Predicted HR Age and PAd w 3 knot splines.png", as(png) name("Graph") replace
restore

preserve
mkspline2 rcage = age, cubic nknots(3) displayknots
assert float(age) == float(rcage1) 
gen mpaaa_rounded = round(mpaaa, 0.01)
mkspline2 rcpaaa = mpaaa_rounded, cubic nknots(3) displayknots            
assert float(mpaaa_rounded) == float(rcpaaa1)
stcox rcage* rcpaaa* male
predict age_paaa_spline_model, xb 
heatplot age_paaa_spline_model age mpaaa_rounded, cut(1(0.5)8) aspectratio(0.7) color(YlOrRd) xlabel(0.5(0.1)1.4, angle(vertical) labsize(3)) xbwidth(0.05) xtitle("PA:AA ratio", size(4)) ybwidth(5) ytitle("Age (years)", size(4)) ylabel(20(10)100 ,labsize(3)) ramp(right format(%3.2f) space(18) subtitle("Predicted" "HR of" "Mortality", size(small) justification(center)) label(1(0.5)8, labsize(3))) p(lcolor(black%10) lwidth(*0.15)) clip     
graph export "Results and Figures/$S_DATE/Predicted HR Age and PAAA w 3 knot splines.png", as(png) name("Graph") replace
restore



/* Modeling 1 large stratified cox-regression by bi-decade, with an interaction. 

I believe this is still insufficient in allowing separate shapes per strata
*/ 

/*

COMPLETING FROM HERE DOWN - NEW PLAN IS TO MODEL WITH XBLC

todo: figure out axes
todo: this seems quite sensitive to the # of knots. why?

*/

preserve
clear
use cleaned_noempi
stset lastfollowupyear, failure(death==1)


/* COMBINED ANALYSIS OF PAD AND PAAA - all adjusted for age and sex*/ 

//Step 1: generate levels (apply to all age strata)
gen mpaaa_rounded = round(mpaaa, 0.01)
//sex_norm_mpad is already rounded.

//Step 2: mksplines for each age strata by metric
//Age < 50
mkspline2 rc_mpaaa_young = mpaaa_rounded if age_bi_decade == 0, cubic nknots(3) displayknots            
assert float(mpaaa_rounded) == float(rc_mpaaa_young1) if age_bi_decade == 0
mkspline2 rc_mpad_young = sex_norm_mpad if age_bi_decade == 0, cubic nknots(3) displayknots            
assert float(sex_norm_mpad) == float(rc_mpad_young1) if age_bi_decade == 0

//Age 50 - 65
mkspline2 rc_mpaaa_mid = mpaaa_rounded if age_bi_decade == 1, cubic nknots(3) displayknots            
assert float(mpaaa_rounded) == float(rc_mpaaa_mid1) if age_bi_decade == 1
mkspline2 rc_mpad_mid = sex_norm_mpad if age_bi_decade == 1, cubic nknots(3) displayknots            
assert float(sex_norm_mpad) == float(rc_mpad_mid1) if age_bi_decade == 1

// Age 65+
mkspline2 rc_mpaaa_old = mpaaa_rounded if age_bi_decade == 2, cubic nknots(3) displayknots            
assert float(mpaaa_rounded) == float(rc_mpaaa_old1) if age_bi_decade == 2
mkspline2 rc_mpad_old = sex_norm_mpad if age_bi_decade == 2, cubic nknots(3) displayknots            
assert float(sex_norm_mpad) == float(rc_mpad_old1) if age_bi_decade == 2


//Step 3: perform regressions and xblc predictions

//Age < 50
stcox rc_mpaaa_young* if age_bi_decade == 0, strata(age_decade)
levelsof mpaaa_rounded if age_bi_decade == 0 & inrange(mpaaa_rounded,0.65,1.2) , local(levels)
xblc rc_mpaaa_young*, covname(mpaaa_rounded) at(`r(levels)') reference(0.77) eform generate(mpaaa_young hr_mpaaa_young lb_mpaaa_young ub_mpaaa_young)

stcox rc_mpad_young* if age_bi_decade == 0, strata(age_decade)
levelsof sex_norm_mpad if age_bi_decade == 0 & inrange(sex_norm_mpad,19,39) , local(levels)
xblc rc_mpad_young*, covname(sex_norm_mpad) at(`r(levels)') reference(25.5) eform generate(mpad_young hr_mpad_young lb_mpad_young ub_mpad_young)

//Age 50-65
stcox rc_mpaaa_mid* if age_bi_decade == 1, strata(age_decade)
levelsof mpaaa_rounded if age_bi_decade == 1 & inrange(mpaaa_rounded,0.65,1.2), local(levels)
xblc rc_mpaaa_mid*, covname(mpaaa_rounded) at(`r(levels)') reference(0.77) eform generate(mpaaa_mid hr_mpaaa_mid lb_mpaaa_mid ub_mpaaa_mid)

stcox rc_mpad_mid* if age_bi_decade == 1, strata(age_decade)
levelsof sex_norm_mpad if age_bi_decade == 1 & inrange(sex_norm_mpad,19 ,39), local(levels)
xblc rc_mpad_mid*, covname(sex_norm_mpad) at(`r(levels)') reference(25.5) eform generate(mpad_mid hr_mpad_mid lb_mpad_mid ub_mpad_mid)

//Age 65+
stcox rc_mpaaa_old* if age_bi_decade == 2, strata(age_decade)
levelsof mpaaa_rounded if age_bi_decade == 2 & inrange(mpaaa_rounded,0.65,1.2), local(levels)
xblc rc_mpaaa_old*, covname(mpaaa_rounded) at(`r(levels)') reference(0.77) eform generate(mpaaa_old hr_mpaaa_old lb_mpaaa_old ub_mpaaa_old)

stcox rc_mpad_old* if age_bi_decade == 2, strata(age_decade)
levelsof sex_norm_mpad if age_bi_decade == 2 & inrange(sex_norm_mpad,19,39), local(levels)
xblc rc_mpad_old*, covname(sex_norm_mpad) at(`r(levels)') reference(25.5) eform generate(mpad_old hr_mpad_old lb_mpad_old ub_mpad_old)

//Step 5:  Create masks: 
//NOTE: these are crafted to ONLY get 1 block - need to be adjusted with data changes
recode hr_mpaaa_young (min/1.08=.)(1.08/1.12 = .)(1.1201/1.24=.)(1.24/1.26 = 1.25)(1.26/1.45=.)(1.45/1.55=1.5)(1.55/1.9=.)(1.9/2.0=2)(2.0/max=.), generate(masked_hr_mpaaa_young)
gen masked_hr_mpaaa_young_label = "HR = " + string(masked_hr_mpaaa_young)
list mpaaa_young hr_mpaaa_young masked_hr_mpaaa_young masked_hr_mpaaa_young_label

recode hr_mpaaa_mid (min/1.08=.)(1.08/1.12 = .)(1.1201/1.24=.)(1.24/1.3 = 1.25)(1.3/1.45=.)(1.45/1.55=1.5)(1.55/1.95=.)(1.95/2.05=2)(2.05/max=.), generate(masked_hr_mpaaa_mid)
gen masked_hr_mpaaa_mid_label = "HR = " + string(masked_hr_mpaaa_mid)
list mpaaa_mid hr_mpaaa_mid masked_hr_mpaaa_mid

recode hr_mpaaa_old (min/1.08=.)(1.08/1.12 = .)(1.1201/1.24=.)(1.24/1.26 = 1.25)(1.26/1.49=.)(1.49/1.51=1.5)(1.51/1.9=.)(1.9/2.0=2)(2.0/max=.), generate(masked_hr_mpaaa_old)
gen masked_hr_mpaaa_old_label = "HR = " + string(masked_hr_mpaaa_old)
list mpaaa_old hr_mpaaa_old masked_hr_mpaaa_old

recode hr_mpad_young (min/1.08=.)(1.08/1.12 = .)(1.1201/1.246=.)(1.246/1.248 = 1.25)(1.248/1.5=.)(1.5/1.51=1.5)(1.51/1.95=.)(1.95/1.96=2)(1.96/max=.), generate(masked_hr_mpad_young)
gen masked_hr_mpad_young_label = "HR = " + string(masked_hr_mpad_young)
list mpad_young hr_mpad_young masked_hr_mpad_young

recode hr_mpad_mid (min/1.08=.)(1.08/1.12 = .)(1.1201/1.25=.)(1.25/1.26 = 1.25)(1.26/1.5=.)(1.5/1.51=1.5)(1.51/1.98=.)(1.98/1.99=2)(1.99/max=.), generate(masked_hr_mpad_mid)
gen masked_hr_mpad_mid_label = "HR = " + string(masked_hr_mpad_mid)
list mpad_mid hr_mpad_mid masked_hr_mpad_mid

recode hr_mpad_old (min/1.08=.)(1.08/1.12 = .)(1.1201/1.249=.)(1.249/1.25 = 1.25)(1.25/1.42=.)(1.42/1.55=1.5)(1.55/1.9=.)(1.9/2.1=2)(2.1/max=.), generate(masked_hr_mpad_old)
gen masked_hr_mpad_old_label = "HR = " + string(masked_hr_mpad_old)
list mpad_old hr_mpad_old masked_hr_mpad_old masked_hr_mpad_old_label

//Step 6: Graph with CIs [two-way]
//Young: forest_green, Mid: gold, Old: Cranberry

//MPAAA
twoway (line lb_mpaaa_young ub_mpaaa_young mpaaa_young, sort lc(forest_green%33 forest_green%33) lp(longdash longdash)) (line hr_mpaaa_young mpaaa_young, sort lc(forest_green) lp(l)) || (line lb_mpaaa_mid ub_mpaaa_mid mpaaa_mid, sort lc(gold%33 gold%33) lp(longdash longdash)) (line hr_mpaaa_mid mpaaa_mid, sort lc(gold) lp(l))  || (line lb_mpaaa_old ub_mpaaa_old mpaaa_old, sort lc(cranberry%33 cranberry%33) lp(longdash longdash)) (line hr_mpaaa_old mpaaa_old, sort lc(cranberry) lp(l)), yscale(log extend) scheme(cleanplots) legend(off) xlabel(0.7(0.05)1.2) xmtick(0.7(0.05)1.2) ylabel(1 2 4, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("Pulmonary Artery to Ascending Aorta Ratio") title("Entire cohort Stratified by Age; adjusted for gender") yline(1, lp("shortdash") lc(gs10)) xline(0.9, lp("shortdash_dot") lc(gs10)) note("Lines added at traditional threshold of abnormal (0.9) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAAA by age strat with CI.png", as(png) name("Graph") replace

//MPAD
twoway (line lb_mpad_young ub_mpad_young mpad_young, sort lc(forest_green%33 forest_green%33) lp(longdash longdash)) (line hr_mpad_young mpad_young, sort lc(forest_green) lp(l)) || (line lb_mpad_mid ub_mpad_mid mpad_mid, sort lc(gold%33 gold%33) lp(longdash longdash)) (line hr_mpad_mid mpad_mid, sort lc(gold) lp(l)) || (line lb_mpad_old ub_mpad_old mpad_old, sort lc(cranberry%33 cranberry%33) lp(longdash longdash)) (line hr_mpad_old mpad_old, sort lc(cranberry) lp(l)), yscale(log extend) scheme(cleanplots) legend(off) xlabel(19(2)39) xmtick(19(2)39) ylabel(0.5 1 2 4, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("Sex-Adjusted PA_d (mm)") title("HR for mortality; Age < 50; adjusted for age(decade) and sex") yline(1, lp("shortdash") lc(gs10)) xline(28, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA_d = 25.5mm) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (28mm) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR MPAD by age strat with CI.png", as(png) name("Graph") replace

/* PA:AA Heatplots */ 
gen young_axis = "Age < 50"
heatplot masked_hr_mpaaa_young young_axis mpaaa_young, ///
 values(label(masked_hr_mpaaa_young_label) angle(90) size(5)) ///
 color(carto Geyser) /// 
 legend(off) /// 
 xbwidth(0.02) ///
 xtitle("") ///
 xlabel(0.7(0.1)1.2, labsize(large)) /// 
 xline(0.77 0.9) ///
 cuts(0.5(0.025)2) ///
 ytitle("") /// 
 ylabel(,labsize(vlarge) angle(90)) 
graph save "temp_mpaaa_young.gph", replace

gen mid_axis = "Age 50-65"
heatplot masked_hr_mpaaa_mid mid_axis mpaaa_mid, ///
 values(label(masked_hr_mpaaa_mid_label) angle(90) size(5)) ///
 color(carto Geyser) ///
 xbwidth(0.02) xtitle("") ///
 xlabel(0.7(0.1)1.2, labsize(large)) ///
 xline(0.77 0.9) ///
 cuts(0.5(0.025)2) ///
 legend(off) ///
 ytitle("") ///
 ylabel(,labsize(vlarge) angle(90)) 
graph save "temp_mpaaa_mid.gph", replace

gen old_axis = "Age > 65"
heatplot masked_hr_mpaaa_old old_axis mpaaa_old, ///
 values(label(masked_hr_mpaaa_old_label) angle(90) size(5)) ///
 color(carto Geyser) /// 
 legend(off) /// 
 xbwidth(0.02) xtitle("") ///
 xlabel(0.7(0.1)1.2, labsize(large)) /// 
 xline(0.77 0.9) ///
 cuts(0.5(0.025)2) ///
 ytitle("") /// 
 ylabel(, labsize(vlarge) angle(90))
graph save "temp_mpaaa_old.gph", replace

graph combine temp_mpaaa_young.gph temp_mpaaa_mid.gph temp_mpaaa_old.gph, ///
 xcommon col(1) ///
 imargin(0 0 0 0) graphregion(margin(l=22 r=22)) ///
 title("PA:AA Thresholds for Age-Stratified Mortality Risk", size(large)) ///
 note("Pulmonary Artery to Ascending Aorta Ratio", pos(6) size(medlarge))
graph export "Results and Figures/$S_DATE/HR by Age strat and PAAA w 3 knot spline.png", as(png) name("Graph") replace

/* PAd Heatplots */ 
heatplot masked_hr_mpad_young young_axis mpad_young, ///
 values(label(masked_hr_mpad_young_label) angle(90) size(5)) ///
 color(carto Geyser) /// 
 legend(off) /// 
 xbwidth(0.8) ///
 xtitle("") ///
 xlabel(19(3)40, labsize(large)) /// 
 xline(25.1 28) ///
 cuts(0.5(0.025)2) ///
 ytitle("") ///
 ylabel(,labsize(vlarge) angle(90)) 
graph save "temp_mpad_young.gph", replace

heatplot masked_hr_mpad_mid mid_axis mpad_mid, ///
 values(label(masked_hr_mpad_mid_label) angle(90) size(5)) ///
 color(carto Geyser) ///
 xbwidth(0.8) xtitle("") ///
 xlabel(19(3)40, labsize(large)) ///
 xline(25.1 28) ///
 cuts(0.5(0.025)2) ///
 legend(off) ///
 ytitle("") ///
 ylabel(,labsize(vlarge) angle(90)) 
graph save "temp_mpad_mid.gph", replace

heatplot masked_hr_mpad_old old_axis mpad_old, ///
 values(label(masked_hr_mpad_old_label) angle(90) size(5)) ///
 color(carto Geyser) /// 
 legend(off) /// 
 xbwidth(0.8) xtitle("") ///
 xlabel(19(3)40, labsize(large)) /// 
 xline(25.1 28) ///
 cuts(0.5(0.025)2) ///
 ytitle("") ///
 ylabel(, angle(90) labsize(vlarge))
graph save "temp_mpad_old.gph", replace

graph combine temp_mpad_young.gph temp_mpad_mid.gph temp_mpad_old.gph, ///
 xcommon col(1) ///
 imargin(0 0 0 0) graphregion(margin(l=22 r=22)) ///
 title("PAd Thresholds for Age-Stratified Mortality Risk", size(large)) ///
 note("Sex-normalized Pulmonary Artery Diameter (mm)", pos(6) size(medlarge))
graph export "Results and Figures/$S_DATE/HR by Age strat and PAd w 3 knot spline.png", as(png) name("Graph") replace


restore





//note: 
//xtitle (PA:AA Ratio - maybe add this to wrapper?)
///   ramp(bottom format(%3.1f) space(18) subtitle("Color Key:" "Hazard Ratio", size(medsmall) justification(center)) label(0.75 1 1.5 2, labsize(3)))

//TODO: 
//[ ] consider a way to add a legend: perhaps using fixed aspect ratios of combine or by adding labels to specific colors [ie. HR 1, HR 2 etc.]
//[ ] tweak colors to be more sensible - https://repec.sowi.unibe.ch/stata/palettes/colors.html perhaps use a custom interpolation, or figure a way to set the diverging value
//[ ] figure how to move vertical line to the front
//[ ] possible to get rid of axes?


//scratch
//   ramp(bottom format(%3.1f) space(18) subtitle("Color Key:" "Hazard Ratio", size(medsmall) justification(center)) label(0.75 1 1.5 2, labsize(3)))










/* WORK IN PROGRESS */ 


/*

heatplot mpad_model_synth i.age_bi_decade sex_norm_mpad, aspectratio(0.7) color(YlOrRd) xlabel(15(05)40, angle(vertical) labsize(3)) xbwidth(2) xtitle("Sex-Adjusted PA_d (mm)", size(4)) ytitle("Age Category", size(4))  ramp(right format(%3.2f) space(18) subtitle("Predicted" "HR of" "Mortality", size(small) justification(center)) label(1(0.5)5.5, labsize(3))) p(lcolor(black%10) lwidth(*0.15)) clip     
*/



preserve 
gen mpaaa_rounded = round(mpaaa, 0.01)
mkspline2 rcpaaa = mpaaa_rounded, cubic nknots(3) displayknots            
assert float(mpaaa_rounded) == float(rcpaaa1)
stcox ib1.age_bi_decade##c.(rcpaaa*) male
levelsof mpaaa_rounded, local(levels)
xblc rcpaaa*, covname(mpaaa_rounded) at(`r(levels)') reference(0.77) eform generate(pa hr lb ub)

//it just writes this to new columns at each of the levels. 

twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line hr pa, sort lc(black) lp(l)) if inrange(pa,0.65,1.3), yscale(log extend) scheme(cleanplots) legend(off) xlabel(0.65(0.05)1.3) xmtick(0.65(0.05)1.3) ylabel(1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("PA:AA") title("HR for mortality; Age 50-65; adjusted for age(decade) and sex") yline(1, lp("shortdash") lc(gs10)) xline(0.9, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA:AA = 0.85) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (0.9) and HR = 1. Dashed line = 95% Confidence interval")

//predict strat_age_paaa_spline_model, hr 
list age male mpaaa_rounded hr

restore


clear
use cleaned_noempi
stset lastfollowupyear, failure(death==1)

/*
heatplot strat_age_paaa_spline_model i.age_bi_decade mpaaa_rounded, cut(0.25(0.05)2) aspectratio(0.7) color(YlOrRd) xlabel(0.5(0.1)1.4, angle(vertical) labsize(3)) xbwidth(0.05) xtitle("PA:AA ratio", size(4)) ytitle("Age (years)", size(4)) ramp(right format(%3.2f) space(18) subtitle("Predicted" "HR of" "Mortality", size(small) justification(center)) label(0.25(0.25)2, labsize(3))) p(lcolor(black%10) lwidth(*0.15)) clip     
graph export "Results and Figures/$S_DATE/Predicted HR Strat-Age and PAAA w 3 knot spline.png", as(png) name("Graph") replace
*/

//TODO: 

//Make separate 4-knot cox models [as above] for each age group and MPAA/PAd
//ensure the predictions are normalized vs the lowest HR in the age group  [perhaps predict each into a separate variable, then subsequently combine them]
//Heat map by bi-decade (y-axis) and PA size (x-axis, many cuts. )


clear
use cleaned_noempi
stset lastfollowupyear, failure(death==1)

/*
------- MPAD Age < 50 , JUST ADJUSTED FOR SEX AND AGE
*/ 
mkspline2 rc_mpad_young = sex_norm_mpad if age_bi_decade == 0, cubic nknots(4) displayknots            
assert float(sex_norm_mpad) == float(rc_mpad_young1) if age_bi_decade == 0
stcox rc_mpad_young* if age_bi_decade == 0, strata(age_decade)
predict mpad_young_model, hr
//Comment out to test
/*
levelsof sex_norm_mpad if age_bi_decade == 0, local(levels)
xblc rc_mpad_young*, covname(sex_norm_mpad) at(`r(levels)') reference(25.5) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,17.5,41.5), yscale(log extend) scheme(cleanplots) legend(off) xlabel(17.5(3)41.5) xmtick(17.5(1)41.5) ylabel(0.5 1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("Sex-Adjusted PA_d (mm)") title("HR for mortality; Age < 50; adjusted for age(decade) and sex") yline(1, lp("shortdash") lc(gs10)) xline(28, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA_d = 25.5mm) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (28mm) and HR = 1. Dashed line = 95% Confidence interval")
*/ 

/*
------- MPAD Age 50-65 , JUST ADJUSTED FOR SEX AND AGE
*/ 
mkspline2 rc_mpad_mid = sex_norm_mpad if age_bi_decade == 1, cubic nknots(4) displayknots
assert float(sex_norm_mpad) == float(rc_mpad_mid1) if age_bi_decade == 1
stcox rc_mpad_mid* if age_bi_decade == 1, strata(age_decade)
predict mpad_mid_model, hr

//Comment out to test
/*
levelsof sex_norm_mpad if age_bi_decade == 1, local(levels)
xblc rc*, covname(sex_norm_mpad) at(`r(levels)') reference(25.5) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,17.5,41.5), yscale(log extend) scheme(cleanplots) legend(off) xlabel(17.5(3)41.5) xmtick(17.5(1)41.5) ylabel(0.5 1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("Sex-Adjusted PA_d (mm)") title("HR for mortality; Age 50-65; adjusted for age(decade) and sex") yline(1, lp("shortdash") lc(gs10)) xline(28, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA_d = 25.5mm) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (28mm) and HR = 1. Dashed line = 95% Confidence interval")
*/ 


/*
------- MPAD Age 65+ , JUST ADJUSTED FOR SEX AND AGE
*/ 
mkspline2 rc_mpad_old = sex_norm_mpad if age_bi_decade == 2, cubic nknots(4) displayknots            
// VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
assert float(sex_norm_mpad) == float(rc_mpad_old1) if age_bi_decade == 2
// MODEL WITH FULL SPLINE
stcox rc_mpad_old* if age_bi_decade == 2, strata(age_decade)
predict mpad_old_model, hr

//Comment out to test 
/*
levelsof sex_norm_mpad if age_bi_decade == 2, local(levels)
xblc rc*, covname(sex_norm_mpad) at(`r(levels)') reference(25.5) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,17.5,41.5), yscale(log extend) scheme(cleanplots) legend(off) xlabel(17.5(3)41.5) xmtick(17.5(1)41.5) ylabel(0.5 1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("Sex-Adjusted PA_d (mm)") title("HR for mortality; Age 65+; adjusted for age(decade) and sex") yline(1, lp("shortdash") lc(gs10)) xline(28, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA_d = 25.5mm) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (28mm) and HR = 1. Dashed line = 95% Confidence interval")
*/ 

egen mpad_model_synth = rmean(mpad_young_model mpad_mid_model mpad_old_model) //only 1 for each so just combines
list mpad_young_model mpad_mid_model mpad_old_model mpad_model_synth 

heatplot mpad_model_synth i.age_bi_decade sex_norm_mpad, aspectratio(0.7) color(YlOrRd) xlabel(15(05)40, angle(vertical) labsize(3)) xbwidth(2) xtitle("Sex-Adjusted PA_d (mm)", size(4)) ytitle("Age Category", size(4))  ramp(right format(%3.2f) space(18) subtitle("Predicted" "HR of" "Mortality", size(small) justification(center)) label(1(0.5)5.5, labsize(3))) p(lcolor(black%10) lwidth(*0.15)) clip     
graph export "Results and Figures/$S_DATE/Predicted HR Age Cat and PAd w 4 knot splines.png", as(png) name("Graph") replace 

// TODO: these seem improperly normalized, then re-add cut(1(0.5)5.5) yscale(0.5(1)3.5)


/*
-------MPAAa Age < 50 , JUST ADJUSTED FOR SEX AND AGE
*/ 
preserve
keep if age_bi_decade == 0
gen mpaaa_rounded = round(mpaaa, 0.01)
mkspline2 rc = mpaaa_rounded, cubic nknots(3) displayknots            
// VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
assert float(mpaaa_rounded) == float(rc1)
// MODEL WITH FULL SPLINE
tab death male
stcox rc* male // c.age //strata(age_decade)
levelsof mpaaa_rounded, local(levels)
xblc rc*, covname(mpaaa_rounded) at(`r(levels)') reference(0.77) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,0.65,1.3), yscale(log extend) scheme(cleanplots) legend(off) xlabel(0.65(0.05)1.3) xmtick(0.65(0.05)1.3) ylabel(1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("PA:AA") title("HR for mortality; Age < 50; adjusted for age(decade) and sex") yline(1, lp("shortdash") lc(gs10)) xline(0.9, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA:AA = 0.85) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (0.9) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAAA Splines - Age less 50 Adj AgeSex.png", as(png) name("Graph") replace
restore


/*
------- Age 50-65 , JUST ADJUSTED FOR SEX AND AGE
*/ 
preserve
keep if age_bi_decade == 1
gen mpaaa_rounded = round(mpaaa, 0.01)
mkspline2 rc = mpaaa_rounded, cubic nknots(4) displayknots            
// VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE
assert float(mpaaa_rounded) == float(rc1)
// MODEL WITH FULL SPLINE
stcox rc* male, strata(age_decade)
levelsof mpaaa_rounded, local(levels)
xblc rc*, covname(mpaaa_rounded) at(`r(levels)') reference(0.77) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,0.65,1.3), yscale(log extend) scheme(cleanplots) legend(off) xlabel(0.65(0.05)1.3) xmtick(0.65(0.05)1.3) ylabel(1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("PA:AA") title("HR for mortality; Age 50-65; adjusted for age(decade) and sex") yline(1, lp("shortdash") lc(gs10)) xline(0.9, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA:AA = 0.85) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (0.9) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAAA Splines - Age 50-65 Adj AgeSex.png", as(png) name("Graph") replace
restore



/*
-------MPAAA Age 65+, JUST ADJUSTED FOR SEX AND AGE
*/ 
preserve
keep if age_bi_decade == 2
gen mpaaa_rounded = round(mpaaa, 0.01)
mkspline2 rc = mpaaa_rounded, cubic nknots(4) displayknots            
assert float(mpaaa_rounded) == float(rc1)
stcox rc* male, strata(age_decade)
levelsof mpaaa_rounded, local(levels)
xblc rc*, covname(mpaaa_rounded) at(`r(levels)') reference(0.77) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) (line or pa, sort lc(black) lp(l)) if inrange(pa,0.65,1.3), yscale(log extend) scheme(cleanplots) legend(off) xlabel(0.65(0.05)1.3) xmtick(0.65(0.05)1.3) ylabel(1 2 4 8, angle(horiz) format(%2.1fc)) ytitle("Hazard Ratio of Mortality") xtitle("PA:AA") title("HR for mortality; Age 65+; adjusted for age(decade) and sex") yline(1, lp("shortdash") lc(gs10)) xline(0.9, lp("shortdash_dot") lc(gs10)) note("Sample mean (PA:AA = 0.85) taken as reference. Restricted cubic splines with 4 knots as described in Harrell (2001)" "Lines added at traditional threshold of abnormal (0.9) and HR = 1. Dashed line = 95% Confidence interval")
graph export "Results and Figures/$S_DATE/HR PAAA Splines - Age 65plus Adj AgeSex.png", as(png) name("Graph") replace
restore
