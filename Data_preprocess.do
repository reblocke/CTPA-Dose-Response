cd "/Users/blocke/Box Sync/Residency Personal Files/Scholarly Work/Locke Research Projects/Pulm Artery Stuff w Scarps/Local Analysis"
//cd "/Users/reblocke/Research/CTPA-Dose-Response"


capture mkdir "Results and Figures"
capture mkdir "Results and Figures/$S_DATE/" //make new folder for figure output if needed
/*capture mkdir "Results and Figures/$S_DATE/Logs/" //new folder for stata logs
local a1=substr(c(current_time),1,2)
local a2=substr(c(current_time),4,2)
local a3=substr(c(current_time),7,2)
local b = "Data_preprocess.do" // do file name
copy "`b'" "Results and Figures/$S_DATE/Logs/(`a1'_`a2'_`a3')`b'"
*/ 

use final_noempi, clear

set scheme cleanplots //cleanplots white_tableau white_w3d //3 options I usually use

missings dropobs, force //drop all rows with no observations. None
missings dropvars, force //drop all columns with no observations
missings report
count // n=990
count if lastfollowupday==0 // n=78 without any follow-up data

/* -----
Data Cleaning 
--------*/ 

label variable age "Age (years)"

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

//traditional cutpoints
label variable enlargedratio "PA:AA Ratio"
label define enlargedratio_lab 0 "Normal PA:AA (<0.9)" 1 "Increased PA:AA (0.9+)"
label values enlargedratio enlargedratio_lab

label variable enlargedpa "PA diameter"
label define enlargedpa_lab 0 "Normal PAd" 1 "Enlarged PAd"
label values enlargedpa enlargedpa_lab

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

//label comorbdities
label variable pulmdisease "Pulmonary Disease"
label variable chf "Congestive Heart Failure"
label variable diabetes3 "Diabetes"
label define diabetes_val 0 "No DM" 1 "Uncomplicated DM" 2 "DM with Complication(s)"
label values diabetes3 diabetes_val 
label variable hypertension "Hypertension"
label variable pulmonarycircdisorder "Pulm Circ. Disorder"
label variable peripheralvascdisorders "Periph Vasc. Disease"
label variable renalfailure "Kidney Disease"

// Split diabetes into complicated vs not. 
tab diabetes3
capture drop diabetes1
gen diabetes1=cond(diabetes3==1,1,0)
tab diabetes3 diabetes1, missing
capture drop diabetes2
gen diabetes2=cond(diabetes3==2,1,0)
tab diabetes3 diabetes2, missing
label variable diabetes1 "Uncomplicated Diabetes"
label variable diabetes2 "Diabetes with Complication(s)"

// Outcomes
label variable death "Death"
label define death_label 0 "Alive or Censored" 1 "Died"
label values death death_label
label variable lastfollowupyear "Follow-up or death(years)"

label variable anyemergency "Any Emergency Visits in Follow-up?"
label variable anyadmission "Any Hospital Admission in Follow-up?"

gen time_of_death = lastfollowupyear if death == 1
label variable time_of_death "Time of Death"
gen time_of_censoring = lastfollowupyear if death != 1
label variable time_of_censoring "Time of longest follow-up alive"

//Categorizations of continuous variables; for graphs
recode age min/30=0 30/40=1 40/50=2 50/60=3 60/70=4 70/max=5, gen(age_decade)
label define age_dec_lab 0 "<30 years" 1 "30-40 years" 2 "40-50 years" 3 "50-60 years" 4 "60-70 years" 5 "70+ years"
label variable age_decade "Age (by decade)"
label values age_decade age_dec_lab 



/*. 
Should we not use these? 
*/ 

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


/* Tertiles for EDA */ 

xtile mpad_tertile = mpad, n(3)  // Divides into 3 equal-sized groups
label define mpad_tertile_lbl 1 "MPAD Tertile 1" 2 "MPAD Tertile 2" 3 "MPAD Tertile 3"
label values mpad_tertile mpad_tertile_lbl
tab mpad_tertile

xtile aa_tertile = ascendingaorta, n(3)  // Divides into 3 equal-sized groups
label define aa_tertile_lbl 1 "AA Tertile 1" 2 "AA Tertile 2" 3 "AA Tertile 3"
label values aa_tertile aa_tertile_lbl
tab aa_tertile

xtile mpaaa_tertile = mpaaa, n(3)  // Divides into 3 equal-sized groups
label define mpaaa_tertile_lbl 1 "MPAAA Tertile 1" 2 "MPAAA Tertile 2" 3 "MPAAA Tertile 3"
label values mpaaa_tertile mpaaa_tertile_lbl
tab mpaaa_tertile








/******************************************************************
 Merge admin CSV (ACCT/EMPI/BMI/Weight/Height) into cleaned_noempi
 - run this after:  use cleaned_noempi
******************************************************************/

*-- Housekeeping: folders already created above in your script.
*   We'll reuse S_DATE and the HH_MM_SS tickers for file names.


local a1 = substr(c(current_time),1,2)
local a2 = substr(c(current_time),4,2)
local a3 = substr(c(current_time),7,2)


*------------------------------------------------------------------
* 1) MASTER: build a str# key for acct_no
*------------------------------------------------------------------
capture drop acct_key
capture confirm numeric variable acct_no
if !_rc {
    gen str80 acct_key = strtrim(stritrim(string(acct_no,"%20.0f")))
}
else {
    gen str80 acct_key = strtrim(stritrim(acct_no))
}
label var acct_key "Account number (merge key, str#)"
order acct_key, after(acct_no)

* (Optional) check for duplicates in master
quietly duplicates report acct_key
if r(N)>0 {
    di as txt "Note: master has duplicate acct_key values. 1:1 merge will fail; switch to m:1 if intended."
}

tempfile MASTER
save `MASTER', replace

*------------------------------------------------------------------
* 2) USING CSV: import and build the same str# key
*------------------------------------------------------------------
import delimited using "ACCT EMPI BMI Weight Height.csv", ///
    varnames(1) case(lower) clear

* Build acct_key to match master semantics
capture drop acct_key
capture confirm numeric variable acct_no
if !_rc {
    gen str80 acct_key = strtrim(stritrim(string(acct_no,"%20.0f")))
}
else {
    gen str80 acct_key = strtrim(stritrim(acct_no))
}
label var acct_key "Account number (merge key, str#)"

* Standardize metric columns and types
rename age  age_csv
rename bmi  bmi_csv
destring age_csv bmi_csv, replace ignore(", ")

* Timestamps for dedup preference
gen double admit_tc  = clock(admit_dts,  "YMDhms")
replace    admit_tc  = clock(admit_dts,  "MDYhms") if missing(admit_tc) & !missing(admit_dts)
format %tcDDmonCCYY_HH:MM admit_tc
gen byte miss_admit  = missing(admit_tc)

gen double weight_tc = clock(weight_dts, "YMDhms")
replace    weight_tc = clock(weight_dts, "MDYhms") if missing(weight_tc) & !missing(weight_dts)
format %tcDDmonCCYY_HH:MM weight_tc
gen byte miss_weight = missing(weight_tc)

* One row per acct_key from CSV:
gsort acct_key +miss_admit -admit_tc +miss_weight -weight_tc
by acct_key: keep if _n==1
drop miss_admit miss_weight

tempfile CSV
save `CSV', replace

*------------------------------------------------------------------
* 3) Merge and audit; CSV values take priority
*------------------------------------------------------------------
use `MASTER', clear

* Preserve originals for audit
capture confirm variable age
if !_rc {
    rename age age_old
    label var age_old "Age (years) — original"
}
else {
    gen age_old = .
    label var age_old "Age (years) — original (missing)"
}
capture confirm variable bmi
if !_rc {
    rename bmi bmi_old
    label var bmi_old "BMI — original"
}
else {
    gen bmi_old = .
    label var bmi_old "BMI — original (missing)"
}

* Use 1:1 if master keys are unique; if duplicates exist in master, change to m:1
merge 1:1 acct_key using `CSV', ///
    keep(1 3) ///
    keepusing(age_csv bmi_csv person_mk empi admit_dts fcilty_id weight height height_dts weight_dts) ///
    nogen
	
* Mismatch flags
gen byte age_mismatch = !missing(age_old, age_csv) & (age_old != age_csv)
gen byte bmi_mismatch = !missing(bmi_old, bmi_csv) & (abs(bmi_old - bmi_csv) > 0.01)
label define yesno 0 "No" 1 "Yes"
label values age_mismatch yesno
label values bmi_mismatch yesno
label var age_mismatch "Age differs (original vs CSV)"
label var bmi_mismatch "BMI differs (>0.01 units)"

* Signed differences and final values (CSV priority)
gen age_diff = age_csv - age_old if !missing(age_csv, age_old)
gen bmi_diff = bmi_csv - bmi_old if !missing(bmi_csv, bmi_old)
gen age = cond(!missing(age_csv), age_csv, age_old)
gen bmi = cond(!missing(bmi_csv), bmi_csv, bmi_old)
label var age "Age (years)"
label var bmi "BMI (kg/m^2)"

gen byte age_from_csv = !missing(age_csv)
gen byte bmi_from_csv = !missing(bmi_csv)
label values age_from_csv yesno
label values bmi_from_csv yesno
label var age_from_csv "Age sourced from CSV?"
label var bmi_from_csv "BMI sourced from CSV?"

* Counts and export of mismatches
quietly count if age_mismatch
local n_age_mis = r(N)
quietly count if bmi_mismatch
local n_bmi_mis = r(N)
di as result "MISMATCH COUNTS  ->  Age: `n_age_mis'   BMI: `n_bmi_mis'"

preserve
keep if age_mismatch | bmi_mismatch
keep acct_no acct_key age_old age_csv age_diff bmi_old bmi_csv bmi_diff male ctdate masterdate admit_dts weight height height_dts weight_dts empi person_mk fcilty_id
order acct_no acct_key age_old age_csv age_diff bmi_old bmi_csv bmi_diff, first
export delimited using "Results and Figures/$S_DATE/age_bmi_mismatches_(`a1'_`a2'_`a3').csv", replace
restore


/* Generate BSA information */

/******************************************************************
BSA from height & weight (robust units -> cm/kg)
Assumes variables: height, weight (may be string or numeric)
Creates: height_cm, weight_kg, bsa_mosteller, bsa_dubois, bsa, flags
******************************************************************/

label var height "Height (cm)"
label var weight "Weight (kg)"


* -- If time tags not defined earlier, define them (for QC export filenames)
capture confirm local a1
if _rc {
    local a1 = substr(c(current_time),1,2)
    local a2 = substr(c(current_time),4,2)
    local a3 = substr(c(current_time),7,2)
}

* Ensure numeric height/weight (handle stray commas/spaces)
capture confirm numeric variable height
if _rc destring height, replace ignore(", ")
capture confirm numeric variable weight
if _rc destring weight, replace ignore(", ")

* Standardize HEIGHT to centimeters, observation-wise
capture drop height_cm height_unit_guess
gen double height_cm = .
gen byte   height_unit_guess = .

* meters -> cm (typical adult 1.2–2.3 m)
replace height_cm = height*100 if inrange(height,1.2,2.3) & missing(height_cm)
replace height_unit_guess = 2   if inrange(height,1.2,2.3)

* inches -> cm (48–90 in)
replace height_cm = height*2.54 if inrange(height,48,90)   & missing(height_cm)
replace height_unit_guess = 3   if inrange(height,48,90)

* centimeters (120–250 cm)
replace height_cm = height      if inrange(height,120,250)  & missing(height_cm)
replace height_unit_guess = 1   if inrange(height,120,250)

* millimeters (rare; 1200–2500 mm)
replace height_cm = height/10   if inrange(height,1200,2500) & missing(height_cm)
replace height_unit_guess = 4   if inrange(height,1200,2500)

label define hu 1 "cm" 2 "m" 3 "in" 4 "mm"
label values height_unit_guess hu
label var height_cm "Height (cm, standardized)"

* Standardize WEIGHT to kilograms, observation-wise
capture drop weight_kg weight_unit_guess
gen double weight_kg = .
gen byte   weight_unit_guess = .

* pounds -> kg (66–700 lb)
replace weight_kg = weight*0.45359237 if inrange(weight,66,700)  & missing(weight_kg)
replace weight_unit_guess = 2         if inrange(weight,66,700)

* kilograms (30–300 kg)
replace weight_kg = weight            if inrange(weight,30,300)   & missing(weight_kg)
replace weight_unit_guess = 1         if inrange(weight,30,300)

* grams -> kg (10,000–300,000 g) [defensive]
replace weight_kg = weight/1000       if inrange(weight,10000,300000) & missing(weight_kg)
replace weight_unit_guess = 3         if inrange(weight,10000,300000)

label define wu 1 "kg" 2 "lb" 3 "g"
label values weight_unit_guess wu
label var weight_kg "Weight (kg, standardized)"

* Optional: quickly flag implausible standardized inputs
gen byte height_implaus = inrange(height_cm,.,.) & !inrange(height_cm,120,230)
gen byte weight_implaus = inrange(weight_kg,.,.) & !inrange(weight_kg,30,300)

* --- BSA formulas (metric) ---
capture drop bsa_mosteller bsa_dubois bsa_haycock bsa
gen double bsa_mosteller = sqrt((height_cm*weight_kg)/3600)          if !missing(height_cm, weight_kg)
gen double bsa_dubois    = 0.007184 * (height_cm^0.725) * (weight_kg^0.425) if !missing(height_cm, weight_kg)
* Optional pediatric-leaning alternative:
* gen double bsa_haycock   = 0.024265 * (height_cm^0.3964) * (weight_kg^0.5378) if !missing(height_cm, weight_kg)

* Final analysis variable: use Mosteller by default
gen double bsa = bsa_mosteller
label var bsa "BSA (m^2)"
label var bsa_mosteller "BSA (m^2) — Mosteller"
label var bsa_dubois    "BSA (m^2) — Du Bois & Du Bois"
* label var bsa_haycock   "BSA (m^2) — Haycock"


format %4.2f bsa bsa_mosteller bsa_dubois

* QC: out-of-range BSA (defensive adult cutpoints 1–3 m^2)
gen byte bsa_implaus = !missing(bsa) & (bsa<1 | bsa>3)
label define yn 0 "No" 1 "Yes"
capture label values bsa_implaus yn
label var bsa_implaus "BSA implausible (<1 or >3 m^2)"

* Summary on screen
di as result "BSA summary (Mosteller):"
quietly count if !missing(bsa)
di as txt "  Non-missing: " as res r(N)
quietly count if bsa_implaus
di as txt "  Implausible: " as res r(N)
sum height_cm weight_kg bsa

* Export QC table of suspect rows
preserve
keep if height_implaus | weight_implaus | bsa_implaus | missing(bsa)
order acct_no height weight height_cm weight_kg bsa bsa_mosteller bsa_dubois height_unit_guess weight_unit_guess
capture mkdir "Results and Figures"
capture mkdir "Results and Figures/$S_DATE"
export delimited using "Results and Figures/$S_DATE/bsa_qc_(`a1'_`a2'_`a3').csv", replace
restore





/* 
Data pre-processing 
Calculate Z-scores

The AA ones work OK, but the PA ones don't directly port to our population well. 
Thus - not used. 

  
 Nevsky (AJR 2011) is the only widely cited CT study in healthy adults that directly supports body‑size–adjusted z‑scores
 ajronline.org/doi/pdf/10.2214/ajr.10.4990
 
 
 
Wolak et al. AA -  outer‑to‑outer. 
 CT findings of
4,387 patients, age 26 to 92 years, free of known
clinical coronary heart disease (CHD), who underwent CAC scanning during the period from July
2004 to March 2007
 
 https://radiology-universe.org/calculator/ascending-aorta-diameter/11005.pdf
 
*/ 


*-----------------------------
* Z-scores: Wolak (AA) & Nevsky, Berger(PA)
* Requires: ascendingaorta (mm), mpad (mm), age (y), male (0/1), bsa (m^2)
*-----------------------------

*--- Guardrails on coding ---
assert inlist(male,0,1)

*--- Derive Wolak stratification factors ---
gen byte ageband = .
replace ageband = 1 if age < 45
replace ageband = 2 if inrange(age,45,54)
replace ageband = 3 if inrange(age,55,64)
replace ageband = 4 if age >= 65
label define ageband 1 "<45" 2 "45-54" 3 "55-64" 4 ">=65"
label values ageband ageband

gen byte bsaband = .
replace bsaband = 1 if bsa < 1.70
replace bsaband = 2 if inrange(bsa,1.70,1.89)
replace bsaband = 3 if inrange(bsa,1.90,2.09)
replace bsaband = 4 if bsa >= 2.10
label define bsaband 1 "<1.70" 2 "1.70-1.89" 3 "1.90-2.09" 4 ">=2.10"
label values bsaband bsaband

*--- Build Wolak normative table (ascending aorta, mm): mean & SD by sex×ageband×bsaband
preserve
tempfile wolak_norm
clear
input byte male byte ageband byte bsaband double(mean_aa sd_aa)

    // FEMALE (male==0)
    0 1 1 28.4 2.7
    0 1 2 30.0 2.2
    0 1 3 29.8 2.6
    0 1 4 31.3 .      // SD not reported in Wolak
    0 2 1 29.6 2.8
    0 2 2 31.4 2.9
    0 2 3 32.5 3.2
    0 2 4 34.4 3.1
    0 3 1 31.1 2.9
    0 3 2 31.8 2.6
    0 3 3 33.0 3.0
    0 3 4 35.4 3.3
    0 4 1 32.5 2.5
    0 4 2 33.4 2.9
    0 4 3 34.3 4.2
    0 4 4 32.8 .      // SD not reported

    // MALE (male==1)
    1 1 1 28.6 2.2
    1 1 2 30.1 3.1
    1 1 3 30.9 2.7
    1 1 4 32.3 3.0
    1 2 1 31.0 3.8
    1 2 2 31.7 3.2
    1 2 3 33.1 3.3
    1 2 4 34.4 3.1
    1 3 1 31.5 2.4
    1 3 2 33.5 3.1
    1 3 3 34.6 3.3
    1 3 4 36.1 3.5
    1 4 1 33.9 2.3
    1 4 2 35.0 3.0
    1 4 3 35.8 3.2
    1 4 4 36.8 2.8
end
label var mean_aa "Wolak mean AscAo (mm)"
label var sd_aa   "Wolak SD AscAo (mm)"
save `wolak_norm', replace
restore

*--- Merge Wolak norms & compute z-score for Ascending Aorta
merge m:1 male ageband bsaband using `wolak_norm', keep(master match) nogen

gen double z_aa_wolak = (ascendingaorta - mean_aa) / sd_aa if !missing(mean_aa, sd_aa, ascendingaorta)
label var z_aa_wolak "Z(Ascending Aorta) per Wolak, sex×age×BSA"

* Optional: if you want to impute missing Wolak SDs for the two female cells with missing SD,
* uncomment the next block to fill sd_aa with the mean SD within the same sex×ageband
/*
bys male ageband: egen double sd_fallback = mean(sd_aa)
replace sd_aa = sd_fallback if missing(sd_aa)
replace z_aa_wolak = (ascendingaorta - mean_aa) / sd_aa if missing(z_aa_wolak) & !missing(sd_aa,mean_aa,ascendingaorta)
drop sd_fallback
*/

*--- Nevsky: BSA-normalized PA diameter (cm/m^2) and z-score
gen double mpad_cm     = mpad/10
gen double mpad_norm   = mpad_cm / bsa     // units: cm per m^2

gen double nev_mean = .
gen double nev_sd   = .
replace nev_mean = 1.32 if male==1   // men
replace nev_sd   = 0.12 if male==1
replace nev_mean = 1.36 if male==0   // women
replace nev_sd   = 0.14 if male==0

gen double z_pa_nevsky = (mpad_norm - nev_mean) / nev_sd if !missing(mpad_norm, nev_mean, nev_sd)
label var z_pa_nevsky "Z(PA, cm/m^2) per Nevsky (sex-specific)"


* --- Berger-based expected mean and z-scores for mPA (ICVTS 2021/2022) ---

* Age groups to match Berger formulae
gen byte agegrp_berger = .
replace agegrp_berger = 1 if age < 45
replace agegrp_berger = 2 if age >= 45 & age < 55
replace agegrp_berger = 3 if age >= 55 & age < 65
replace agegrp_berger = 4 if age >= 65
label define agegrp_berger 1 "<45" 2 "45-54" 3 "55-64" 4 ">=65"
label values agegrp_berger agegrp_berger

* Predicted mean mPA by age group and BSA (Berger piecewise linear models)
gen double berg_mean = .
replace berg_mean = 4.31*bsa + 22.56 if agegrp_berger==1
replace berg_mean = 4.70*bsa + 23.37 if agegrp_berger==2
replace berg_mean = 5.71*bsa + 20.43 if agegrp_berger==3
replace berg_mean = 4.60*bsa + 24.27 if agegrp_berger==4
label var berg_mean "Berger mean mPA (mm) given age-group & BSA"

* BSA bins (to map stratum-specific SDs from Table 4)
gen byte bsa_bin_berger = .
replace bsa_bin_berger = 1 if bsa < 1.70
replace bsa_bin_berger = 2 if bsa >= 1.70 & bsa < 1.90
replace bsa_bin_berger = 3 if bsa >= 1.90 & bsa < 2.10
replace bsa_bin_berger = 4 if bsa >= 2.10
label define bsa_bin_berger 1 "<1.7" 2 "1.7–<1.9" 3 "1.9–<2.1" 4 ">=2.1"
label values bsa_bin_berger bsa_bin_berger

* SD lookup from Berger Table 4 (age x BSA)
gen double berg_sd = .

* <45 y
replace berg_sd = 3.67 if agegrp_berger==1 & bsa_bin_berger==1
replace berg_sd = 3.58 if agegrp_berger==1 & bsa_bin_berger==2
replace berg_sd = 4.27 if agegrp_berger==1 & bsa_bin_berger==3
replace berg_sd = 3.86 if agegrp_berger==1 & bsa_bin_berger==4

* 45–54 y
replace berg_sd = 4.84 if agegrp_berger==2 & bsa_bin_berger==1
replace berg_sd = 4.49 if agegrp_berger==2 & bsa_bin_berger==2
replace berg_sd = 4.35 if agegrp_berger==2 & bsa_bin_berger==3
replace berg_sd = 4.25 if agegrp_berger==2 & bsa_bin_berger==4

* 55–64 y
replace berg_sd = 3.43 if agegrp_berger==3 & bsa_bin_berger==1
replace berg_sd = 4.42 if agegrp_berger==3 & bsa_bin_berger==2
replace berg_sd = 3.71 if agegrp_berger==3 & bsa_bin_berger==3
replace berg_sd = 4.33 if agegrp_berger==3 & bsa_bin_berger==4

* >=65 y
replace berg_sd = 4.40 if agegrp_berger==4 & bsa_bin_berger==1
replace berg_sd = 5.43 if agegrp_berger==4 & bsa_bin_berger==2
replace berg_sd = 5.00 if agegrp_berger==4 & bsa_bin_berger==3
replace berg_sd = 6.22 if agegrp_berger==4 & bsa_bin_berger==4

label var berg_sd "Berger SD (mm) for mPA by age-group & BSA bin"

* Berger z-score (preferred: stratum-specific SD)
gen double z_pa_berger = (mpad - berg_mean) / berg_sd
label var z_pa_berger "z-score mPA vs Berger (age+BSA-specific SD)"

* Optional: global-SD variant (uses overall SD=4.6 mm reported by Berger)
gen double z_pa_berger_globalSD = (mpad - berg_mean) / 4.6
label var z_pa_berger_globalSD "z-score mPA vs Berger (global SD=4.6)"

* --- Sanity check summary (add Berger vars to your existing display) ---
display "—  Sanity check  —"
summ z_aa_wolak z_pa_nevsky z_pa_berger mean_aa sd_aa nev_mean nev_sd berg_mean berg_sd, detail






mdesc
codebook
save cleaned_noempi, replace
di as result "Saved merged dataset: cleaned_noempi.dta"



