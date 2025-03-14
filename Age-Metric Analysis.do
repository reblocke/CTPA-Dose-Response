cd "/Users/blocke/Box Sync/Residency Personal Files/Scholarly Work/Locke Research Projects/Pulm Artery Stuff w Scarps/Local Analysis"

capture mkdir "Results and Figures"
capture mkdir "Results and Figures/$S_DATE/" //make new folder for figure output if needed
capture mkdir "Results and Figures/$S_DATE/Logs/" //new folder for stata logs
local a1=substr(c(current_time),1,2)
local a2=substr(c(current_time),4,2)
local a3=substr(c(current_time),7,2)
local b = "BL CombofilewithoutEMPI.do" // do file name
copy "`b'" "Results and Figures/$S_DATE/Logs/(`a1'_`a2'_`a3')`b'"

use final_noempi, clear

set scheme cleanplots //cleanplots white_tableau white_w3d //3 options I usually use

missings dropobs, force //drop all rows with no observations. None
missings dropvars, force //drop all columns with no observations
missings report
mdesc
codebook
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

// [ ] remove? 
gen known_assoc_comorb = 0
replace known_assoc_comorb = 1 if (pulmdisease == 1 | pulmonarycircdisorders == 1 | chf == 1) 
tab chf known_assoc_comorb, missing //sanity checks
tab pulmdisease known_assoc_comorb, missing
tab pulmonarycircdisorders known_assoc_comorb, missing
label variable known_assoc_comorb "Known Comorbidity w PA->Mortality Association?"
label define known_assoc_comorb_lab 0 "No Pulm Dz/CHF/Pulm Vasc Dz" 1 "Any of Pulm Dz/CHF/Pulm Vasc Dz"
label values known_assoc_comorb known_assoc_comorb_lab

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

save cleaned_noempi, replace

clear
use cleaned_noempi

/* Analysis */ 

/* 
Table 1
*/ 

table1_mc, by (age_decade) ///
vars( ///
male bin %4.1f \ /// 
obesity_calc bin %4.1f \ /// 
pulmdisease bin %4.1f \ /// 
chf bin %4.1f \ /// 
diabetes3 cat %4.1f \ /// 
hypertension bin %4.1f \ ///
pulmonarycircdisorders bin %4.1f \ /// 
peripheralvascdisorders bin %4.1f \ /// 
renalfailure bin %4.1f \ /// 
mpad conts %4.1f \ ///
enlargedpa cat %4.1f \ ///
ascendingaorta conts %4.1f \ ///
mpaaa conts %4.2f \ ///
enlargedratio cat %4.1f \ /// 
anyemergency bin %4.1f \ /// 
anyadmission bin %4.1f \ /// 
death bin %4.1f \ /// 
time_of_death conts %4.1f \ /// 
time_of_death conts %4.1f \ /// 
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol total(before) ///
saving("Results and Figures/$S_DATE/Table 1 PA enlargement by Age.xlsx", replace)


/* Raw Data */ 
colorpalette Spectral, n(6) nograph  // [ ] TODO: use same CET C6 pallette somehow?
local color6 `"`r(p1)'"' // Reverse order
local color5 `"`r(p2)'"'
local color4 `"`r(p3)'"'
local color3 `"`r(p4)'"'
local color2 `"`r(p5)'"'
local color1 `"`r(p6)'"'

twoway ///
    (scatter ascendingaorta mpad if age_decade == 0 & male == 0, mcolor("`color1'") msymbol(O)) || /// <30 years, Female
    (scatter ascendingaorta mpad if age_decade == 0 & male == 1, mcolor("`color1'") msymbol(S)) || /// <30 years, Male
    (scatter ascendingaorta mpad if age_decade == 1 & male == 0, mcolor("`color2'") msymbol(O)) || /// 30-40 years, Female
    (scatter ascendingaorta mpad if age_decade == 1 & male == 1, mcolor("`color2'") msymbol(S)) || /// 30-40 years, Male
    (scatter ascendingaorta mpad if age_decade == 2 & male == 0, mcolor("`color3'") msymbol(O)) || /// 40-50 years, Female
    (scatter ascendingaorta mpad if age_decade == 2 & male == 1, mcolor("`color3'") msymbol(S)) || /// 40-50 years, Male
    (scatter ascendingaorta mpad if age_decade == 3 & male == 0, mcolor("`color4'") msymbol(O)) || /// 50-60 years, Female
    (scatter ascendingaorta mpad if age_decade == 3 & male == 1, mcolor("`color4'") msymbol(S)) || /// 50-60 years, Male
    (scatter ascendingaorta mpad if age_decade == 4 & male == 0, mcolor("`color5'") msymbol(O)) || /// 60-70 years, Female
    (scatter ascendingaorta mpad if age_decade == 4 & male == 1, mcolor("`color5'") msymbol(S)) || /// 60-70 years, Male
    (scatter ascendingaorta mpad if age_decade == 5 & male == 0, mcolor("`color6'") msymbol(O)) || /// 70+ years, Female
    (scatter ascendingaorta mpad if age_decade == 5 & male == 1, mcolor("`color6'") msymbol(S)), /// 70+ years, Male
    ///
    legend(order(1 "<30 Female" 2 "<30 Male" 3 "30-40 Female" 4 "30-40 Male" ///
                 5 "40-50 Female" 6 "40-50 Male" 7 "50-60 Female" 8 "50-60 Male" ///
                 9 "60-70 Female" 10 "60-70 Male" 11 "70+ Female" 12 "70+ Male")) ///
    xtitle("Pulmonary Artery Diameter (mm)") ytitle("Ascending Aorta Diameter (mm)") ///
    title("Scatterplot of PA vs Ascending Aorta by Age Group & Sex") ///
    xlabel(, labsize(medlarge)) ylabel(, labsize(medlarge)) ///
    scheme(white_w3d)
	



//ridgelines plots to visualize overall change by age. 

/* PAd */ 
ridgeline mpad, by(age_decade) yline ylw(0.2) overlap(1.7) ylc(blue) ylp(dot) ////
	labpos(right) bwid(1.2) laboffset(-3) showstats xlabel(15(5)40) ///
	palette(CET C6) alpha(75) xtitle("Pulmonary Artery Diameter (mm)")
	
* Ridgeline plot for Females (male == 0)
ridgeline mpad if male == 0, by(age_decade) yline ylw(0.2) overlap(1.7) ///
    ylc(blue) ylp(dot) labpos(right) bwid(1.2) laboffset(-3) ///
    showstats xlabel(15(5)40) palette(CET C6) alpha(75) ///
    xtitle("Pulmonary Artery Diameter (mm)") ///
    title("PA Diameter by Age Group - Female") ///
    saving(ridgeline_female, replace)

* Ridgeline plot for Males (male == 1)
ridgeline mpad if male == 1, by(age_decade) yline ylw(0.2) overlap(1.7) ///
    ylc(red) ylp(dot) labpos(right) bwid(1.2) laboffset(-3) ///
    showstats xlabel(15(5)40) palette(CET C6) alpha(75) ///
    xtitle("Pulmonary Artery Diameter (mm)") ///
    title("PA Diameter by Age Group - Male") ///
    saving(ridgeline_male, replace)

graph combine ridgeline_female.gph ridgeline_male.gph, ///
    title("Ridgeline Plot of PA Diameter by Age Group & Sex")
	

/* Asc Aorta */ 
ridgeline ascendingaorta, by(age_decade) yline ylw(0.2) overlap(1.7) ylc(blue) ylp(dot) ////
	labpos(right) bwid(1.2) laboffset(-3) showstats xlabel(15(5)45) ///
	palette(CET C6) alpha(75) xtitle("Ascending Aorta Diameter (mm)")

* Ridgeline plot for Females (male == 0) - Ascending Aorta
ridgeline ascendingaorta if male == 0, by(age_decade) yline ylw(0.2) overlap(1.7) ///
    ylc(blue) ylp(dot) labpos(right) bwid(1.2) laboffset(-3) ///
    showstats xlabel(15(5)45) palette(CET C6) alpha(75) ///
    xtitle("Ascending Aorta Diameter (mm)") ///
    title("Ascending Aorta by Age Group - Female") ///
    saving(ridgeline_aorta_female, replace)

* Ridgeline plot for Males (male == 1) - Ascending Aorta
ridgeline ascendingaorta if male == 1, by(age_decade) yline ylw(0.2) overlap(1.7) ///
    ylc(red) ylp(dot) labpos(right) bwid(1.2) laboffset(-3) ///
    showstats xlabel(15(5)45) palette(CET C6) alpha(75) ///
    xtitle("Ascending Aorta Diameter (mm)") ///
    title("Ascending Aorta by Age Group - Male") ///
    saving(ridgeline_aorta_male, replace)

graph combine ridgeline_aorta_female.gph ridgeline_aorta_male.gph, ///
    title("Ridgeline Plot of Ascending Aorta by Age & Sex")
	
	
/* MPA:AA */ 
ridgeline mpaaa, by(age_decade) yline ylw(0.2) overlap(1.7) ylc(blue) ylp(dot) ////
	labpos(right) bwid(0.05) laboffset(-0.05) showstats xlabel(0.6(.1)1.3) ///
	palette(CET C6) alpha(75) xtitle("PA:AA Ratio")

* Ridgeline plot for Females (male == 0) - PA:AA Ratio
ridgeline mpaaa if male == 0, by(age_decade) yline ylw(0.2) overlap(1.7) ///
    ylc(blue) ylp(dot) labpos(right) bwid(0.05) laboffset(-0.05) ///
    showstats xlabel(0.6(.1)1.3) palette(CET C6) alpha(75) ///
    xtitle("PA:AA Ratio") ///
    title("PA:AA Ratio by Age Group - Female") ///
    saving(ridgeline_mpaaa_female, replace)

* Ridgeline plot for Males (male == 1) - PA:AA Ratio
ridgeline mpaaa if male == 1, by(age_decade) yline ylw(0.2) overlap(1.7) ///
    ylc(red) ylp(dot) labpos(right) bwid(0.05) laboffset(-0.05) ///
    showstats xlabel(0.6(.1)1.3) palette(CET C6) alpha(75) ///
    xtitle("PA:AA Ratio") ///
    title("PA:AA Ratio by Age Group - Male") ///
    saving(ridgeline_mpaaa_male, replace)

graph combine ridgeline_mpaaa_female.gph ridgeline_mpaaa_male.gph, ///
    title("Ridgeline Plot of PA:AA Ratio by Age & Sex")
	


	
/* -------------------
Quantile Regression 
-------------------*/ 
	
/* PA size and age */ 

/* All individuals */ 
/* First - shared quantiles */ 

/* Run Bootstrap Quantile Regressions for the Full Population */
bsqreg mpad c.age, quantile(50) reps(500) // 50th Percentile
bsqreg mpad c.age, quantile(90) reps(500) // 90th Percentile
bsqreg mpad c.age, quantile(95) reps(500) // 95th Percentile

/* Extract Colors from Spectral Palette */
colorpalette Spectral, n(6) nograph  
local color6 `"`r(p1)'"' // Reverse order
local color5 `"`r(p2)'"'
local color4 `"`r(p3)'"'
local color3 `"`r(p4)'"'
local color2 `"`r(p5)'"'
local color1 `"`r(p6)'"'

/* Run Quantile Regressions for the Full Population */
qreg mpad c.age, quantile(50) // 50th Percentile
predict mpad_q50_all, xb

qreg mpad c.age, quantile(90) // 90th Percentile
predict mpad_q90_all, xb

qreg mpad c.age, quantile(95) // 95th Percentile
predict mpad_q95_all, xb

/* Scatter Plot with Both Genders + Quantile Regression Lines */
twoway ///
    (scatter mpad age if age_decade == 0 & male == 0, mcolor("`color1'") msymbol(O)) || /// <30 years, Female
    (scatter mpad age if age_decade == 0 & male == 1, mcolor("`color1'") msymbol(S)) || /// <30 years, Male
    (scatter mpad age if age_decade == 1 & male == 0, mcolor("`color2'") msymbol(O)) || /// 30-40 years, Female
    (scatter mpad age if age_decade == 1 & male == 1, mcolor("`color2'") msymbol(S)) || /// 30-40 years, Male
    (scatter mpad age if age_decade == 2 & male == 0, mcolor("`color3'") msymbol(O)) || /// 40-50 years, Female
    (scatter mpad age if age_decade == 2 & male == 1, mcolor("`color3'") msymbol(S)) || /// 40-50 years, Male
    (scatter mpad age if age_decade == 3 & male == 0, mcolor("`color4'") msymbol(O)) || /// 50-60 years, Female
    (scatter mpad age if age_decade == 3 & male == 1, mcolor("`color4'") msymbol(S)) || /// 50-60 years, Male
    (scatter mpad age if age_decade == 4 & male == 0, mcolor("`color5'") msymbol(O)) || /// 60-70 years, Female
    (scatter mpad age if age_decade == 4 & male == 1, mcolor("`color5'") msymbol(S)) || /// 60-70 years, Male
    (scatter mpad age if age_decade == 5 & male == 0, mcolor("`color6'") msymbol(O)) || /// 70+ years, Female
    (scatter mpad age if age_decade == 5 & male == 1, mcolor("`color6'") msymbol(S)) || /// 70+ years, Male
    (line mpad_q50_all age, lcolor(black) lpattern(solid) lwidth(medthick)) || /// 50th quantile
    (line mpad_q90_all age, lcolor(gs8) lpattern(solid) lwidth(medthick)) || /// 90th quantile
    (line mpad_q95_all age, lcolor(gs12) lpattern(solid) lwidth(medthick)), /// 95th quantile
    legend(order(1 "<30 Female" 2 "<30 Male" 3 "30-40 Female" 4 "30-40 Male" ///
                 5 "40-50 Female" 6 "40-50 Male" 7 "50-60 Female" 8 "50-60 Male" ///
                 9 "60-70 Female" 10 "60-70 Male" 11 "70+ Female" 12 "70+ Male" ///
                 13 "50th Quantile" 14 "90th Quantile" 15 "95th Quantile")) ///
    xtitle("Age (years)") ytitle("Pulmonary Artery Diameter (mm)") ///
    title("Quantile Regression of MPAD by Age (Both Sexes)") ///
    xlabel(, labsize(medlarge)) ylabel(, labsize(medlarge)) ///
    scheme(white_w3d)


/* -----
Sex-Stratified Quantile Regression Scatter Plots
-----*/ 

/* Run Bootstrap Quantile Regressions for the Men */
bsqreg mpad c.age if male == 1, quantile(50) reps(500) // 50th Percentile
bsqreg mpad c.age if male == 1, quantile(90) reps(500) // 90th Percentile
bsqreg mpad c.age if male == 1, quantile(95) reps(500) // 95th Percentile

/* Extract Colors from Spectral Palette */
colorpalette Spectral, n(6) nograph  
local color6 `"`r(p1)'"' // Reverse order
local color5 `"`r(p2)'"'
local color4 `"`r(p3)'"'
local color3 `"`r(p4)'"'
local color2 `"`r(p5)'"'
local color1 `"`r(p6)'"'

/* -----
Generate Quantile Regression Predictions for Males
-----*/ 
qreg mpad c.age if male == 1, quantile(50) // Male
predict mpad_q50_male if male == 1, xb
qreg mpad c.age if male == 1, quantile(90) // Male
predict mpad_q90_male if male == 1, xb
qreg mpad c.age if male == 1, quantile(95) // Male
predict mpad_q95_male if male == 1, xb

/* -----
Male Scatter Plot
-----*/ 
twoway ///
    (scatter mpad age if age_decade == 0 & male == 1, mcolor("`color1'") msymbol(S)) || /// <30 years
    (scatter mpad age if age_decade == 1 & male == 1, mcolor("`color2'") msymbol(S)) || /// 30-40 years
    (scatter mpad age if age_decade == 2 & male == 1, mcolor("`color3'") msymbol(S)) || /// 40-50 years
    (scatter mpad age if age_decade == 3 & male == 1, mcolor("`color4'") msymbol(S)) || /// 50-60 years
    (scatter mpad age if age_decade == 4 & male == 1, mcolor("`color5'") msymbol(S)) || /// 60-70 years
    (scatter mpad age if age_decade == 5 & male == 1, mcolor("`color6'") msymbol(S)) || /// 70+ years
    (line mpad_q50_male age, lcolor(black) lpattern(solid) lwidth(medthick)) || /// 50th percentile line
    (line mpad_q90_male age, lcolor(gs8) lpattern(solid) lwidth(medthick)) || /// 90th percentile line
    (line mpad_q95_male age, lcolor(gs12) lpattern(solid) lwidth(medthick)), /// 95th percentile line
    ///
    legend(order(7 "50th Quantile" 8 "90th Quantile" 9 "95th Quantile") pos(6) row(1)) /// Only quantile labels in legend
    xtitle("Age (years)") ytitle("Pulmonary Artery Diameter (mm)", margin(medlarge)) ///
    title("Men") ///
    xlabel(, labsize(medlarge)) ylabel(15(5)40, nogrid labsize(medlarge)) /// Adjusted Y-axis range
    scheme(white_w3d) yscale(range(15 40)) ///
    saving(male_qreg_plot, replace)


/* Run Bootstrap Quantile Regressions for the Women */
bsqreg mpad c.age if male == 0, quantile(50) reps(500) // 50th Percentile
bsqreg mpad c.age if male == 0, quantile(90) reps(500) // 90th Percentile
bsqreg mpad c.age if male == 0, quantile(95) reps(500) // 95th Percentile

/* -----
Generate Quantile Regression Predictions for Females
-----*/ 
qreg mpad c.age if male == 0, quantile(50) // Female
predict mpad_q50_female if male == 0, xb
qreg mpad c.age if male == 0, quantile(90) // Female
predict mpad_q90_female if male == 0, xb
qreg mpad c.age if male == 0, quantile(95) // Female
predict mpad_q95_female if male == 0, xb

/* -----
Female Scatter Plot
-----*/ 
twoway ///
    (scatter mpad age if age_decade == 0 & male == 0, mcolor("`color1'") msymbol(O)) || /// <30 years
    (scatter mpad age if age_decade == 1 & male == 0, mcolor("`color2'") msymbol(O)) || /// 30-40 years
    (scatter mpad age if age_decade == 2 & male == 0, mcolor("`color3'") msymbol(O)) || /// 40-50 years
    (scatter mpad age if age_decade == 3 & male == 0, mcolor("`color4'") msymbol(O)) || /// 50-60 years
    (scatter mpad age if age_decade == 4 & male == 0, mcolor("`color5'") msymbol(O)) || /// 60-70 years
    (scatter mpad age if age_decade == 5 & male == 0, mcolor("`color6'") msymbol(O)) || /// 70+ years
    (line mpad_q50_female age, lcolor(black) lpattern(solid) lwidth(medthick)) || /// 50th percentile line
    (line mpad_q90_female age, lcolor(gs8) lpattern(solid) lwidth(medthick)) || /// 90th percentile line
    (line mpad_q95_female age, lcolor(gs12) lpattern(solid) lwidth(medthick)), /// 95th percentile line
    ///
    legend(order(7 "50th Quantile" 8 "90th Quantile" 9 "95th Quantile") pos(6) row(1)) /// Only quantile labels in legend
    xtitle("Age (years)") ytitle("Pulmonary Artery Diameter (mm)", margin(medlarge)) ///
    title("Women") ///
    xlabel(, labsize(medlarge)) ylabel(15(5)40, nogrid labsize(medlarge)) /// Adjusted Y-axis range
    scheme(white_w3d) yscale(range(15 40)) ///
    saving(female_qreg_plot, replace)

/* -----
Combine the Two Plots into a Single Figure
-----*/ 
graph combine male_qreg_plot.gph female_qreg_plot.gph, ///
    title("Sex-Stratified Quantile Regression of MPAD by Age") ///
    ycommon xcommon ///
    iscale(*1.2) ///
    rows(1) ///
    graphregion(margin(2 2 2 2)) ///
    xsize(9) ysize(5)
	
	
	
/* Ascending Aorta */ 


/* -----
Quantile Regression Analysis for Ascending Aorta
-----*/

/* Run Bootstrap Quantile Regressions for the Full Population */
bsqreg ascendingaorta c.age, quantile(50) reps(500) // 50th Percentile
bsqreg ascendingaorta c.age, quantile(90) reps(500) // 90th Percentile
bsqreg ascendingaorta c.age, quantile(95) reps(500) // 95th Percentile

/* Extract Colors from Spectral Palette */
colorpalette Spectral, n(6) nograph  
local color6 `"`r(p1)'"' // Reverse order
local color5 `"`r(p2)'"'
local color4 `"`r(p3)'"'
local color3 `"`r(p4)'"'
local color2 `"`r(p5)'"'
local color1 `"`r(p6)'"'

/* Run Quantile Regressions for the Full Population */
qreg ascendingaorta c.age, quantile(50) // 50th Percentile
predict ascendingaorta_q50_all, xb

qreg ascendingaorta c.age, quantile(90) // 90th Percentile
predict ascendingaorta_q90_all, xb

qreg ascendingaorta c.age, quantile(95) // 95th Percentile
predict ascendingaorta_q95_all, xb

/* Scatter Plot with Both Genders + Quantile Regression Lines */
twoway ///
    (scatter ascendingaorta age if age_decade == 0 & male == 0, mcolor("`color1'") msymbol(O)) || /// <30 years, Female
    (scatter ascendingaorta age if age_decade == 0 & male == 1, mcolor("`color1'") msymbol(S)) || /// <30 years, Male
    (scatter ascendingaorta age if age_decade == 1 & male == 0, mcolor("`color2'") msymbol(O)) || /// 30-40 years, Female
    (scatter ascendingaorta age if age_decade == 1 & male == 1, mcolor("`color2'") msymbol(S)) || /// 30-40 years, Male
    (scatter ascendingaorta age if age_decade == 2 & male == 0, mcolor("`color3'") msymbol(O)) || /// 40-50 years, Female
    (scatter ascendingaorta age if age_decade == 2 & male == 1, mcolor("`color3'") msymbol(S)) || /// 40-50 years, Male
    (scatter ascendingaorta age if age_decade == 3 & male == 0, mcolor("`color4'") msymbol(O)) || /// 50-60 years, Female
    (scatter ascendingaorta age if age_decade == 3 & male == 1, mcolor("`color4'") msymbol(S)) || /// 50-60 years, Male
    (scatter ascendingaorta age if age_decade == 4 & male == 0, mcolor("`color5'") msymbol(O)) || /// 60-70 years, Female
    (scatter ascendingaorta age if age_decade == 4 & male == 1, mcolor("`color5'") msymbol(S)) || /// 60-70 years, Male
    (scatter ascendingaorta age if age_decade == 5 & male == 0, mcolor("`color6'") msymbol(O)) || /// 70+ years, Female
    (scatter ascendingaorta age if age_decade == 5 & male == 1, mcolor("`color6'") msymbol(S)) || /// 70+ years, Male
    (line ascendingaorta_q50_all age, lcolor(black) lpattern(solid) lwidth(medthick)) || /// 50th quantile
    (line ascendingaorta_q90_all age, lcolor(gs8) lpattern(solid) lwidth(medthick)) || /// 90th quantile
    (line ascendingaorta_q95_all age, lcolor(gs12) lpattern(solid) lwidth(medthick)), /// 95th quantile
    legend(order(13 "50th Quantile" 14 "90th Quantile" 15 "95th Quantile")) ///
    xtitle("Age (years)") ytitle("Ascending Aorta Diameter (mm)") ///
    title("Quantile Regression of Ascending Aorta by Age (Both Sexes)") ///
    xlabel(, labsize(medlarge)) ylabel(15(5)50, nogrid labsize(medlarge)) ///
    scheme(white_w3d)
	
	
	
	
/* -----
Sex-Stratified Quantile Regression Scatter Plots for Ascending Aorta
-----*/ 

/* Run Bootstrap Quantile Regressions for the Men */
bsqreg ascendingaorta c.age if male == 1, quantile(50) reps(500) // 50th Percentile
bsqreg ascendingaorta c.age if male == 1, quantile(90) reps(500) // 90th Percentile
bsqreg ascendingaorta c.age if male == 1, quantile(95) reps(500) // 95th Percentile

/* Extract Colors from Spectral Palette */
colorpalette Spectral, n(6) nograph  
local color6 `"`r(p1)'"' // Reverse order
local color5 `"`r(p2)'"'
local color4 `"`r(p3)'"'
local color3 `"`r(p4)'"'
local color2 `"`r(p5)'"'
local color1 `"`r(p6)'"'

/* -----
Generate Quantile Regression Predictions for Males
-----*/ 
qreg ascendingaorta c.age if male == 1, quantile(50) // Male
predict ascendingaorta_q50_male if male == 1, xb
qreg ascendingaorta c.age if male == 1, quantile(90) // Male
predict ascendingaorta_q90_male if male == 1, xb
qreg ascendingaorta c.age if male == 1, quantile(95) // Male
predict ascendingaorta_q95_male if male == 1, xb

/* -----
Male Scatter Plot
-----*/ 
twoway ///
    (scatter ascendingaorta age if age_decade == 0 & male == 1, mcolor("`color1'") msymbol(S)) || /// <30 years
    (scatter ascendingaorta age if age_decade == 1 & male == 1, mcolor("`color2'") msymbol(S)) || /// 30-40 years
    (scatter ascendingaorta age if age_decade == 2 & male == 1, mcolor("`color3'") msymbol(S)) || /// 40-50 years
    (scatter ascendingaorta age if age_decade == 3 & male == 1, mcolor("`color4'") msymbol(S)) || /// 50-60 years
    (scatter ascendingaorta age if age_decade == 4 & male == 1, mcolor("`color5'") msymbol(S)) || /// 60-70 years
    (scatter ascendingaorta age if age_decade == 5 & male == 1, mcolor("`color6'") msymbol(S)) || /// 70+ years
    (line ascendingaorta_q50_male age, lcolor(black) lpattern(solid) lwidth(medthick)) || /// 50th percentile line
    (line ascendingaorta_q90_male age, lcolor(gs8) lpattern(solid) lwidth(medthick)) || /// 90th percentile line
    (line ascendingaorta_q95_male age, lcolor(gs12) lpattern(solid) lwidth(medthick)), /// 95th percentile line
    ///
    legend(order(7 "50th Quantile" 8 "90th Quantile" 9 "95th Quantile") pos(6) row(1)) /// Only quantile labels in legend
    xtitle("Age (years)") ytitle("Ascending Aorta Diameter (mm)", margin(medlarge)) ///
    title("Men") ///
    xlabel(, labsize(medlarge)) ylabel(15(5)50, nogrid labsize(medlarge)) /// Adjusted Y-axis range
    scheme(white_w3d) yscale(range(15 50)) ///
    saving(male_qreg_plot, replace)


/* Run Bootstrap Quantile Regressions for the Women */
bsqreg ascendingaorta c.age if male == 0, quantile(50) reps(500) // 50th Percentile
bsqreg ascendingaorta c.age if male == 0, quantile(90) reps(500) // 90th Percentile
bsqreg ascendingaorta c.age if male == 0, quantile(95) reps(500) // 95th Percentile

/* -----
Generate Quantile Regression Predictions for Females
-----*/ 
qreg ascendingaorta c.age if male == 0, quantile(50) // Female
predict ascendingaorta_q50_female if male == 0, xb
qreg ascendingaorta c.age if male == 0, quantile(90) // Female
predict ascendingaorta_q90_female if male == 0, xb
qreg ascendingaorta c.age if male == 0, quantile(95) // Female
predict ascendingaorta_q95_female if male == 0, xb

/* -----
Female Scatter Plot
-----*/ 
twoway ///
    (scatter ascendingaorta age if age_decade == 0 & male == 0, mcolor("`color1'") msymbol(O)) || /// <30 years
    (scatter ascendingaorta age if age_decade == 1 & male == 0, mcolor("`color2'") msymbol(O)) || /// 30-40 years
    (scatter ascendingaorta age if age_decade == 2 & male == 0, mcolor("`color3'") msymbol(O)) || /// 40-50 years
    (scatter ascendingaorta age if age_decade == 3 & male == 0, mcolor("`color4'") msymbol(O)) || /// 50-60 years
    (scatter ascendingaorta age if age_decade == 4 & male == 0, mcolor("`color5'") msymbol(O)) || /// 60-70 years
    (scatter ascendingaorta age if age_decade == 5 & male == 0, mcolor("`color6'") msymbol(O)) || /// 70+ years
    (line ascendingaorta_q50_female age, lcolor(black) lpattern(solid) lwidth(medthick)) || /// 50th percentile line
    (line ascendingaorta_q90_female age, lcolor(gs8) lpattern(solid) lwidth(medthick)) || /// 90th percentile line
    (line ascendingaorta_q95_female age, lcolor(gs12) lpattern(solid) lwidth(medthick)), /// 95th percentile line
    ///
    legend(order(7 "50th Quantile" 8 "90th Quantile" 9 "95th Quantile") pos(6) row(1)) /// Only quantile labels in legend
    xtitle("Age (years)") ytitle("Ascending Aorta Diameter (mm)", margin(medlarge)) ///
    title("Women") ///
    xlabel(, labsize(medlarge)) ylabel(15(5)50, nogrid labsize(medlarge)) /// Adjusted Y-axis range
    scheme(white_w3d) yscale(range(15 50)) ///
    saving(female_qreg_plot, replace)

/* -----
Combine the Two Plots into a Single Figure
-----*/ 
graph combine male_qreg_plot.gph female_qreg_plot.gph, ///
    title("Sex-Stratified Quantile Regression of Ascending Aorta by Age") ///
    ycommon xcommon ///
    iscale(*1.2) ///
    rows(1) ///
    graphregion(margin(2 2 2 2)) ///
    xsize(9) ysize(5)


	
	

/* -----
PA:AA Ratio
-----*/

/* Run Bootstrap Quantile Regressions for the Full Population */
bsqreg mpaaa c.age, quantile(50) reps(500) // 50th Percentile
bsqreg mpaaa c.age, quantile(90) reps(500) // 90th Percentile
bsqreg mpaaa c.age, quantile(95) reps(500) // 95th Percentile

/* Extract Colors from Spectral Palette */
colorpalette Spectral, n(6) nograph  
local color6 `"`r(p1)'"' // Reverse order
local color5 `"`r(p2)'"'
local color4 `"`r(p3)'"'
local color3 `"`r(p4)'"'
local color2 `"`r(p5)'"'
local color1 `"`r(p6)'"'

/* Run Quantile Regressions for the Full Population */
qreg mpaaa c.age, quantile(50) // 50th Percentile
predict mpaaa_q50_all, xb

qreg mpaaa c.age, quantile(90) // 90th Percentile
predict mpaaa_q90_all, xb

qreg mpaaa c.age, quantile(95) // 95th Percentile
predict mpaaa_q95_all, xb

/* Scatter Plot with Both Genders + Quantile Regression Lines */
twoway ///
    (scatter mpaaa age if age_decade == 0 & male == 0, mcolor("`color1'") msymbol(O)) || /// <30 years, Female
    (scatter mpaaa age if age_decade == 0 & male == 1, mcolor("`color1'") msymbol(S)) || /// <30 years, Male
    (scatter mpaaa age if age_decade == 1 & male == 0, mcolor("`color2'") msymbol(O)) || /// 30-40 years, Female
    (scatter mpaaa age if age_decade == 1 & male == 1, mcolor("`color2'") msymbol(S)) || /// 30-40 years, Male
    (scatter mpaaa age if age_decade == 2 & male == 0, mcolor("`color3'") msymbol(O)) || /// 40-50 years, Female
    (scatter mpaaa age if age_decade == 2 & male == 1, mcolor("`color3'") msymbol(S)) || /// 40-50 years, Male
    (scatter mpaaa age if age_decade == 3 & male == 0, mcolor("`color4'") msymbol(O)) || /// 50-60 years, Female
    (scatter mpaaa age if age_decade == 3 & male == 1, mcolor("`color4'") msymbol(S)) || /// 50-60 years, Male
    (scatter mpaaa age if age_decade == 4 & male == 0, mcolor("`color5'") msymbol(O)) || /// 60-70 years, Female
    (scatter mpaaa age if age_decade == 4 & male == 1, mcolor("`color5'") msymbol(S)) || /// 60-70 years, Male
    (scatter mpaaa age if age_decade == 5 & male == 0, mcolor("`color6'") msymbol(O)) || /// 70+ years, Female
    (scatter mpaaa age if age_decade == 5 & male == 1, mcolor("`color6'") msymbol(S)) || /// 70+ years, Male
    (line mpaaa_q50_all age, lcolor(black) lpattern(solid) lwidth(medthick)) || /// 50th quantile
    (line mpaaa_q90_all age, lcolor(gs8) lpattern(solid) lwidth(medthick)) || /// 90th quantile
    (line mpaaa_q95_all age, lcolor(gs12) lpattern(solid) lwidth(medthick)), /// 95th quantile
    legend(order(13 "50th Quantile" 14 "90th Quantile" 15 "95th Quantile")) ///
    xtitle("Age (years)") ytitle("PA:AA Ratio") ///
    title("Quantile Regression of PA:AA Ratio by Age (Both Sexes)") ///
    xlabel(, labsize(medlarge)) ylabel(0.6(0.1)1.3, nogrid labsize(medlarge)) ///
    scheme(white_w3d)
	

	
/* -----
Sex-Stratified Quantile Regression Scatter Plots for PA:AA Ratio
-----*/ 

/* Run Bootstrap Quantile Regressions for the Men */
bsqreg mpaaa c.age if male == 1, quantile(50) reps(500) // 50th Percentile
bsqreg mpaaa c.age if male == 1, quantile(90) reps(500) // 90th Percentile
bsqreg mpaaa c.age if male == 1, quantile(95) reps(500) // 95th Percentile

/* Extract Colors from Spectral Palette */
colorpalette Spectral, n(6) nograph  
local color6 `"`r(p1)'"' // Reverse order
local color5 `"`r(p2)'"'
local color4 `"`r(p3)'"'
local color3 `"`r(p4)'"'
local color2 `"`r(p5)'"'
local color1 `"`r(p6)'"'

/* -----
Generate Quantile Regression Predictions for Males
-----*/ 
qreg mpaaa c.age if male == 1, quantile(50) // Male
predict mpaaa_q50_male if male == 1, xb
qreg mpaaa c.age if male == 1, quantile(90) // Male
predict mpaaa_q90_male if male == 1, xb
qreg mpaaa c.age if male == 1, quantile(95) // Male
predict mpaaa_q95_male if male == 1, xb

/* -----
Male Scatter Plot
-----*/ 
twoway ///
    (scatter mpaaa age if age_decade == 0 & male == 1, mcolor("`color1'") msymbol(S)) || /// <30 years
    (scatter mpaaa age if age_decade == 1 & male == 1, mcolor("`color2'") msymbol(S)) || /// 30-40 years
    (scatter mpaaa age if age_decade == 2 & male == 1, mcolor("`color3'") msymbol(S)) || /// 40-50 years
    (scatter mpaaa age if age_decade == 3 & male == 1, mcolor("`color4'") msymbol(S)) || /// 50-60 years
    (scatter mpaaa age if age_decade == 4 & male == 1, mcolor("`color5'") msymbol(S)) || /// 60-70 years
    (scatter mpaaa age if age_decade == 5 & male == 1, mcolor("`color6'") msymbol(S)) || /// 70+ years
    (line mpaaa_q50_male age, lcolor(black) lpattern(solid) lwidth(medthick)) || /// 50th percentile line
    (line mpaaa_q90_male age, lcolor(gs8) lpattern(solid) lwidth(medthick)) || /// 90th percentile line
    (line mpaaa_q95_male age, lcolor(gs12) lpattern(solid) lwidth(medthick)), /// 95th percentile line
    ///
    legend(order(7 "50th Quantile" 8 "90th Quantile" 9 "95th Quantile") pos(6) row(1)) /// Only quantile labels in legend
    xtitle("Age (years)") ytitle("PA:AA Ratio", margin(medlarge)) ///
    title("Men") ///
    xlabel(, labsize(medlarge)) ylabel(0.6(0.1)1.3, nogrid labsize(medlarge)) /// Adjusted Y-axis range
    scheme(white_w3d) yscale(range(0.6 1.3)) ///
    saving(male_qreg_plot, replace)


/* Run Bootstrap Quantile Regressions for the Women */
bsqreg mpaaa c.age if male == 0, quantile(50) reps(500) // 50th Percentile
bsqreg mpaaa c.age if male == 0, quantile(90) reps(500) // 90th Percentile
bsqreg mpaaa c.age if male == 0, quantile(95) reps(500) // 95th Percentile

/* -----
Generate Quantile Regression Predictions for Females
-----*/ 
qreg mpaaa c.age if male == 0, quantile(50) // Female
predict mpaaa_q50_female if male == 0, xb
qreg mpaaa c.age if male == 0, quantile(90) // Female
predict mpaaa_q90_female if male == 0, xb
qreg mpaaa c.age if male == 0, quantile(95) // Female
predict mpaaa_q95_female if male == 0, xb

/* -----
Female Scatter Plot
-----*/ 
twoway ///
    (scatter mpaaa age if age_decade == 0 & male == 0, mcolor("`color1'") msymbol(O)) || /// <30 years
    (scatter mpaaa age if age_decade == 1 & male == 0, mcolor("`color2'") msymbol(O)) || /// 30-40 years
    (scatter mpaaa age if age_decade == 2 & male == 0, mcolor("`color3'") msymbol(O)) || /// 40-50 years
    (scatter mpaaa age if age_decade == 3 & male == 0, mcolor("`color4'") msymbol(O)) || /// 50-60 years
    (scatter mpaaa age if age_decade == 4 & male == 0, mcolor("`color5'") msymbol(O)) || /// 60-70 years
    (scatter mpaaa age if age_decade == 5 & male == 0, mcolor("`color6'") msymbol(O)) || /// 70+ years
    (line mpaaa_q50_female age, lcolor(black) lpattern(solid) lwidth(medthick)) || /// 50th percentile line
    (line mpaaa_q90_female age, lcolor(gs8) lpattern(solid) lwidth(medthick)) || /// 90th percentile line
    (line mpaaa_q95_female age, lcolor(gs12) lpattern(solid) lwidth(medthick)), /// 95th percentile line
    ///
    legend(order(7 "50th Quantile" 8 "90th Quantile" 9 "95th Quantile") pos(6) row(1)) /// Only quantile labels in legend
    xtitle("Age (years)") ytitle("PA:AA Ratio", margin(medlarge)) ///
    title("Women") ///
    xlabel(, labsize(medlarge)) ylabel(0.6(0.1)1.3, nogrid labsize(medlarge)) /// Adjusted Y-axis range
    scheme(white_w3d) yscale(range(0.6 1.3)) ///
    saving(female_qreg_plot, replace)

/* -----
Combine the Two Plots into a Single Figure
-----*/ 
graph combine male_qreg_plot.gph female_qreg_plot.gph, ///
    title("Sex-Stratified Quantile Regression of PA:AA Ratio by Age") ///
    ycommon xcommon ///
    iscale(*1.2) ///
    rows(1) ///
    graphregion(margin(2 2 2 2)) ///
    xsize(9) ysize(5)
	

/* ------------------

Survival Analysis

--------------------*/

stset lastfollowupyear, failure(death==1)	

/* ------------------

Non Z-score APPROACH 

--------------------*/ 

/* First, split into tertiles */ 

/* MPAD */ 

* Generate tertile cutoffs
xtile mpad_tertile = mpad, n(3)  // Divides into 3 equal-sized groups

* Label the tertile groups
label define mpad_tertile_lbl 1 "MPAD Tertile 1" 2 "MPAD Tertile 2" 3 "MPAD Tertile 3"
label values mpad_tertile mpad_tertile_lbl

* Verify the tertile distribution
tab mpad_tertile

sts test mpad_tertile, logrank 
sts graph, by(mpad_tertile) tmax(10) ci surv ///
 plotopts(lwidth(thick)) ///
 risktable(0(1)10, order(1 "MPA Tertile 1" 2 "MPA Tertile 2" 3 "MPA Tertile 3") title("Number Eligible", size(medium)) size(medsmall)) ///
 xlabel(0(1)10, labsize(medlarge)) ///
 ylabel(,labsize(medlarge)) ///
 xtitle("Year of Followup", size(medlarge)) ///
 ytitle("Survival", size(medlarge)) ///
 legend(order(2 "MPA Tertile 1" 4 "MPA Tertile 2" 6 "MPA Tertile 3") position(6) ring(0) rows(1) size(medsmall)) ///
 title("Survival by Main PA Diameter Tertile", size(large))

stcox i.mpad_tertile
estat concordance
stbrier i.mpad_tertile, bt(12.4819)


/* Add age and sex */ 

stcox i.mpad_tertile c.age i.male
estat concordance
stbrier i.mpad_tertile, bt(12.4819)


stcox i.mpad_tertile c.age male
stcurve, surv at(mpad_tertile=1 male=0.3737 age=52) /// at mean values
          at(mpad_tertile=2 male=0.3737 age=52) ///
          at(mpad_tertile=3 male=0.3737 age=52) ///
          range(0 10) ///
          xlabel(0(1)10, labsize(medlarge)) ///
          ylabel(0(.1)1,labsize(medlarge)) ///
          xtitle("Year of Followup", size(medlarge)) ///
          ytitle("Predicted Survival", size(medlarge)) ///
          title("Predicted Survival by PA Diameter Tertile, Adjusted for Age & Sex", size(large)) ///
          legend(order(1 "MPA Tertile 1" 2 "MPA Tertile 2" 3 "MPA Tertile 3") ///
                 position(6) ring(0) rows(1) size(medsmall)) 

					
					
					
					
					
/* -----
PA:AA Ratio (MPAAA) 
-----*/ 

* Generate tertile cutoffs for MPAAA
xtile mpaaa_tertile = mpaaa, n(3)  // Divides into 3 equal-sized groups

* Label the tertile groups
label define mpaaa_tertile_lbl 1 "MPAAA Tertile 1" 2 "MPAAA Tertile 2" 3 "MPAAA Tertile 3"
label values mpaaa_tertile mpaaa_tertile_lbl

* Verify the tertile distribution
tab mpaaa_tertile

/* -----
Survival Analysis by MPAAA Tertile
-----*/ 

* Log-rank test for survival differences across MPAAA tertiles
sts test mpaaa_tertile, logrank 

* Kaplan-Meier Survival Curve by MPAAA Tertile
sts graph, by(mpaaa_tertile) tmax(10) ci surv ///
 plotopts(lwidth(thick)) ///
 risktable(0(1)10, order(1 "MPAAA Tertile 1" 2 "MPAAA Tertile 2" 3 "MPAAA Tertile 3") title("Number Eligible", size(medium)) size(medsmall)) ///
 xlabel(0(1)10, labsize(medlarge)) ///
 ylabel(,labsize(medlarge)) ///
 xtitle("Year of Followup", size(medlarge)) ///
 ytitle("Survival", size(medlarge)) ///
 legend(order(2 "MPAAA Tertile 1" 4 "MPAAA Tertile 2" 6 "MPAAA Tertile 3") position(6) ring(0) rows(1) size(medsmall)) ///
 title("Survival by PA:AA Ratio Tertile", size(large))

/* -----
Cox Proportional Hazards Model by MPAAA Tertile
-----*/ 
stcox i.mpaaa_tertile
estat concordance

/* -----
Brier Score for MPAAA Tertile Model
-----*/ 
stbrier i.mpaaa_tertile, bt(12.4819)

//adding age
stcox i.mpaaa_tertile c.age i.male
estat concordance

/* -----
Brier Score for MPAAA Tertile Model
-----*/ 
stbrier i.mpaaa_tertile c.age i.male, bt(12.4819)





/* --------- 

Age-specific quantile (tertile) approach

----------*/ 

/* Method = use conditional quantile regression to generate age-sex specific
quantile estimates, then split into tertiles based on that and repeat */ 


/* -------
MPA
--------*/


/* -----
Define Colors for MPAD Predicted Tertiles
-----*/ 

colorpalette Spectral, n(3) nograph  
local color3 `"`r(p1)'"'  // MPAD Tertile 3 (High)
local color2 `"`r(p2)'"'  // MPAD Tertile 2 (Middle)
local color1 `"`r(p3)'"'  // MPAD Tertile 1 (Low)

/* -----
Generate Quantile Regression Predictions for Males
-----*/ 

qreg mpad c.age if male == 1, quantile(0.333) // Male
predict mpad_33_male if male == 1, xb
qreg mpad c.age if male == 1, quantile(0.666) // Male
predict mpad_66_male if male == 1, xb
qreg mpad c.age if male == 1, quantile(95) // Male
predict mpad_95_male if male == 1, xb

/* -----
Male Scatter Plot (Color by MPAD Predicted Tertile)
-----*/ 
twoway ///
    (scatter mpad age if male == 1 & mpad_tertile == 1, mcolor("`color1'") msymbol(S)) || /// MPAD Tertile 1 (Low)
    (scatter mpad age if male == 1 & mpad_tertile == 2, mcolor("`color2'") msymbol(S)) || /// MPAD Tertile 2 (Middle)
    (scatter mpad age if male == 1 & mpad_tertile == 3, mcolor("`color3'") msymbol(S)) || /// MPAD Tertile 3 (High)
    (line mpad_33_male age, lcolor(black) lpattern(solid) lwidth(medthick)) || /// 50th percentile line
    (line mpad_66_male age, lcolor(gs8) lpattern(dash) lwidth(medthick)) || /// 90th percentile line
    (line mpad_95_male age, lcolor(gs12) lpattern(longdash) lwidth(medthick)), /// 95th percentile line
    ///
    legend(order(4 "1st-2nd Tertile" 5 "2nd-3rd Tertile" 6 "95th Percentile" ///
                 1 "MPAD Pred Tertile 1" 2 "MPAD Pred Tertile 2" 3 "MPAD Pred Tertile 3") ///
           pos(6) row(2)) ///
    xtitle("Age (years)") ytitle("Pulmonary Artery Diameter (mm)", margin(medlarge)) ///
    title("Men") ///
    xlabel(, labsize(medlarge)) ylabel(15(5)40, nogrid labsize(medlarge)) /// Adjusted Y-axis range
    scheme(white_w3d) yscale(range(15 40)) ///
    saving(male_qreg_plot, replace)


/* -----
Generate Quantile Regression Predictions for Females
-----*/ 

qreg mpad c.age if male == 0, quantile(0.333) // Female
predict mpad_33_female if male == 0, xb
qreg mpad c.age if male == 0, quantile(0.666) // Female
predict mpad_66_female if male == 0, xb
qreg mpad c.age if male == 0, quantile(95) // Female
predict mpad_95_female if male == 0, xb

/* -----
Female Scatter Plot (Color by MPAD Predicted Tertile)
-----*/ 
twoway ///
    (scatter mpad age if male == 0 & mpad_tertile == 1, mcolor("`color1'") msymbol(O)) || /// MPAD Tertile 1 (Low)
    (scatter mpad age if male == 0 & mpad_tertile == 2, mcolor("`color2'") msymbol(O)) || /// MPAD Tertile 2 (Middle)
    (scatter mpad age if male == 0 & mpad_tertile == 3, mcolor("`color3'") msymbol(O)) || /// MPAD Tertile 3 (High)
    (line mpad_33_female age, lcolor(black) lpattern(solid) lwidth(medthick)) || /// 50th percentile line
    (line mpad_66_female age, lcolor(gs8) lpattern(dash) lwidth(medthick)) || /// 90th percentile line
    (line mpad_95_female age, lcolor(gs12) lpattern(longdash) lwidth(medthick)), /// 95th percentile line
    ///
    legend(order(4 "1st-2nd Tertile" 5 "2nd-3rd Tertile" 6 "95th Quantile" ///
                 1 "MPAD Tertile 1" 2 "MPAD Tertile 2" 3 "MPAD Tertile 3") ///
           pos(6) row(2)) ///
    xtitle("Age (years)") ytitle("Pulmonary Artery Diameter (mm)", margin(medlarge)) ///
    title("Women") ///
    xlabel(, labsize(medlarge)) ylabel(15(5)40, nogrid labsize(medlarge)) /// Adjusted Y-axis range
    scheme(white_w3d) yscale(range(15 40)) ///
    saving(female_qreg_plot, replace)


/* -----
Combine the Two Plots into a Single Figure
-----*/ 
graph combine male_qreg_plot.gph female_qreg_plot.gph, ///
    title("Sex-Stratified Quantile Regression of MPAD by Predicted Tertile") ///
    ycommon xcommon ///
    iscale(*1.2) ///
    rows(1) ///
    graphregion(margin(2 2 2 2)) ///
    xsize(9) ysize(5)




/* -----
Step 1: Compute 33.3rd and 66.6th Percentiles Using Quantile Regression
-----*/
* Run Quantile Regression Separately for Males
qreg mpad c.age if male == 1, quantile(0.333)
predict mpad_q33_male if male == 1, xb

qreg mpad c.age if male == 1, quantile(0.666)
predict mpad_q66_male if male == 1, xb

* Run Quantile Regression Separately for Females
qreg mpad c.age if male == 0, quantile(0.333)
predict mpad_q33_female if male == 0, xb

qreg mpad c.age if male == 0, quantile(0.666)
predict mpad_q66_female if male == 0, xb

/* -----
Step 2: Categorize Individuals Into Predicted Tertiles
-----*/

* Create an empty variable for tertile classification
gen mpad_pred_tertile = .

* Assign individuals into tertiles based on their predicted quantile cutoffs
replace mpad_pred_tertile = 1 if male == 1 & mpad < mpad_q33_male
replace mpad_pred_tertile = 2 if male == 1 & mpad >= mpad_q33_male & mpad < mpad_q66_male
replace mpad_pred_tertile = 3 if male == 1 & mpad >= mpad_q66_male

replace mpad_pred_tertile = 1 if male == 0 & mpad < mpad_q33_female
replace mpad_pred_tertile = 2 if male == 0 & mpad >= mpad_q33_female & mpad < mpad_q66_female
replace mpad_pred_tertile = 3 if male == 0 & mpad >= mpad_q66_female

/* -----
Step 3: Label the Categories
-----*/

label define mpad_pred_tertile_lbl 1 "MPA Pred Tertile 1" 2 "MPAD Pred Tertile 2" 3 "MPAD Pred Tertile 3"
label values mpad_pred_tertile mpad_pred_tertile_lbl

* Verify the tertile distribution
tab mpad_pred_tertile mpad_tertile
kappaetc mpad_pred_tertile mpad_tertile


tab3way mpad_pred_tertile mpad_tertile age_decade 

/* Repeat analyses with predicted tertiles */ 

* Log-rank test for survival differences across MPAAA tertiles
sts test mpad_pred_tertile, logrank 

* Kaplan-Meier Survival Curve by MPAAA Tertile
sts graph, by(mpad_pred_tertile) tmax(10) ci surv ///
 plotopts(lwidth(thick)) ///
 risktable(0(1)10, order(1 "Age-Sex MPA Tertile 1" 2 "Age-Sex MPA Tertile 2" 3 "Age-Sex MPA Tertile 3") title("Number Eligible", size(medium)) size(medsmall)) ///
 xlabel(0(1)10, labsize(medlarge)) ///
 ylabel(,labsize(medlarge)) ///
 xtitle("Year of Followup", size(medlarge)) ///
 ytitle("Survival", size(medlarge)) ///
 legend(order(2 "Age-Sex MPA Tertile 1" 4 "Age-Sex MPA Tertile 2" 6 "Age-Sex MPA Tertile 3") position(6) ring(0) rows(3) size(medsmall)) ///
 title("Survival by Age-Sex MPA Tertile", size(large))

/* -----
Cox Proportional Hazards Model by MPAD Tertile
-----*/ 
stcox i.mpad_pred_tertile
estat concordance

/* -----
Brier Score for MPAD Tertile Model
-----*/ 
stbrier i.mpad_tertile, bt(12.4819)


/* -----
Cox Proportional Hazards Model by MPAD Tertile
-----*/ 
stcox i.mpad_pred_tertile c.age i.male
estat concordance

/* -----
Brier Score for MPAD Tertile Model
-----*/ 
stbrier i.mpad_tertile, bt(12.4819)






/* ---------
MPA:AA 
----------*/ 

/* -----
Define Colors for MPAAA Predicted Tertiles
-----*/ 

colorpalette Spectral, n(3) nograph  
local color3 `"`r(p1)'"'  // MPAAA Tertile 3 (High)
local color2 `"`r(p2)'"'  // MPAAA Tertile 2 (Middle)
local color1 `"`r(p3)'"'  // MPAAA Tertile 1 (Low)

/* -----
Generate Quantile Regression Predictions for Males
-----*/ 

qreg mpaaa c.age if male == 1, quantile(0.333) // Male
predict mpaaa_33_male if male == 1, xb
qreg mpaaa c.age if male == 1, quantile(0.666) // Male
predict mpaaa_66_male if male == 1, xb
qreg mpaaa c.age if male == 1, quantile(0.95) // Male
predict mpaaa_95_male if male == 1, xb

/* -----
Male Scatter Plot (Color by MPAAA Predicted Tertile)
-----*/ 
twoway ///
    (scatter mpaaa age if male == 1 & mpaaa_tertile == 1, mcolor("`color1'") msymbol(S)) || /// MPAAA Tertile 1 (Low)
    (scatter mpaaa age if male == 1 & mpaaa_tertile == 2, mcolor("`color2'") msymbol(S)) || /// MPAAA Tertile 2 (Middle)
    (scatter mpaaa age if male == 1 & mpaaa_tertile == 3, mcolor("`color3'") msymbol(S)) || /// MPAAA Tertile 3 (High)
    (line mpaaa_33_male age, lcolor(black) lpattern(solid) lwidth(medthick)) || /// 33rd percentile line
    (line mpaaa_66_male age, lcolor(gs8) lpattern(dash) lwidth(medthick)) || /// 66th percentile line
    (line mpaaa_95_male age, lcolor(gs12) lpattern(longdash) lwidth(medthick)), /// 95th percentile line
    ///
    legend(order(4 "1st-2nd Age-Sex Tertile" 5 "2nd-3rd Age-Sex Tertile" 6 "95th Age-Sex Percentile" ///
                 1 "MPAAA Tertile 1" 2 "MPAAA Tertile 2" 3 "MPAAA Tertile 3") ///
           pos(6) row(2)) ///
    xtitle("Age (years)") ytitle("PA:AA Ratio", margin(medlarge)) ///
    title("Men") ///
    xlabel(, labsize(medlarge)) ylabel(0.6(0.1)1.3, nogrid labsize(medlarge)) /// Adjusted Y-axis range
    scheme(white_w3d) yscale(range(0.6 1.3)) ///
    saving(male_qreg_plot, replace)


/* -----
Generate Quantile Regression Predictions for Females
-----*/ 

qreg mpaaa c.age if male == 0, quantile(0.333) // Female
predict mpaaa_33_female if male == 0, xb
qreg mpaaa c.age if male == 0, quantile(0.666) // Female
predict mpaaa_66_female if male == 0, xb
qreg mpaaa c.age if male == 0, quantile(0.95) // Female
predict mpaaa_95_female if male == 0, xb

/* -----
Female Scatter Plot (Color by MPAAA Predicted Tertile)
-----*/ 
twoway ///
    (scatter mpaaa age if male == 0 & mpaaa_tertile == 1, mcolor("`color1'") msymbol(O)) || /// MPAAA Tertile 1 (Low)
    (scatter mpaaa age if male == 0 & mpaaa_tertile == 2, mcolor("`color2'") msymbol(O)) || /// MPAAA Tertile 2 (Middle)
    (scatter mpaaa age if male == 0 & mpaaa_tertile == 3, mcolor("`color3'") msymbol(O)) || /// MPAAA Tertile 3 (High)
    (line mpaaa_33_female age, lcolor(black) lpattern(solid) lwidth(medthick)) || /// 33rd percentile line
    (line mpaaa_66_female age, lcolor(gs8) lpattern(dash) lwidth(medthick)) || /// 66th percentile line
    (line mpaaa_95_female age, lcolor(gs12) lpattern(longdash) lwidth(medthick)), /// 95th percentile line
    ///
    legend(order(4 "1st-2nd Age-Sex Tertile" 5 "2nd-3rd Age-Sex Tertile" 6 "95th Age-sex Quantile" ///
                 1 "MPAAA Tertile 1" 2 "MPAAA Tertile 2" 3 "MPAAA Tertile 3") ///
           pos(6) row(2)) ///
    xtitle("Age (years)") ytitle("PA:AA Ratio", margin(medlarge)) ///
    title("Women") ///
    xlabel(, labsize(medlarge)) ylabel(0.6(0.1)1.3, nogrid labsize(medlarge)) /// Adjusted Y-axis range
    scheme(white_w3d) yscale(range(0.6 1.3)) ///
    saving(female_qreg_plot, replace)


/* -----
Combine the Two Plots into a Single Figure
-----*/ 
graph combine male_qreg_plot.gph female_qreg_plot.gph, ///
    title("Sex-Stratified Quantile Regression of PA:AA Ratio by Predicted Tertile") ///
    ycommon xcommon ///
    iscale(*1.2) ///
    rows(1) ///
    graphregion(margin(2 2 2 2)) ///
    xsize(9) ysize(5)




/* -----
Step 1: Compute 33.3rd and 66.6th Percentiles Using Quantile Regression for MPAAA
-----*/

* Run Quantile Regression Separately for Males
qreg mpaaa c.age if male == 1, quantile(0.333)
predict mpaaa_q33_male if male == 1, xb

qreg mpaaa c.age if male == 1, quantile(0.666)
predict mpaaa_q66_male if male == 1, xb

* Run Quantile Regression Separately for Females
qreg mpaaa c.age if male == 0, quantile(0.333)
predict mpaaa_q33_female if male == 0, xb

qreg mpaaa c.age if male == 0, quantile(0.666)
predict mpaaa_q66_female if male == 0, xb

/* -----
Step 2: Categorize Individuals Into Predicted Tertiles for MPAAA
-----*/

* Create an empty variable for tertile classification
gen mpaaa_pred_tertile = .

* Assign individuals into tertiles based on their predicted quantile cutoffs
replace mpaaa_pred_tertile = 1 if male == 1 & mpaaa < mpaaa_q33_male
replace mpaaa_pred_tertile = 2 if male == 1 & mpaaa >= mpaaa_q33_male & mpaaa < mpaaa_q66_male
replace mpaaa_pred_tertile = 3 if male == 1 & mpaaa >= mpaaa_q66_male

replace mpaaa_pred_tertile = 1 if male == 0 & mpaaa < mpaaa_q33_female
replace mpaaa_pred_tertile = 2 if male == 0 & mpaaa >= mpaaa_q33_female & mpaaa < mpaaa_q66_female
replace mpaaa_pred_tertile = 3 if male == 0 & mpaaa >= mpaaa_q66_female

/* -----
Step 3: Label the Categories for MPAAA Predicted Tertiles
-----*/

label define mpaaa_pred_tertile_lbl 1 "PA:AA Pred Tertile 1" 2 "PA:AA Pred Tertile 2" 3 "PA:AA Pred Tertile 3"
label values mpaaa_pred_tertile mpaaa_pred_tertile_lbl

* Verify the tertile distribution
tab mpaaa_pred_tertile mpaaa_tertile
kappaetc mpaaa_pred_tertile mpaaa_tertile

tab3way mpaaa_pred_tertile mpaaa_tertile age_decade
	
	

	
	
/* Survival analysis */ 

/* -----
Log-Rank Test for Survival Differences Across PA:AA Predicted Tertiles
-----*/ 
sts test mpaaa_pred_tertile, logrank 

/* -----
Kaplan-Meier Survival Curve by PA:AA Predicted Tertile
-----*/ 
sts graph, by(mpaaa_pred_tertile) tmax(10) ci surv ///
 plotopts(lwidth(thick)) ///
 risktable(0(1)10, order(1 "Age-Sex PA:AA Tertile 1" 2 "Age-Sex PA:AA Tertile 2" 3 "Age-Sex PA:AA Tertile 3") title("Number Eligible", size(medium)) size(medsmall)) ///
 xlabel(0(1)10, labsize(medlarge)) ///
 ylabel(,labsize(medlarge)) ///
 xtitle("Year of Followup", size(medlarge)) ///
 ytitle("Survival", size(medlarge)) ///
 legend(order(2 "Age-Sex PA:AA Tertile 1" 4 "Age-Sex PA:AA Tertile 2" 6 "Age-Sex PA:AA Tertile 3") position(6) ring(0) rows(3) size(medsmall)) ///
 title("Survival by Age-Sex PA:AA Tertile", size(large))
	
	
	
	
/* -----
Cox Proportional Hazards Model by PA:AA Predicted Tertile
-----*/ 
	stcox i.mpaaa_pred_tertile
	estat concordance

/* -----
Brier Score for PA:AA Tertile Model
-----*/ 
stbrier i.mpaaa_tertile, bt(12.4819)


/* -----
Cox Proportional Hazards Model by PA:AA Predicted Tertile with Age and Sex Adjustment
-----*/ 
stcox i.mpaaa_pred_tertile c.age i.male
estat concordance

/* -----
Brier Score for PA:AA Tertile Model
-----*/ 
stbrier i.mpaaa_tertile, bt(12.4819)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	//Old analysis 






/*PA Size*/ 

/* Without Age-sex adjustment */ 

preserve 
mkspline2 rc = sex_norm_mpad, cubic nknots(4) displayknots            
assert float(sex_norm_mpad) == float(rc1)
stbrier c.age rc*, bt(12.4819)
stcox c.age rc*
estat ic
estat concordance
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
restore

preserve 
* Generate splines for sex-normalized PA diameter
mkspline2 rc = sex_norm_mpad, cubic nknots(4) displayknots            
assert float(sex_norm_mpad) == float(rc1)

* Fit the Cox proportional hazards model
stcox c.age rc*

* Extract levels of sex_norm_mpad within range (19-39)
levelsof sex_norm_mpad if inrange(sex_norm_mpad, 19, 39), local(levels)

* Generate predicted hazard ratios relative to reference value (25.1)
xblc rc*, covname(sex_norm_mpad) at(`r(levels)') reference(25.1) eform generate(sex_norm_mpad hazard lb ub)

* Store predicted hazard ratio as new variable
gen mpad_hazard_pr = hazard




/* AA Size */ 
/* Without Age-sex adjustment */ 

preserve 
mkspline2 rc = ascendingaorta, cubic nknots(4) displayknots            
assert float(ascendingaorta) == float(rc1)
stbrier c.age rc*, bt(12.4819)
stcox c.age rc*
estat ic
estat concordance
levelsof ascendingaorta if inrange(ascendingaorta, 18, 50), local(levels)
xblc rc*, covname(ascendingaorta) at(`r(levels)') reference(30) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) ///
 (line or pa, sort lc(black) lp(l)) if inrange(pa,18,50), ///
 yscale(log extend) ///
 scheme(cleanplots) ///
 legend(off) ///
 xlabel(18(3)50, labsize(large)) ///
 xmtick(18(3)50) ///
 ylabel(0.125 0.25 0.5 1 2 4, angle(horiz) format(%2.1fc) labsize(large)) ///
 ytitle("Hazard Ratio of Mortality", size(large)) ///
 xtitle("Asc Aorta Diameter (mm)", size(large)) ///
 title("Mortality Risk by AA Diameter", size(vlarge)) ///
 yline(1, lp("shortdash") lc(gs10)) 
graph export "Results and Figures/$S_DATE/HR AA Splines - Whole Cohort Adj AgeSex.png", as(png) name("Graph") replace
restore

/* MPAAA */ 
/* Without Age-sex adjustment */ 
stcox mpaaa 
estat ic
estat concordance
stbrier mpaaa, bt(12.4819)


preserve 
gen mpaaa_rounded = round(mpaaa, 0.01)
mkspline2 rc = mpaaa_rounded, cubic nknots(4) displayknots            
assert float(mpaaa_rounded) == float(rc1)
stbrier c.age rc*, bt(12.4819)
stcox rc*
estat ic
estat concordance
levelsof mpaaa_rounded if inrange(mpaaa_rounded, 0.47, 1.4), local(levels)
xblc rc*, covname(mpaaa_rounded) at(`r(levels)') reference(0.77) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) ///
 (line or pa, sort lc(black) lp(l)) if inrange(pa,0.47,1.4), ///
 yscale(log extend) ///
 scheme(cleanplots) ///
 legend(off) ///
 xlabel(0.5(.1)1.4, labsize(large)) ///
 xmtick(0.5(0.05)1.4) ///
 ylabel(1 2 4, angle(horiz) format(%2.1fc) labsize(large)) ///
 ytitle("Hazard Ratio of Mortality", size(large)) ///
 xtitle("PA:AA", size(large)) ///
 title("PA:AA Mortality Risk", size(vlarge)) ///
 yline(1, lp("shortdash") lc(gs10)) ///
 xline(28, lp("shortdash_dot") lc(gs10))
graph export "Results and Figures/$S_DATE/HR PAAA Splines - Whole Cohort.png", as(png) name("Graph") replace
restore







/* ------------------
Z-SCORE APPROACH

Note: this is similar to the "LMS Method" used for things like PFTs - but without the skew part since the distributions don't seem very skewed. 

Overall aim: 
-- model age, sex specific PAd z-score - modeling both the changing mean and the changing spread. 
-- then, create cox-HR directly from that (may or may not be non-linear)

Downside? Heteroscedasticity

the z‐scores represent "how far an individual's PA measurement is from the expected mean for someone of the same age and sex in this sample," - it's a valid internal standardization, but not necessarily representative of population norms

------------------ */ 


/* PA diameter */ 

/* Generate within-sample z-score */ 
regress mpad c.age i.male
predict mpad_hat
egen mpad_hat_z = std(mpad_hat) 
label variable mpad_hat_z "Age, Sex-adjusted PAd Z-score"

stcox mpad_hat_z
estat ic
estat concordance
stbrier mpad_hat_z, bt(12.4819) //not sure why brier is higher here?  - need to look into this. c-statistic much better.

hist mpad_hat_z

/* Notably, it's not much different than (log)linear */ 
preserve 
gen mpad_hat_z_rounded = round(mpad_hat_z, 0.05)
mkspline2 rc = mpad_hat_z_rounded, cubic nknots(4) displayknots            
assert float(mpad_hat_z_rounded) == float(rc1)
stbrier rc*, bt(12.4819)
stcox rc* 
estat ic
estat concordance 
levelsof mpad_hat_z_rounded if inrange(mpad_hat_z_rounded, -2, 2.65), local(levels)
xblc rc*, covname(mpad_hat_z_rounded) at(`r(levels)') reference(0) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) ///
 (line or pa, sort lc(black) lp(l)) if inrange(pa,-2,2.65), ///
 yscale(log extend) ///
 scheme(cleanplots) ///
 legend(off) ///
 xlabel(-2(0.5)2.65, labsize(large)) ///
 xmtick(-2(0.25)2.65) ///
 ylabel(0.125 0.25 0.5 1 2 4, angle(horiz) format(%2.1fc) labsize(large)) ///
 ytitle("Hazard Ratio of Mortality", size(large)) ///
 xtitle("Sex-Normalized Pulmonary Artery Diameter (mm)", size(large)) ///
 title("Age-adjusted Mortality Risk by PA Diameter", size(vlarge)) ///
 yline(1, lp("shortdash") lc(gs10)) ///
 xline(28, lp("shortdash_dot") lc(gs10))
graph export "Results and Figures/$S_DATE/HR PAd Z-score Splines - Whole Cohort Adj AgeSex.png", as(png) name("Graph") replace
restore



/* Ascending Aorta */ 

/* Generate within-sample z-score */ 
regress ascendingaorta c.age i.male
predict aa_hat
egen aa_hat_z = std(aa_hat) 
label variable aa_hat_z "Age, Sex-adjusted Asc Aorta Z-score"


/* With Age & Sex Adjustment */ 
stcox aa_hat_z
estat ic
estat concordance
stbrier aa_hat_z, bt(12.4819) //not sure why brier is higher here?  - need to look into this. c-statistic much better.

hist aa_hat_z

/* Very (log) linear */ 
preserve 
gen aa_hat_z_rounded = round(aa_hat_z, 0.05)
mkspline2 rc = aa_hat_z_rounded, cubic nknots(4) displayknots            
assert float(aa_hat_z_rounded) == float(rc1)
stbrier rc*, bt(12.4819)
stcox rc*
estat ic
estat concordance
levelsof aa_hat_z_rounded if inrange(aa_hat_z_rounded, -2, 2.5), local(levels)
xblc rc*, covname(aa_hat_z_rounded) at(`r(levels)') reference(-0.1) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) ///
 (line or pa, sort lc(black) lp(l)) if inrange(pa,-2,2.5), ///
 yscale(log extend) ///
 scheme(cleanplots) ///
 legend(off) ///
 xlabel(-2(0.5)2.5, labsize(large)) ///
 xmtick(-2(0.25)2.5) ///
 ylabel(0.125 0.25 0.5 1 2 4, angle(horiz) format(%2.1fc) labsize(large)) ///
 ytitle("Hazard Ratio of Mortality", size(large)) ///
 xtitle("Asc Aorta Diameter z-score (mm)", size(large)) ///
 title("Mortality Risk by AA Diameter Z-score", size(vlarge)) ///
 yline(1, lp("shortdash") lc(gs10)) 
graph export "Results and Figures/$S_DATE/HR AA z-score Splines - Whole Cohort Adj AgeSex.png", as(png) name("Graph") replace
restore



/* Ratio */ 


/* Generate within-sample z-score */ 
regress mpaaa c.age i.male
predict mpaaa_hat
egen mpaaa_hat_z = std(mpaaa_hat)
label variable mpaaa_hat_z "Age, Sex-adjusted PA:AA Z-score"


/* With Age-sex adjustment */ 
stcox mpaaa_hat_z
estat ic
estat concordance
stbrier mpaaa_hat_z, bt(12.4819)

preserve 
gen mpaaa_hat_z_rounded = round(mpaaa_hat_z, 0.02)
mkspline2 rc = mpaaa_hat_z_rounded, cubic nknots(4) displayknots            
assert float(mpaaa_hat_z_rounded) == float(rc1)
stbrier rc*, bt(12.4819)
stcox rc*
estat ic
estat concordance
levelsof mpaaa_hat_z_rounded if inrange(mpaaa_hat_z_rounded, -2.5, 2), local(levels)
xblc rc*, covname(mpaaa_hat_z_rounded) at(`r(levels)') reference(-1) eform generate(pa or lb ub)
twoway (line lb ub pa, sort lc(black black) lp(longdash longdash)) ///
 (line or pa, sort lc(black) lp(l)) if inrange(pa,-2.51,2.01), ///
 yscale(log extend) ///
 scheme(cleanplots) ///
 legend(off) ///
 xlabel(-2.5(.5)2, labsize(large)) ///
 xmtick(-2.5(0.25)2) ///
 ylabel(1 2 4, angle(horiz) format(%2.1fc) labsize(large)) ///
 ytitle("Hazard Ratio of Mortality", size(large)) ///
 xtitle("PA:AA z-score", size(large)) ///
 title("PA:AA z-score Mortality Risk", size(vlarge)) ///
 yline(1, lp("shortdash") lc(gs10)) ///
 xline(28, lp("shortdash_dot") lc(gs10))
graph export "Results and Figures/$S_DATE/HR PAAA z score Splines - Whole Cohort Adj AgeSex.png", as(png) name("Graph") replace
restore


//TODO: bivariate splines predicted heatmap like in the HCO3 paper? z-scored and not. 
//splines

stcox c.mpad_hat_z c.aa_hat_z
estat ic
estat concordance
stbrier mpaaa_hat_z, bt(12.4819)


stcox c.mpad_hat_z##c.aa_hat_z
estat ic
estat concordance
stbrier mpaaa_hat_z, bt(12.4819)





/* To do a version that does variable std deviation: 
* Regress PA on age (linearly) and sex
regress PA c.age i.sex

* Obtain residuals
predict double e, resid
gen double abs_e = abs(e)

* Model log(std dev)
gen double ln_abs_e = ln(abs_e)
regress ln_abs_e c.age i.sex

* Predict the fitted (linear) mean from Step 1
predict double meanhat, xb

* Predict log-SD from Step 2
predict double ln_sdhat, xb
gen double sdhat = exp(ln_sdhat)

* Compute the heteroskedasticity-adjusted z-score
gen double z_score = (PA - meanhat) / sdhat
*/ 




/* -----------------
QUANTILE REGRESSION Approach

Overall aim: 
-- model age, sex specific PAd quantile
-- then, create splined relationship between  (need this to capture presumed non-linear relationship)

Why? you can adjust for the mortality rate on the back-end, or swap it to the front end
------------------ */


sqreg mpad male age, q(0.1 0.5 0.9)
estat coefplot

predict q10, xb q(0.1)
predict q50, xb q(0.5)
predict q90, xb q(0.9)

twoway ///
    scatter mpad age, msymbol(O) mcolor(black) || /// Scatter plot
    line q10 age, lcolor(red) lpattern(dash) || /// 10th quantile line
    line q50 age, lcolor(blue) lpattern(solid) || /// 50th quantile line (median)
    line q90 age, lcolor(green) lpattern(dash), /// 90th quantile line
    xlabel(, labsize(medlarge)) ylabel(, labsize(medlarge)) ///
    xtitle("Age (years)") ytitle("PA Size (mm)") ///
    title("Quantile Regression: PA Diameter vs Age") ///
    legend(order(2 "10th Quantile" 3 "50th Quantile (Median)" 4 "90th Quantile"))

	
	
/* Examples for plots


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

*/ 
	

	
sqreg mpaaa c.age, q(0.1 0.5 0.9)
	
	
regress mpad c.age i.male
sqreg mpad c.age i.male, q(0.1 0.5 0.9)


qreg mpad age, quantile(0.1)
predict q10, xb

qreg mpad age, quantile(0.5)
predict q50, xb

qreg mpad age, quantile(0.9)
predict q90, xb

twoway (scatter mpad age) ///
       (line q10 age, lcolor(red) lpattern(dash)) ///
       (line q50 age, lcolor(blue) lpattern(solid)) ///
       (line q90 age, lcolor(green) lpattern(dash)), ///
       xlabel(, labsize(medlarge)) ylabel(, labsize(medlarge)) ///
       xtitle("Age (years)") ytitle("PA Size (mm)") ///
       title("Quantile Regression: PA Diameter vs Age") ///
       legend(order(2 "10th Quantile" 3 "50th Quantile (Median)" 4 "90th Quantile"))



/*
	Title:		Quantile Regression STATA Code
	Author:		Aayush Khadka, Jilly Hebert, and Anusha Vable
	Institution: University of California, San Francisco
*/
	
	****************************************************************************
	clear all
	
	** Setting directories and loading data
	
	
	****************************************************************************
	
	** Mean Model
	
	****************************************************************************
	
	** Create OLS results sheet
	putexcel set `saveloc'qr_results, modify sheet("OLS")
	putexcel A1 = "quantile"	
	putexcel B1 = "coef" //Saves school year estimate
	putexcel C1 = "lci"	
	putexcel D1 = "uci"	
	
	
	** Bootstrap 95% CI's
	regress sbp c.schlyrs c.age c.age2 i.gender i.race c.rameduc c.rafeduc //
	i.southern i.year, vce(boostrap, reps(500))
	
	** Other modeling option
	*regress sbp c.schlyrs c.age c.age2 i.female i.black i.latinx //
	*c.rameduc c.rafeduc i.southern i.y08 i.y10 i.y12 i.y14 i.y16 i.y18, //
	*vce(bootstrap, reps(500))
	
	
	** Extracting results
	matrix param = r(table)	
	
	** Storing results
	putexcel A2 = -0.10
	putexcel B2 = param[1,1]
	putexcel C2 = param[5,1]
	putexcel D2 = param[6,1]
	
	
	****************************************************************************
	
	** Conditional Quantile Regression (CQR)
	
	****************************************************************************
	
	** Create CQR results sheet
	putexcel set `saveloc'qr_results.xlsx, modify sheet("CQR")
	putexcel A1 = "quantile"			
	putexcel B1 = "coef" //Saves school year estimate
	putexcel C1 = "lci"	
	putexcel D1 = "uci"	
	
	** Creating a counter
	local c = 2
	
	**Bootstrap 95% CIs
	forval i=0.1(0.01)0.9 {
		
		bsqreg sbp c.schlyrs c.age c.age2 i.gender i.race //
		c.rameduc c.rafeduc i.southern i.year, quantile(`i') reps(500)
		
		** Other modeling option
		*bsqreg sbp c.schlyrs c.age c.age2 i.female i.black i.latinx //
		*c.rameduc c.rafeduc i.southern i.y08 i.y10 i.y12 i.y14 i.y16 i.y18, //
		*quantile(`i') reps(500)
		
		** Extracting results
		matrix param = r(table)
		
		** Storing results
		putexcel A`c' = `i'*100
		putexcel B`c' = param[1,1]
		putexcel C`c' = param[5,1]
		putexcel D`c' = param[6,1]
		
		** Updating counter
		local c = `c' + 1		
		
	}
	
	
	****************************************************************************
	
	** Unconditional Quantile Regression (UQR)
	
	****************************************************************************
	
	** Create CQR results sheet
	putexcel set `saveloc'qr_results.xlsx, modify sheet("UQR")
	putexcel A1 = "quantile"			
	putexcel B1 = "coef" //Saves school year estimate
	putexcel C1 = "lci"	
	putexcel D1 = "uci"	
	
	** Creating a counter
	local c = 2
	
	**Bootstrap 95% CI's
	forval i=10(1)90 {
		
		rifhdreg sbp c.schlyrs c.age c.age2 i.gender i.race c.rameduc //
		c.rafeduc i.southern i.year, rif(q(`i')) vce(bootstrap, reps(500))
		
		** Other modeling option
		*rifhdreg sbp c.schlyrs c.age c.age2 i.female i.black i.latinx //
		*c.rameduc c.rafeduc i.southern i.y08 i.y10 i.y12 i.y14 i.y16 i.y18, //
		*rif(q(`i')) vce(bootstrap, reps(500))
		
		** Extracting results
		matrix param = r(table)
		
		** Storing results
		putexcel A`c' = `i'
		putexcel B`c' = param[1,1]
		putexcel C`c' = param[5,1]
		putexcel D`c' = param[6,1]
		
		** Updating counter
		local c = `c' + 1		
		
	}
	
	
	
	   