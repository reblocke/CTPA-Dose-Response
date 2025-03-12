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


//Categorizations of continuous variables

//[ ] remove? 
recode age min/30=0 30/40=1 40/50=2 50/60=3 60/70=4 70/max=5, gen(age_decade)
label define age_dec_lab 0 "<30 years" 1 "30-40 years" 2 "40-50 years" 3 "50-60 years" 4 "60-70 years" 5 "70+ years"
label variable age_decade "Age (by decade)"
label values age_decade age_dec_lab 


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
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol total(before) ///
saving("Results and Figures/$S_DATE/Table 1 PA enlargement by Age.xlsx", replace)


//ridgelines plots to visualize overall change by age. (probably want to account for sex somehow)
ridgeline mpad, by(age_decade) yline ylw(0.2) overlap(1.7) ylc(blue) ylp(dot) ////
	labpos(right) bwid(1.2) laboffset(-3) showstats xlabel(15(5)40) ///
	palette(CET C6) alpha(75) xtitle("Pulmonary Artery Diameter (mm)")

ridgeline ascendingaorta, by(age_decade) yline ylw(0.2) overlap(1.7) ylc(blue) ylp(dot) ////
	labpos(right) bwid(1.2) laboffset(-3) showstats xlabel(20(5)45) ///
	palette(CET C6) alpha(75) xtitle("Ascending Aorta Diameter (mm)")
	
ridgeline mpaaa, by(age_decade) yline ylw(0.2) overlap(1.7) ylc(blue) ylp(dot) ////
	labpos(right) bwid(0.05) laboffset(-0.05) showstats xlabel(0.6(.1)1.3) ///
	palette(CET C6) alpha(75) xtitle("PA:AA Ratio")
	
// [ ] Maybe make a scatter plot with age as color and MPA on one axes, aa on the other? 


twoway ///
    (scatter ascendingaorta mpad if age_decade == 0, mcolor(blue) msymbol(O) ) || /// <30 years
    (scatter ascendingaorta mpad if age_decade == 1, mcolor(red) msymbol(O) ) || /// 30-40 years
    (scatter ascendingaorta mpad if age_decade == 2, mcolor(green) msymbol(O) ) || /// 40-50 years
    (scatter ascendingaorta mpad if age_decade == 3, mcolor(orange) msymbol(O) ) || /// 50-60 years
    (scatter ascendingaorta mpad if age_decade == 4, mcolor(purple) msymbol(O) ) || /// 60-70 years
    (scatter ascendingaorta mpad if age_decade == 5, mcolor(black) msymbol(O) ), /// 70+ years
    ///
    legend(order(1 "<30 years" 2 "30-40 years" 3 "40-50 years" 4 "50-60 years" 5 "60-70 years" 6 "70+ years")) ///
    xtitle("Pulmonary Artery Diameter (mm)") ytitle("Ascending Aorta Diameter (mm)") ///
    title("Scatterplot of PA vs Ascending Aorta by Age Group") ///
    xlabel(, labsize(medlarge)) ylabel(, labsize(medlarge)) ///
    scheme(white_w3d)

	//todo: consider adding z-score or quantile lines to this? 


/* ------------------
Z-SCORE APPROACH

Overall aim: 
-- model age, sex specific PAd z-score - modeling both the changing mean and the changing spread. 
-- then, create cox-HR directly from that (may or may not be non-linear)

Downside? Heteroscedasticity

------------------ */ 


/*
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
	
	
	
	   