clear all
set more off
set trace off

cd /home/vegayon/Dropbox/usc/clases/2019_PM511B/assignments/04

use vitals, clear
global outregopts tex(fragment) label dec(3)

// Previous variables created
gen hibp = sbp > 140 if sbp != .
gen bmi = 703 * weight / height^2 if (weight != . & sbp != .)
gen overweight = bmi >= 25 if bmi != .

lab var hibp "High Blood Pressure (yes/no)"
lab var overw "Overweight (yes/no)"
lab var age "Age (years)"

egen     age_mean = mean(age)
gen age_cent = age-age_mean
lab var age_cent "Age (mean centered)"

gen age_cent60 = age - 60
lab var age_cent60 "Age centered at 60"


// Program to predict disease for individuals with 0 or 5 in s
prog def predictedprobs, rclass
	preserve
	quietly {
		drop _all
		set obs 2
		gen age_cent60 = 70*(_n == 1) + 50*(_n == 2) - 60
		gen overweight = _n == 1
		predict pred_prob
		noi list

	}
	local p0 = round(pred_prob[1], .01)
	local p5 = round(pred_prob[2], .01)
	local rr = round(pred_prob[2]/pred_prob[1], .01)
	di "P(D|S=0): `p0', P(D|S=5): `p5', RR: `rr'"
	restore
	
	// This , together with the 'rclass' option in the program definition,
	// allows the program returning stuff, in this case, locals. We use this
	// to add them to the final table by calling "r(p0)", "r(p5)", and "r(rr)".
	return add
	return local p0 : di %03.2f `p0'
	return local p5 : di %03.2f `p5'
	return local rr : di %03.2f `rr'
end

/*******************************************************************************
PART 1
*******************************************************************************/

* a
glm hibp age_cent60, link(logit) family(binomial)
outreg2 using part1.tex, replace stats(coef se ci tstat) $outregopts

glm hibp overweight, link(logit) family(binomial)
outreg2 using part1.tex, append stats(coef se ci tstat) $outregopts

* b
glm hibp age_cent60 overweight, link(logit) family(binomial)
outreg2 using part1.tex, append stats(coef se ci tstat) $outregopts

predictedprobs

/*******************************************************************************
PART 2
*******************************************************************************/

* a
sort age
egen age_quant = cut(age), group(4) lab
lab var age_quant "Age quantile"

logit hibp i.age_quant
mat def betas = e(b)
bysort age: egen age_midpoint = mean(age_quant)
lab var age_mid "Midpoint of the Age quantile"

gen beta = .
forval i = 1/4 {
	replace beta = betas[1, `i'] if age_quant == (`i' - 1)
}
// replace beta = 0 if beta == .
lab var beta "Logit coefficient"

twoway scatter beta age_midpoint || line beta age_midpoint
graph export part2_hibp_vs_age_quant.eps, replace

* b
lowess hibp age, addplot(scatter hibp age) legend(off) logit
graph export part2_hibp_vs_age_lowess.eps, replace

* c
fp <age>, replace: logit hibp <age>
fp <age>, dim(1) replace: logit hibp <age>

/*******************************************************************************
PART 3
*******************************************************************************/

gen age_cat = 0 + ///
	(age >= 45) + ///
	(age >= 50) + ///
	(age >= 55) + ///
	(age >= 60) + ///
	(age >= 65) + ///
	(age >= 70) + ///
	(age >= 75) + ///
	(age >= 80) if age != .


local labs = `"lab def agecatlab 0 "< 45""'
forval i = 1/7 {
	local labs = `"`labs' `i' "< `=`i'*5 + 45'""'
}
local labs = `"`labs' 8 ">= 80""'

`labs'
lab values age_cat agecatlab 

// Checking proper labeling
tabstat age, by(age_cat) s(min max)

// Program to get the stats
cap program drop goftest
prog def goftest, rclass
	// GOF commands
	estat gof
	return local gofpearson : di round(r(chi2), 0.01)
	estat gof, group(10) table
	return local gofhosmer : di round(r(chi2), 0.01)
	
	// AIC
	estimates stats
	tempname mat1 AIC BIC
	mat def `mat1' = r(S)
	mat def `AIC' = `mat1'[1, "AIC"]
	mat def `BIC' = `mat1'[1, "BIC"]
	return local aic : di round(`=`AIC'[1,1]', 0.01)
	return local bic : di round(`=`BIC'[1,1]', 0.01)
end

logit hibp age_cent60 overw
goftest
outreg2 using part3.tex, replace stats(coef se ci tstat) $outregopts ///
	addtext( ///
		Pearson-Chi2, `=r(gofpearson)', ///
		Hosmer-Lemeshow-Chi2, `=r(gofhosmer)', ///
		AIC, `=r(aic)', ///
		BIC, `=r(bic)' ///
	)

logit hibp i.age_cat overw
goftest
outreg2 using part3.tex, append stats(coef se ci tstat) $outregopts ///
	addtext( ///
		Pearson-Chi2, `=r(gofpearson)', ///
		Hosmer-Lemeshow-Chi2, `=r(gofhosmer)', ///
		AIC, `=r(aic)', ///
		BIC, `=r(bic)' ///
	)


/*******************************************************************************
PART 4
*******************************************************************************/

// Higher AIC and BIC
logit hibp age_cent60

// No need of collinearity test, since we know these are mutually exclusive categories.
// TBD

