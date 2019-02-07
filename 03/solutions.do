clear all
set more off
set trace off


cd /home/vegayon/Dropbox/usc/clases/2019_PM511B/assignments

// An awful way to read-in data
input str40 Snoring HeartD1 HeartD0 PropYes Linear Logit Probit
	"Never" 24 1355 0.017 0.017 0.021 0.020
	"Occasional" 35 603 0.055 0.057 0.044 0.046
	"Nearly every night" 21 192 0.099 0.096 0.093 0.095
	"Every night" 30 224 0.118 0.116 0.132 0.131
	end
	
// Turning the variable into a labeled variable
gen SnoringFactor = ///
	2*(Sno == "Occasional") + ///
	4*(Sno == "Nearly every night") + ///
	5*(Sno == "Every night")
	

lab def snoringlab 0 Never 2 Occasional 4 "Nearly every night" ///
	5 "Every night"
lab values SnoringFactor snoringlab

// Are we getting the right thing?
table Snoring SnoringF

// Droping the wrong
drop Snoring
rename SnoringF snoring
label var snoring Snoring

// Reshaping the data
reshape long HeartD, i(snoring) j(disease)
rename HeartD nobs
label var disease "Heart Disease"

// Checking
table snoring dise [fw=nobs]

/*******************************************************************************
Question 1
*******************************************************************************/

// Program to capture the coefficients in a mata matrix
prog def getcoefs
	args globname
	
	// Getting from the "ereturn list"
	mata: `globname' = st_matrix("e(b)")
end

// Program to predict disease for individuals with 0 or 5 in s
prog def predictedprobs, rclass
	preserve
	quietly {
		drop _all
		set obs 2
		gen snoring = 0 if _n == 1
		replace snoring = 5 if _n == 2
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

// i. Linear regression (use the REG command and compare estimates)
glm disease snoring [fw=nobs], family(gaussian) link(identity)
predictedprobs // Calculating predicted probabilities and RR
outreg2 using 03/question1.tex, replace stats(coef se ci) tex(fragment) ///
	addtext(Family, Gaussian, Link, Identity, $ E[y|X] $, $ X\hat\beta $, ///
	$ P[Disease=1|X=0] $, `=r(p0)', $ P[Disease=1|X=5] $, `=r(p5)', ///
	Risk Ratio, `=r(rr)') label dec(2)
getcoefs mat0

reg disease snoring [fw=nobs]
getcoefs mat1

// Getting the same?
mata: mat0 == mat1


// ii. Binary regression (random component = binomial family)
//     1. Linear probability model (link = identity)
glm disease snoring [fw=nobs], family(binomial) link(identity)
predictedprobs // Calculating predicted probabilities and RR
outreg2 using 03/question1.tex, append stats(coef se ci) tex(fragment) ///
	addtext(Family, Binomail, Link, Identity, $ E[y|X] $, $ X\hat\beta $, ///
	$ P[Disease=1|X=0] $, `=r(p0)', $ P[Disease=1|X=5] $, `=r(p5)', ///
	Risk Ratio, `=r(rr)') label dec(2)


//     2. Logit probability model (link = logit)
glm disease snoring [fw=nobs], family(binomial) link(logit)
predictedprobs // Calculating predicted probabilities and RR
outreg2 using 03/question1.tex, append stats(coef se ci) tex(fragment) ///
	addtext(Family, Binomial, Link, Logit, $ E[y|X] $, $\mbox{Logit}^{-1}[X\hat\beta] $, ///
	$ P[Disease=1|X=0] $, `=r(p0)', $ P[Disease=1|X=5] $, `=r(p5)', ///
	Risk Ratio, `=r(rr)') label dec(2)

//     3. Probit probability model (link = probit)
glm disease snoring [fw=nobs], family(binomial) link(probit)
predictedprobs // Calculating predicted probabilities and RR
outreg2 using 03/question1.tex, append stats(coef se ci) tex(fragment) ///
	addtext(Family, Binomial, Link, Probit, $ E[y|X] $, $\Phi^{-1}[X\hat\beta] $, ///
	$ P[Disease=1|X=0] $, `=r(p0)', $ P[Disease=1|X=5] $, `=r(p5)', ///
	Risk Ratio, `=r(rr)') label dec(2)


	

/*******************************************************************************
Question 2
*******************************************************************************/

// We will need outreg2
// ssc install outreg2

cd /home/vegayon/Dropbox/usc/clases/2019_PM511B/assignments/

// Loading the data and generating the binary variable hibp
use 01/vitals, clear

// Previous variables created
gen hibp = sbp > 140 if sbp != .
gen bmi = 703 * weight / height^2 if (weight != . & sbp != .)
gen overweight = bmi >= 25 if bmi != .

// Part a.
probit hibp age
outreg2 using 03/question2.tex, replace stats(coef se ci) tex(fragment)

probit hibp overweight
outreg2 using 03/question2.tex, append stats(coef se ci) tex(fragment)

// Part b.
probit hibp age overweight
outreg2 using 03/question2.tex, append stats(coef se ci) tex(fragment)

// Part e.
// i. A person aged 70 who is overweight 
// ii. A person aged 50 who is not overweight

preserve
drop _all
set obs 2
gen overweight = 1
replace overweight = 0 if _n == 2

gen age = 70
replace age = 50 if _n == 2
predict prob_hibp

list

restore
