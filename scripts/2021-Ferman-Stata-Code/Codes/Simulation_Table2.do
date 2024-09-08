
* Note: this code generates the results for panels A to C. To generate the other panels, we just have to change the distribution of y for the treated. 

clear * 

cd "/Users/bferman/Dropbox/Projects/Matching with a finite Nt and many Nc/Codes and data/Table 2"

do mata/mata_sign_changes.do


forvalues round=1(1)50 {

foreach N_1 in 5 10 25 50  { 

foreach N_0 in   500  {

foreach power in 0  { 

foreach M in  1 4 10  {  

foreach k in 1 8 999 {

timer clear
timer on 1

set seed `round'


clear
clear matrix
set matsize 200
mat R=J(100,5,.)


forvalues row=1(1)100 {

clear

qui {



************** Sample from Prova Brasil distribution

clear
local N = `N_1'+`N_0'
set obs `N'

gen T = _n<=`N_1'

gen x = rnormal()

gen y = rnormal()

if `k'!=999 {

gen p = normal(x)
gen TE = (invchi2(`k',p) - `k')/sqrt(2*`k')

}

if `k'==999 {

gen TE = x

}


replace y = y + TE + `power'/10 if T==1



************** Run matching estimator and calculate AI standard error.


nnmatch y T x, tc(att) m(`M') population robust(2)  keep(temp/temp`round'_`row'_`N_1'_`N_0'_`power'_`M'_k_`k', replace)

gen n=_n

mat B=e(b)
mat V=e(V)

mat R[`row',1]=B[1,1]   // Bias without bias-correction
mat R[`row',2]=B[1,1]/sqrt(V[1,1]) // Test without bias-correction


************** Organize data on the nearest neighbors.

use "temp/temp`round'_`row'_`N_1'_`N_0'_`power'_`M'_k_`k'.dta", clear

cap erase "temp/temp`round'_`row'_`N_1'_`N_0'_`power'_`M'_k_`k'.dta"


expand 2 if _n==1 | id[_n]!=id[_n-1]

drop T
gen T=_n>`N_1'*`M'

replace index=id if _n>`N_1'*`M'
replace y=y_0 if _n<=`N_1'*`M'

keep y id index T
gsort id -T  index



***** RI - sign changes

bysort id: egen temp = mean(y) if T==0 
gsort -T id index

egen y_hat = mean(temp), by(id)
drop temp
replace y_hat=. if T==0

gen tau_temp=y-y_hat if T==1

egen tau = mean(tau_temp), by(id)

gsort id -T

drop if T==1
drop y_hat tau_temp T y 

bysort id: gen M=_n

reshape wide index, i(id) j(M)

summ tau 
local tau=0.999 * (r(mean)/r(sd))^2

mkmat index*, matrix(M)

	mat B=J(200,2,.)
	
	forvalues B=1(1)200 {		
		
		gen temp=runiform()
		gen sign = 1 if temp>0.5
		replace sign = -1 if temp<=0.5
		mkmat sign, matrix(sign)

		gen g = .

		mata: sign_changes()
	
		gen tau_p = tau*g
		
		summ tau_p	
		mat B[`B',1]=(r(mean)/r(sd))^2	
		
		drop  temp g tau_p sign
	
	}
	
	
	
	preserve
	
	clear
	set more off
	svmat B
	
	gen p = `tau'<B1
	
	summ p
	mat R[`row',4] = r(mean)	
	
	restore
	
	drop tau



}

di "row=`row' power=`power' N1=`N_1' N0=`N_0' Round=`round' M=`M' ,  Chi treatment k =`k'"

}



clear
set more off 
svmat R

rename R1 bias
la var bias "Bias without correction"

rename R2 t_AI 
la var t_AI "t-stat AI without correction"

rename R3 p_RI_p
la var p_RI_p "p-value RI permutation test"

rename R4 p_RI_s
la var p_RI_s "p-value RI sign changes"

rename R5 p_RI_p_st
la var p_RI_s "p-value RI permutation test standardized"



gen reject5_AI =  t_AI<-1.96 | t_AI>1.96 if t_AI!=.

gen reject5_RI_p =  p_RI_p<0.05 if p_RI_p!=.

gen reject5_RI_p_st =  p_RI_p_st<0.05 if p_RI_p_st!=.


gen reject5_RI_s =  p_RI_s<0.05 if p_RI_s!=.



gen N1=`N_1'
gen N0=`N_0'
gen R=`round'
gen power=`power'

gen M=`M'


timer off 1
timer list 1
gen time = r(t1)/60


aorder



saveold Results/N1_`N_1'_N0_`N_0'_power_`power'_M`M'_chi_treatment`k'_`round', replace version(12)  

}
}
}
}
}
}
	
	
	
	
	
	
	
	
	





