
clear * 

cd "/Users/bferman/Dropbox/Projects/Matching with a finite Nt and many Nc/Codes and data/Table 2"

gen k=.

forvalues round=1(1)50 {

foreach N_1 in 5 10 25 50   { // 5 10 25 50

foreach N_0 in   1000  { //  50 500

foreach power in  0  { // 20 40 60

foreach M in  1 4 10  {  // 1 2 4 10

foreach k in  1 8  999 {

cap append using  "Results/N1_`N_1'_N0_`N_0'_power_`power'_M`M'_chi_treatment`k'_`round'"

replace k = `k' if k==.
}
}
}
}
}
}

bysort M: tab N1 k

	
mat R = J(20,8,.)

prog define Store, rclass
args row col M k
 {
 
forvalues n = 1(1)4 {
summ rej_AI  if M==`M' & N==`n' & k==`k'
mat R[`row'+`n'-1,`col']=r(mean)

summ rej_RI  if M==`M' & N==`n' & k==`k'
mat R[`row'+`n'-1,`col'+1]=r(mean)
}

			
}
end	
	
	
gen rej_AI = abs(t_AI)>invnormal(0.95)
gen rej_RI = p_RI_s<=0.1


egen N = group(N1)

Store 1 1 1 999

Store 1 4 4 999

Store 1 7 10 999
	
	
Store 6 1 1 8

Store 6 4 4 8

Store 6 7 10 8
		
	
Store 11 1 1 1

Store 11 4 4 1

Store 11 7 10 1
		
	
	
mat li R
	
	
