mata:
void sign_changes()
{

M=st_matrix("M")
sign=st_matrix("sign")

g = J(rows(M),1,0)

for (i=1;i<=rows(M);i++) {

if (g[i]==0){


	I = M[i,.]
	S = J(rows(M),cols(M),0)

	I_before=1
	I_after=2	
	

	while (I_before!=I_after) {
	
	I_before = sum(I)
	
	for (c=1;c<=cols(I);c++) {
	
		for (r=1;r<=rows(I);r++) {
	
		S_temp = M:==I[r,c]
		
		S=S+S_temp
		
		}
	
	}
	
	S2 = rowsum(S):>0

	I = select( M:*S2,rowsum(M:*S2):!=0)	
	
	I_after = sum(I)
	
	
	}

	g = g + S2:*sign[i]
	
} 

}

g

st_store(.,"g",g)

}
end


/*


clear

use  "/Users/bferman/Desktop/Matrix.dta"

gen sign = 1 if temp>0.5
replace sign = -1 if temp<0.5

mkmat index*, matrix(M)
mkmat sign, matrix(sign)

gen g=.

mata: sign_changes()
