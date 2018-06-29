# Extract the two dimensionalise matrix from the TWO-LEVEL-MATRIX 
#[. . . .]                      Hint : Two will be non-diagonal and two diagonal(padding-x,y are same) 
#[. * * .]			
#[. * * .]
#[. . . .]

import numpy as np
ctr = 0
k,l = 0,0
def twoD(U,dim):
	U_tilde = np.identity(2,dtype='complex')
	ctr = 0
	for x in range(dim):
		for y in range(dim):
			if(x != y and np.allclose(U[x][y],0) == False ):
				if(ctr == 0):
					i = [x,y]
					ctr = ctr + 1
				else :  j = [x,y] 
	if(ctr!=0):
		k = i[1]+1 # columns of non-trivial matrix 
		l = j[1]+1 # columns of non-trivial matrix 	
	if(ctr != 0):	
		if(i[1] > j[1]): i,j = j,i # swap 
		#print("asdA",i,j)
		# print("\n\n\n",j[0],j[1],i[0],i[1])		
		x1 = U[j[0]][i[1]] 
		x2 = U[j[0]][j[1]]
		x3 = U[i[0]][i[1]]
		x4 = U[i[0]][j[1]]
		#print(x1,x2,x3,x4)		
		U_tilde = np.array([[x1,x2],[x3,x4]])
		#print(U_tilde)		
#.................................................................
	if(ctr == 0):
		diag = []
		cols = []
		flag = 0
		for d in range(dim):
			if(np.allclose(U[d][d],1) == False):
				flag = flag + 1 
				cols.append(d)
				diag.append(U[d][d])

		if(len(diag) == 2):
			k,l = cols[0],cols[1]
			U_tilde = np.array([[diag[0],0],[0,diag[1]]])
		if(len(diag) == 1): 
			U_tilde = np.identity(2,dtype='complex')

			if(cols[0] == 0): 
				k,l = 1,dim
				U_tilde[0][0] = diag[0]
			if(cols[0] != 0): 
				k,l = 1,cols[0]+1
				U_tilde[1][1] = diag[0] 
		
		if(len(diag) == 0):
			U_tilde = np.identity(2,dtype='complex')
			k,l = 1,dim			
#.....................................................................

	return U_tilde,k,l	
# checked .............................................................................................................
