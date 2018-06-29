# Decomposition of a matrix into two level forms 
# Given a unitary matrix U decompose it into U = V1V2...Vk

import numpy as np 
import math 

# VSET = [] # set of Vi's.........................................................................add

def star(x): 
	return x.conjugate() 

def Iden(l):
	return np.identity(l,dtype='complex')

def level_up(U,dim,l): # level up the matrix(l) to dim = dim 
	Y = Iden(dim)
	pad = dim-l
	for i in range(l):
		for j in range(l):
			Y[i+pad][j+pad] = U[i][j]
	return Y 	
	
def level_down(U,dim,l): # one-level-down 
	Y = Iden(l-1) # lower one-level down 
	i = 1 
	j = 1 
	while(i < l):
		j=1
		while(j < l):
			Y[i-1][j-1] = U[i][j]			
			j=j+1
		i=i+1
	return Y

def norm(a,b):
	return math.sqrt((abs(a)**2 + abs(b)**2))
#......................................................DEBUG ENDS >........................................................
def TLD(VSET,X,dim,l): # Gate U , dimension U ,level = l 
#	if(np.allclose(X,Iden(dim))== True):#matrix converts into I 
#		return VSET 
	if(l == 1):
		#print(np.allclose(X,Iden(dim)))
		Y = Iden(dim)
		Y[dim-1][dim-1] = star(X[dim-1][dim-1]) 
		VSET.append(Y) # last crucial step
		return VSET
	i = 1 
	pad = dim-l	
	while(i <= l-1): # for one-level down
		if(np.allclose(X[i+pad][0+pad],0) == True):
			i = i + 1
			continue 
		V = Iden(dim) # identity of dim=l
		a = X[0+pad][0+pad]
		b = X[i+pad][0+pad]	
		V[0+pad][0+pad] = star(a)/norm(a,b) 
		V[i+pad][0+pad] = b/norm(a,b)
		V[0+pad][i+pad] = star(b)/norm(a,b)
		V[i+pad][i+pad] = -1.0*a/norm(a,b)
		
		X = np.matmul(V,X)
		
		#print(X)
		#print(V)			
		VSET.append(V)
		i = i + 1
	#print("LEVEL:",l)
	#X = level_down(X,dim,l)
	#print(X)
	return TLD(VSET,X,dim,l-1)
def dcmp(U,dim):
	VSET = [] # reset VSET 
	return TLD(VSET,U,dim,dim)
# ..............................................CHECKED...............................................................
#X = U 
#print("Original MATRIX U \n\n",X,"\n")
#VSET = []
# print(TLD(X,4,4))
