# Program to get the non-trivial elements 
# Caution : Matrix starts from 1 

def MOD(m,n): # speacial mod 
	if(m % n == 0): return n 
	else : return m % n 

def ntvl(k,l,n): # k,l columns of non-trivial elements of n qubits 
	d = 2**n # total no of qubits possible 
	
	g1 = '' 
	gm = ''
 
	r=1
	while(r <= n):
		u = (n-r+1)
		if(MOD(k,2**u) <= 2**(u-1)): g1 = g1 + '0'
		if(MOD(k,2**u) > 2**(u-1)): g1 = g1 + '1'
		if(MOD(l,2**u) <=  2**(u-1)): gm = gm + '0'
		if(MOD(l,2**u) > 2**(u-1)): gm = gm + '1'
		r = r + 1
	return g1,gm 


#.......................................................DEBUGGED.......................................................................
