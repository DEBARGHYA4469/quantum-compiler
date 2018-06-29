from two_level_decompose import * 
from random import randint 
import numpy as np 
import math

def e_i(x): # complexexponent
	return (math.cos(x) + 1j*math.sin(x))  

def unitary(alpha,beta,gama,delta):               # Generate unitary matrix 
	mat = np.array([[1,0],[0,1]])
	a1 = e_i(-1.0*beta/2.0)
	a2 = e_i(beta/2.0)
	b1 = e_i(delta/2.0)
	b2 = e_i(-1.0*delta/2.0)
	c = e_i(alpha)
	A_mat = np.array([[a1,0],[0,a2]])
	B_mat = np.array([ [math.cos(gama/2.0),math.sin(gama/2.0)*-1], [math.sin(gama/2.0),math.cos(gama/2.0)] ])
	C_mat = np.array([[b2,0],[0,b1]])
	mat = np.matmul(A_mat,np.matmul(B_mat,np.matmul(C_mat,mat)))*c
	return mat

for i in range(100):
	Z = np.identity(4,dtype='complex')
	u1 = unitary(randint(1,100),randint(1,1000),randint(1,100),randint(1,100))
	u2 = unitary(randint(1,100),randint(1,1000),randint(1,100),randint(1,100))
	U = np.kron(u1,u2)
	VSET = dcmp(U,4)
	Z = U
	for V in VSET:
		Z = np.matmul(V,Z)
	if(np.allclose(Z,np.identity(4)) == True):print("Success!")
	else :print("Failure",len(VSET))

