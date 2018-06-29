from random import randint 
import numpy as np 
from zy_decomposition import * 
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

def test_zy_decompose():
	for i in range(1000):
		U = unitary(randint(0,100), randint(0,100) , randint(0,100) , randint(0,100) ) 
		#if(np.allclose(np.matmul(U,U.transpose().conjugate()),np.identity(2)) == True): print("Yeah:",i)
		#else:break 
		a,b,g,d = zy_decompose(U,0) 
		if(np.allclose((math.cos(a)+1j*math.sin(a))*np.matmul(rotateZ(b),np.matmul(rotateY(g),rotateZ(d))),U) == True):
			print("Yeah:",i)
		else : 
			return U
test_zy_decompose()	
