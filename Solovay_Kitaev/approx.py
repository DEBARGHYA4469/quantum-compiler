import pickle 
import numpy as np
from numpy import linalg as LA 
import math,cmath
from random import randint  
from Gates import *
from util import dagger

PI = math.acos(-1.0)

def trace_of(U):
	return (abs(U[0][0] + U[1][1]))

def square_root_matrix(U):
	W,V = LA.eig(U)
	D = np.array([[cmath.sqrt(W[0]),0],[0,cmath.sqrt(W[1])]]) # diagonal matrix 
	V_dagger = V.transpose().conjugate() # compute conjugate-transpose 	
	sqR = np.matmul(V,np.matmul(D,V_dagger))
	return sqR 

def trace_norm(U,V):
	M = U-V 
	M_dagger = (U-V).transpose().conjugate()
	norm = square_root_matrix(np.matmul(M_dagger,M))
	return trace_of(norm)
	

def multiply(seq):
	M = I
	l = len(seq)
	for i in range(l):
		if(seq[l-i-1]=="H"): M = np.matmul(H,M)				
		if(seq[l-i-1]=="T"): M = np.matmul(T,M)
		if(seq[l-i-1]=="S"): M = np.matmul(S,M)
		if(seq[l-i-1]=="s"): M = np.matmul(s,M)
		if(seq[l-i-1]=="t"): M = np.matmul(t,M)
		if(seq[l-i-1]=="h"): M = np.matmul(h,M)			
	return M 

def e_i(x): # complexexponent
	return (math.cos(x) + 1j*math.sin(x))  


def unitary(alpha,beta,gama,delta):
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

#checked ..................................................................











