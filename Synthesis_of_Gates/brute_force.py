import math,cmath 
import numpy as np 
import kdtree 
from Gates import *
from numpy import linalg as LA 

PI = math.acos(-1.0)  
basic_gates_sequences = []
l0 = 13 

#Creating a class to add a payload to the tuple object to be sent to kdtree 
class Sequence(object): # class to hold the record of the matrix in SU(2) and its payload string 
	def __init__(self,r1,c1,r2,c2,payload):
		self.coords = (r1,c1,r2,c2)
		self.payload = payload 
	
	def __len__(self): 
		return len(self.coords)

	def __getitem__(self,i):
		return self.coords[i]
	def __repr__(self):
		return 'Sequence({},{},{},{},{})'.format(self.coords[0],self.coords[1],self.coords[2],self.coords[3],self.payload)

def GenerateUtil(gates,prefix,n,l): 
	
	if(l==0):
		k = len(basic_gates_sequences)  
		# tag the matrix to it
		mat = np.array([[1,0],[0,1]])		
		print(3**13-k) # count down timer  
		for i in range(l0):
			if(prefix[l0-i-1] == "H"): mat = np.matmul(H,mat) 
			if(prefix[l0-i-1] == "T"): mat = np.matmul(T,mat)	
			if(prefix[l0-i-1] == "S"): mat = np.matmul(S,mat)		 					
		
			# Unroll the matrix and tuple it 
		r1 = mat[0][0].real
		c1 = mat[0][0].imag
		r2 = mat[0][1].real
		c2 = mat[0][1].imag		
			# basic_gates_sequences[k] = (prefix,mat)
						
		pts = Sequence(r1,c1,r2,c2,prefix) # a,bz,by,bx 
		basic_gates_sequences.append(pts) 
				
		return 

	for i in range(n):
		newPrefix = prefix + gates[i]	
		GenerateUtil(gates,newPrefix,n,l-1)


def Generate_Sequences(gates,l):
	n=len(gates)
	GenerateUtil(gates,"",n,l)
	

def multiply(seq):
	M = I
	l = len(seq)
	for i in range(l):
		if(seq[l-i-1]=="H"): M = np.matmul(H,M)				
		if(seq[l-i-1]=="T"): M = np.matmul(T,M)
		if(seq[l-i-1]=="S"): M = np.matmul(S,M)	
	return M 
 
def get_tree():
	known_gates = ['H','T','S'] # Hadamard,Identity,Phase-gate 
	Generate_Sequences(known_gates,l0)   
	print(basic_gates_sequences[0]) 

	KD_Tree  = kdtree.create(basic_gates_sequences) 
	return KD_Tree

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

#....................................................................................................................................
import math
from util import dagger
import numpy as np
from NNSearch import *

#tree = get_tree()

print("Loading done...!")

def rotateZ(t):
	Z = np.array([[math.cos(t/2)-1j*math.sin(t/2) , 0],[0 ,math.cos(t/2) + 1j*math.sin(t/2)]])
	return Z 
		
#from random import randint
#for i in range(100):
#	U  = rotateZ(randint(1,100))
#	pause = input()
#	V  = multiply(basic_approx(tree,U,0))
#	print(basic_approx(tree,U,0))	
#	print(np.round(U,decimals=3))
#	print(np.round(V,decimals=3))
#	print("ACCURACY",trace_norm(U,V))
#........................................................................................................................................



