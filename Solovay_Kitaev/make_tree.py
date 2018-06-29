import pickle 
import math,cmath 
import numpy as np 
import kdtree 
import dill as pickle  # dill is important for dumping lambda functions
from Gates import *

PI = math.acos(-1.0)  
basic_gates_sequences = []
l0 = 13 # CRASHES AFTER GATE SEQUENCE OF LENGTH 10 
 
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
			#if(prefix[l0-i-1] == "I"): mat = np.matmul(I,mat)	
			if(prefix[l0-i-1] == "S"): mat = np.matmul(S,mat)
			#if(prefix[l0-i-1] == "s"): mat = np.matmul(s,mat)
			#if(prefix[l0-i-1] == "t"): mat = np.matmul(t,mat)		 					
		
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
	
#known_gates = ['H','I','T','S','s','t'] # Hadamard,Identity,Phase-gate and their adjoints..
known_gates = ['H','T','S'] # Hadamard,Identity,Phase-gate 

 

Generate_Sequences(known_gates,l0)   
print(basic_gates_sequences[0]) 

KD_Tree  = kdtree.create(basic_gates_sequences) # create the KDTree 
print(KD_Tree.search_nn([0,0,0,0,1])) # nearest neighbour search check ! 


pickle_out = open("kdtree.pickle","wb") 
pickle.dump(KD_Tree,pickle_out) # dump into pickle file 
pickle_out.close() #close 

#checked ...............................................................................
