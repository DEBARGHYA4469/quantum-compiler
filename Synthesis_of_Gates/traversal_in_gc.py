# Traversal in Gray code 
import numpy as np
import math 
from two_level_decompose import *
from two_dimensionalise import *

#H2 = np.array([[1,1],[1.,-1]])/math.sqrt(2)
#H4 = np.kron(H2,H2) # kronecker product

def dagger(U):
	return U.conjugate().transpose()

def get_matrix(U,dim,l):# call only once .........
	VSET = dcmp(U,dim) # set of all two-level-matrices
	Z = np.identity(dim)
	Vk = []  
	for V in VSET:
		V = dagger(V)
		Vk.append(V) 
		Z=np.matmul(Z,V)

# ........Vk contains the set of all two level matrix multiplying to U .............................

	V_tilde = [] 
	ntv_cols = []  
	for V in Vk : 
		V_tilde.append(twoD(V,4)[0]) # only matrix...  
		ntv_cols.append([twoD(V,4)[1],twoD(V,4)[2]])
	return Vk,V_tilde,ntv_cols 


#------------------------------------------------------TRAVERSAL ---------------------------------------------------------
#
# U = V1 x V2 x V3 x V4 x .....................x Vk  
#
# <---------   V1  ------------------------------> <--------------------V2-------------------------->               .......Vk........>
# g11 g21 g31 g41 ............gm-1,1 C~U gm-1,1 ........g12 g22............gm...........................................
# [           Call this part a BLOCK 1           ]
#
# g(i,j)'s are represented by n qubit tofolli gate. Stimulate this by 3-Q Tofollis 
# Controlled-U by SWAP and Tofolli operations ............
#
# 1. Get the set of matrices Vi's
# 2. Get the two-dimensional form of it 
# 3. get the set of non-trivial qubits g1 and gm from Vk
# 4. Get the Grey code sequence 
# 5. Traverse g1 to gm-1 
# 6. Apply controlled-U 
# 7. Traverse gm-1 to gm 
# 8. Go to step 1 unless set replenishes
#----------------------------------------------------TRAVERSAL------------------------------------------------------------

# Code for the main part 
from qiskit import ClassicalRegister,QuantumRegister
from qiskit import available_backends,execute,register
from qiskit import get_backend,QuantumJob   
from nontrivial import *
from graycode import *
 
def traverse(U,dim,nqubit):
	Vk,Vtilde,ntv_cols = get_matrix(U,dim,dim) 
	# Vk contains 2-level matrix,Vtilde : 2x2 
	# ntv_cols are non-trivial cols for gray code sequences
	
	for block in range(len(Vk)): # each blocks   
		g1,gm = ntvl(ntv_cols[block][0],ntv_cols[block][1],nqubit)
		Gcode_sq = []
		g1_list = []
		gm_list = []
		for k in range(nqubit):
			g1_list.append(g1[k])
			gm_list.append(gm[k])
		Gcode_sq = graycodes(g1_list,gm_list) # all the codes 
       # CODES generated ..................
		# TRAVERSE FROM g1 to gm-1
						
					
		#for i in range(len(G_code_sq)):
		#	print(G_code_sq)
		print("hello",g1,gm,Gcode_sq)		





