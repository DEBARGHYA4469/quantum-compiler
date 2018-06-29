# Main routine for synthesis of n-qubit unitary 
# Written by Debarghya Kundu 
# Email : debarghya4469@iitg.ernet.in 

from two_level_decompose import *
from two_dimensionalise import *
from block import * 
import numpy as np 
import math
from random import randint
from graycode import *
from nontrivial import *
from qiskit import ClassicalRegister,QuantumRegister
from qiskit import available_backends,execute,register
from qiskit import get_backend,QuantumJob   

def dagger(U):
	return U.conjugate().transpose()

def get_Vk(U,dim,l):# call only once .........
	VSET = dcmp(U,dim) # set of all two-level-matrices
	Vk = []  
	for V in VSET:
		V = dagger(V)
		Vk.append(V) 
	V_tilde = [] 
	ntv_cols = []  
	for V in Vk : 
		V_tilde.append(twoD(V,dim)[0]) # only matrix...  
		ntv_cols.append([twoD(V,dim)[1],twoD(V,dim)[2]])
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

def bitflip(gi,gj,n):  # checked..................................
	for i in range(n):
		if(gi[i] != gj[i]):
			return i 

def synth(U,dim,nqubit): 
		
	q = QuantumRegister(2*nqubit-1)
	c = ClassicalRegister(2*nqubit-1)
	ckt = QuantumCircuit(q,c)
	

	Vk,V_tilde,ntv_cols = get_Vk(U,dim,dim) # all the two-level matrices
	for block in range(len(Vk)):
		#print("\n--------------------------NEW BLOCK----------------------------------\n")
		#print(np.round(Vk[block],decimals=2))
		#print(np.round(V_tilde[block],decimals=2))		
		
		g1,gm = ntvl(ntv_cols[block][0],ntv_cols[block][1],nqubit) # correct 
		Gcode_sq = [] 
		g1_list = [] 
		gm_list = [] 
		for k in range(nqubit):
			g1_list.append(g1[k])	
			gm_list.append(gm[k])	
		Gcode_sq = graycodes(g1_list,gm_list) # all the gray codes for Vk
		
		m = len(Gcode_sq) # g1----->gm
		#print(g1,gm,Gcode_sq,m) # correct 
		j = 1 
#....................................................DEBUGGED........................................................................
		while(j<=m-2): # all the first m-2 tofollis
			#print("\n.........Each tofollis...................\n")
			bflip=bitflip(Gcode_sq[j-1],Gcode_sq[j],nqubit)
			#print("bitflip:-->",bflip)

			for bi in range(nqubit): 
				if(bflip==bi):continue
				if(Gcode_sq[j][bi]=='0'): 
					ckt.x(q[bi]) # X
			#		print("x is done on:",bi)
			
			if(bflip != nqubit-1):#not last bit flip	                                                          
				ckt = swaper(ckt,bflip,nqubit-1,q)
			#	print("swap",bflip,nqubit-1)

			ckt = nQ_tofolli(ckt,nqubit,q)	
			#print("Toffoli")					
			
			if(bflip != nqubit-1):#not last bit flip	                                                          
				ckt = swaper(ckt,bflip,nqubit-1,q)
			#	print("swap",bflip,nqubit-1)

			for bi in range(nqubit): 
				if(bflip==bi):continue
				if(Gcode_sq[j][bi]=='0'):
					ckt.x(q[bi]) # X
			#		print("x is done on:",bi)
			j = j + 1 

		# Apply controlled-U for change from gm-1 to gm
#---------------------------------------------------------------------------------------------------------------------		
		bflip = bitflip(Gcode_sq[m-2],Gcode_sq[m-1],nqubit)
		for bi in range(nqubit): 
			if(bflip==bi):continue
			if(Gcode_sq[j][bi]=='0'): 
				ckt.x(q[bi]) # X

		if(bflip != nqubit-1):#not last bit flip	                                                          
			ckt = swaper(ckt,bflip,nqubit-1,q)

		ckt = Control_U(ckt,V_tilde[block],q,nqubit)		

		if(bflip != nqubit-1):#not last bit flip	                                                          
			ckt = swaper(ckt,bflip,nqubit-1,q)

		for bi in range(nqubit): 
			if(bflip==bi):continue
			if(Gcode_sq[j][bi]=='0'):
				ckt.x(q[bi]) # X
#----------------------------------------------------------------------------------------------------------------------
		# Reverse the states ......................................

		j = m-2 
		while(j >= 1): # all the first m-2 tofollis
			bflip=bitflip(Gcode_sq[j-1],Gcode_sq[j],nqubit)
			#print("bitflip:",j,"-->",bflip)
			for bi in range(nqubit): 
				if(bflip==bi):continue
				if(Gcode_sq[j][bi]=='0'): 
					ckt.x(q[bi]) # X
			
			if(bflip != nqubit-1):#not last bit flip	                                                          
				ckt = swaper(ckt,bflip,nqubit-1,q)
			
			ckt = nQ_tofolli(ckt,nqubit,q)						
			
			if(bflip != nqubit-1):#not last bit flip	                                                          
				ckt = swaper(ckt,bflip,nqubit-1,q)

			for bi in range(nqubit): 
				if(bflip==bi):continue
				if(Gcode_sq[j][bi]=='0'):
					ckt.x(q[bi]) # X
			j = j - 1

	for i in range(2*nqubit-1):
		ckt.measure(q[i],c[i])
	get_result(ckt)

def kronecker(str): # find the kroncker product 
	z=0	
	if(str[0]=='0'):
		z=np.array([[1],[0]])
	if(str[0]=='1'):
		z=np.array([[0],[1]])
	for i in range(len(str)-1):
		if(str[i+1]=='0'): z = np.kron(z,np.array([[1],[0]]))
		if(str[i+1]=='1'): z = np.kron(z,np.array([[0],[1]])) 
	print(z)	

#u1 = np.identity(2)
#u2 = np.array([[1,1],[1,-1]])/math.sqrt(2.0)
#u= np.kron(u1,u2)
#u= np.kron(u,u1)
#u = np.array([[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,1,0],[0,0,0,0,0,1,0,0],[0,0,0,0,1,0,0,0],[0,0,0,1,0,0,0,0],[0,0,1,0,0,0,0,0],[0,1,0,0,0,0,0,0],[1,0,0,0,0,0,0,0]])
#u = np.identity(16)
#for i in range(16):
#	u[i][i]=0
#	u[i][15-i]=1
#print(u)
#kk = np.array([[1],[0],[0],[0],[0],[0],[0],[0]])
#print(np.matmul(u,kk))
#synth(u,16,4)


# test .....................2018 


#Vk,V_tilde,ntv_cols = get_Vk(u,8,8) # all the two-level matrices

#for V in range(len(Vk)):
#	print(np.round(Vk[V],decimals=2),"\n")
#	print(np.round(V_tilde[V],decimals=2),"\n")
#	print("gcodes-->",ntv_cols[V][0],ntv_cols[V][1])

#print(np.allclose(np.identity(8),Vk[len(Vk)-1]))






