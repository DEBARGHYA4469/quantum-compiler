from qiskit import ClassicalRegister,QuantumRegister,QuantumJob
from qiskit import available_backends,execute,register,get_backend
from qiskit.tools.visualization import plot_histogram
from qiskit import QuantumCircuit 
from zy_decomposition import *
from control_U2 import *
#.................................................................
def swaper(ckt,ctr,tgt,q): # swap control and target bit
	ckt.cx(q[ctr],q[tgt])	
	ckt.h(q[ctr])
	ckt.h(q[tgt])	
	ckt.cx(q[ctr],q[tgt])
	ckt.h(q[ctr])
	ckt.h(q[tgt])
	ckt.cx(q[ctr],q[tgt])
	return ckt 
#.................................................................

#.........................................................n-Q Toffolli...................................................
def nQ_tofolli(ckt,n,q): # always take first n-1 as ctr,n as tgt  
	i = 0 	
	anc = n # ancillas are from n to 2*n-1 
	while(i < n-2):
		if(i==0):		
			ckt.ccx(q[0],q[1],q[anc])
		else : 
			ckt.ccx(q[i+1],q[anc+i-1],q[anc+i])
		i = i + 1
	i = i - 1
	# Perform the controlled operations 
	ckt.cx(q[anc+i],q[n-1]) # Controlled-U operation
	
	# Free the ancillas....
	while(i >= 0):
		if(i==0):		
			ckt.ccx(q[0],q[1],q[anc])
		else : 
			ckt.ccx(q[i+1],q[anc+i-1],q[anc+i])
		i = i - 1 
	return ckt 
#.......................................................n-Q Tofolli ......................................................	

#................................................Controlled-U operations .................................................
def Control_U(ckt,V_tilde,q,n): # n is no of qubits 		
	i = 0 	
	anc = n # ancillas are from n to 2*n-1 
	while(i < n-2):
		if(i==0):		
			ckt.ccx(q[0],q[1],q[anc])
		else : 
			ckt.ccx(q[i+1],q[anc+i-1],q[anc+i])
		i = i + 1
	i = i - 1
	# Perform the controlled operations 
	ckt = CU2(V_tilde,ckt,q,n) #  Controlled-U operation
	
	# Free the ancillas....
	while(i >= 0):
		if(i==0):		
			ckt.ccx(q[0],q[1],q[anc])
		else : 
			ckt.ccx(q[i+1],q[anc+i-1],q[anc+i])
		i = i - 1 
	return ckt 

#................................................Controlled-U operations .................................................


def test_block():	
	q   =  QuantumRegister(5+4) # Quantum bits 
	c   =  ClassicalRegister(9) #Classical Register
	U   =  np.array([[1,1],[1,-1]])/math.sqrt(2.0)
	ckt = QuantumCircuit(q,c)
	for j in range(5):
		ckt.x(q[j])
	ckt.x(q[4])
	ckt = Control_U(ckt,U,q,5)
	
	for i in range(9):
		ckt.measure(q[i],c[i])
	get_result(ckt)



