# Implement the control_U2 gate 
# Written by Debarghya Kundu 
# Email : debarghya4469@iitg.ernet.in 

from zy_decomposition import * 
from qiskit import ClassicalRegister,QuantumRegister,QuantumJob
from qiskit import available_backends,execute,register,get_backend
from qiskit.tools.visualization import plot_histogram
from qiskit import QuantumCircuit 

def Rzz(ckt,theta,q,n): # Phase-Shift gate 
	ckt.u1(theta,q[n-1]) # Rzz for last qubit is required 	
	ckt.x(q[n-1])
	ckt.u1(theta,q[n-1])
	ckt.x(q[n-1])
	return ckt

def Rz(ckt,theta,q,n): #rotation along -z axis 
	ckt.u1(theta,q[n-1])
	ckt = Rzz(ckt,-1.0*theta/2.0,q,n)
	return ckt 

def Ry(ckt,theta,q,n): #rotation along -y axis 
	ckt.u3(theta,0,0,q[n-1])
	return ckt 

def CU2(U,ckt,q,n):
	alpha,beta,gama,delta = zy_decompose(U) 
	ckt = Rz(ckt,(delta-beta)/2.0,q,n) # Apply :: C 
	ckt.cx(q[2*n-3],q[n-1]) # ancilla as ctr,tgt
	ckt = Rz(ckt,-1.0*(beta+delta)/2.0,q,n)  # Apply :: B
	ckt = Ry(ckt,-1.0*gama/2.0,q,n) # Apply :: B
	ckt.cx(q[2*n-3],q[n-1]) # ancilla as ctr,tgt
	ckt = Ry(ckt,gama/2.0,q,n) # Apply :: A
	ckt = Rz(ckt,beta,q,n)  # Apply :: A
	ckt.u1(alpha,q[2*n-3]) # phase to ancilla -> ref: N&C 	
	return ckt 

def get_result(ckt):
	shots = 10000
	backend = 'local_qasm_simulator'
	job = execute(ckt,backend=backend,shots=shots)
	results = job.result()
	answer  = results.get_counts()
	plot_histogram(answer)
	
def test_main():
	q = QuantumRegister(3+1)
	c = ClassicalRegister(4)
	a = math.cos(3.0)
	b = math.sin(3.0)
	U = np.array([[a,-1.0*b.conjugate()],[b,a.conjugate()]])
	ckt = QuantumCircuit(q,c)
	# ckt.x(q[3])
	ckt = CU2(U,ckt,q,3)
	for i in range(4):
		ckt.measure(q[i],c[i])

	get_result(ckt)	
		
	
	
	
	


