# IMPLEMENTATION OF GROVER'S SEARCH ALGORITHM TO SEARCH A DATABASE CONTAINING N(here N=4) ELEMENTS BY QUERING ONLY O(sqrt(M/N)) 
# WHERE M IS THE NO OF SOLUTIONS. HERE M < N/2 .

#from connection import * 
import numpy as np  
from qiskit import QuantumCircuit,ClassicalRegister,QuantumRegister,QuantumJob
from qiskit import available_backends,execute,register,get_backend
from qiskit.tools.visualization import plot_histogram

#backend = get_qc() # get the quantum computer with least jobs queued up 

def Swap(grov,x1,x2,q): # Swap lines x1 and x2 
	grov.cx(q[x1],q[x2])
	grov.h(q[x1])
	grov.h(q[x2])
	grov.cx(q[x1],q[x2])
	grov.h(q[x1])
	grov.h(q[x2])
	grov.cx(q[x1],q[x2])
	return grov 
	
def Tofolli(grov,q): # implementation of Tofolli Gate 
	grov.h(q[2])
	grov.cx(q[1],q[2])
	grov.tdg(q[2])
	grov.cx(q[0],q[2])
	grov.t(q[2])
	grov.cx(q[1],q[2])
	grov.tdg(q[2])
	grov.cx(q[0],q[2])
	grov.t(q[1])
	grov.t(q[2])
	grov.h(q[2])
	grov = Swap(grov,1,2,q)
	grov.cx(q[0],q[2])
	grov.t(q[0])
	grov.tdg(q[2])
	grov.cx(q[0],q[2])
	grov = Swap(grov,1,2,q)
	return grov 	

def phase(grov,q):
	grov.x(q[0])
	grov.x(q[1])
	grov.h(q[1])
	grov.cx(q[0],q[1])
	grov.h(q[1])
	grov.x(q[0])
	grov.x(q[1])
	return grov							
def Hadamard_transform(grov,q):
	grov.h(q[0])
	grov.h(q[1])
	return grov
	
def grover_operator(grov,q): #grover_operator
	grov = Tofolli(grov,q)
	grov = Hadamard_transform(grov,q)
	grov = phase(grov,q)	
	grov = Hadamard_transform(grov,q)
	return grov 
	

N=2 # =log(N) where database of size N is searched 
M=1 # no of solutions 

q = QuantumRegister(N+1)
c = ClassicalRegister(N) 

circuit_Name = "Grover's Search!"
grov = QuantumCircuit(q,c)

grov.h(q[0])
grov.h(q[1])# Prepare psi : equally probable state
grov.x(q[2])
grov.h(q[2])

for k in range(int((N/M)**0.5)): # no of rotation required 
	grov = grover_operator(grov,q)

grov.h(q[2])
grov.measure(q[0],c[0])
grov.measure(q[1],c[1])

backend = "local_qasm_simulator"
shots = 1000
job = execute(grov,backend=backend,shots=shots)
results=job.result()
answer = results.get_counts()

plot_histogram(answer)
