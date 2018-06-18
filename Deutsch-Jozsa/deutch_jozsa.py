#from connection import * 
import numpy as np 
import random 
import matplotlib.pyplot as plt 
from qiskit import QuantumCircuit,ClassicalRegister,QuantumRegister,QuantumJob
from qiskit import available_backends,execute,register,get_backend
from qiskit.tools.visualization import plot_histogram,circuit_drawer

#backend = get_qc() # get the quantum computer with least jobs queued up 
def flip_x(arr,binlength,l):
	if(l == 0):
		newstr = "" 
		for b in range(binlength):newstr = newstr + str(arr[b])			
		print("f(",newstr,")=1")
	else:
		 
		for pv in range(2):
			acopy = []
			for i in range(binlength): acopy.append(arr[i]) 
			for k in range(binlength):
				if(acopy[k] == -1): 
					acopy[k] = pv 
					break
			flip_x(acopy,binlength,l-1)

def mark_ones(j,shift,pos_ones,arr):
	for i in range(j): 		
		arr[pos_ones[shift+i]] = 1
	return arr 

def rec_util(j,binlength,pos_ones,shift):
	if(shift+j-1 >= len(pos_ones)):
		return
	arr = []
	for i in range(binlength): arr.append(-1) 
	arr = mark_ones(j,shift,pos_ones,arr)
	flip_x(arr,binlength,binlength-j)
	rec_util(j,binlength,pos_ones,shift+1)

def ones_in_binary(r): # Decimal to binary conversion
	num = r
	bry = 0 
	m = [] 
	j = 0 
	f=1
	while(num>0):         
		bry = bry + (num%2)*f
		if(num%2 == 1):
			m.append(j) 
		num = num//2
		f=f*10
		j=j+1

	for i in range(len(m)):
		m[i] = j -1- m[i]
	pos_ones = [] 
	for i in range(len(m)):
		pos_ones.append(m[len(m)-i-1])
	print(pos_ones)
	return pos_ones,j
    
def print_fx(r,n): 
	j=1 
	pos_ones,binlength = ones_in_binary(r)#position of the ones and total length of binary string
	binlength = n
	while(j <= len(pos_ones)): #only odd results in flip by CNOT
		shift = 0 # to traverse the pos_ones
		rec_util(j,binlength,pos_ones,shift)
		j=j+2 
print("\n................DEUTSCH JOZSA ALGORITHM....................\n")

print("\nInput the size of the domain of the function..\n")
n = int(input())

print("\nInput 0 if you want a constant function,1 otherwise\n")

oracle_Type  = int(input()) # constant   or balanced 
if(oracle_Type == 0): # constant function 
	print("What is the constant value 0/1 ?") 
	oracle_val = int(input())  
else: 
	print("The function is balanced...Generating a random f(x)")
	fx_generator = np.random.randint(1,2**n) # this variable creates the quantum circuit of the f(x)
	print("\nYour Balanced Function is .........................")
	print(fx_generator)
	print_fx(fx_generator,n) # print the function 
	print("\nn Rest all values are 0s")

print("\nTHIS FUNCTION WILL REMAIN HIDDEN AND ONLY BE QUERIED TO ITS EQUIVALENT QUANTUM CIRCUIT")

#Print the quantum circuit for your random balanced oracle

file = open("oracle.qasm","a") # append mode 
file.write("\n")

# Interfacing with quantum computer 
qr = QuantumRegister(n+1) # query string 
cr = ClassicalRegister(n) # first register

circuitName = "Deutsch_Joszsa"
djCircuit = QuantumCircuit(qr,cr)
for i in range(n+1):
	file.write(" qubit q"+str(i)+"\n")

#superposition of all input queries by hadamard transform 
for i in range(n):
	file.write(" h q"+str(i)+"\n")
	djCircuit.h(qr[i])

#FLip the second register and apply hadamard 
djCircuit.x(qr[n])
file.write(" X q"+str(n)+"\n")
djCircuit.h(qr[n])
file.write(" h q"+str(n)+"\n")

djCircuit.barrier()
# IMPLEMENTING THE ORACLE BUT ITS VALUES ARE HIDDEN 
if(oracle_Type == 0): # if constant
	if(oracle_val==1):
		djCircuit.x(qr[n])
		file.write(" X q"+str(n)+"\n")
	else:
		djCircuit.iden(qr[n])
else:
	for i in range(n):
		if(fx_generator & (1 << i)):
			djCircuit.cx(qr[i],qr[n])
			file.write(" cnot q"+str(i)+",q"+str(n)+"\n")		

djCircuit.barrier()

#apply hadamard after querying the oracle 
for i in range(n):
	djCircuit.h(qr[i])
	file.write(" h q"+str(i)+"\n")
#measurements
for i in range(n):
	djCircuit.measure(qr[i],cr[i])
	file.write(" measure q"+str(i)+"\n")
#circuit_drawer(djCircuit)
file.close()



backend = "local_qasm_simulator"
shots = 1000
job = execute(djCircuit,backend=backend,shots=shots)
results=job.result()
answer = results.get_counts()

plot_histogram(answer)
