
import os 

# run the script : 
os.system("python3 make_tree.py")

import math
from util import dagger
from recur_sk import *
import numpy as np



print("\n-------------------SOLOVAY KITAEV ALGORITM TEST----------------------------\n") 
print("1.Press 1 for generating random unitary matrix \n 2.Press 2 to input a unitary matrix")

def deleteI(str):
	newstr= ""
	for i in range(len(str)):
		if(str[i] == "I"):continue
		newstr = newstr + str[i]
	return newstr		

choice = int(input()) 
if(choice == 1):
	print("Enter the accuracy needed to get the approximate circuit")
	e = input()
	# e0 = (1./32)
	U =generate_SU2()
	# n = math.ceil(math.log(ln(1/(e*c*c)/ln(1./(c*c*e0))))/ln 1.5)
	circuit = solovay_kitaev(U,3) # check with length 3  
	Accuracy = trace_norm(multiply(circuit),U) # accuracy obtained 
	print(deleteI(circuit))
	print("ACCURACY",Accuracy)
if(choice == 2): 
	print("Enter the matrix element one after another:\n")
	print("U[0][0]")
	u1 = float(input())
	print("U[0][1]")
	u2 = float(input())	
	print("U[1][0]")
	u3 = float(input())
	print("U[1][1]")
	u4 = float(input())
	U = np.array([[u1,u2],[u3,u4]])
	I = np.array([[1,0],[0,1]])
	if(isEqual(np.matmul(U,dagger(U)),I) == False): 
		print("NOT A UNITARY MATRIX")
	else : 
		print("Enter the accuracy needed to get the approximate circuit")
		e = input()
		# e0 = (1./32)
		# n = math.ceil(math.log(ln(1/(e*c*c)/ln(1./(c*c*e0))))/ln 1.5)
		circuit = solovay_kitaev(U,3) # check with length 3  
		Accuracy = trace_norm(multiply(circuit),U) # accuracy obtained 
		print(deleteI(circuit))
		print("ACCURACY",Accuracy)






