from synthesis import * 
import numpy as np 
import math 

PI = math.acos(-1.0)

H = np.array([[1 , 1],[1 , -1]])/math.sqrt(2) # Hadamard  
T = np.array([[math.cos(PI/8)-1j*math.sin(PI/8) , 0],[0 , math.cos(PI/8)+1j*math.sin(PI/8)]]) # Phase Gate 
I = np.array([[1,0],[0,1]]) # Identity matrix	
S = np.matmul(T,T) # S gates 
X = np.array([[0,1],[1,0]])

from util import *
from random import random  
# test one_qubit gates  
def test_U():
	U = unitary(random(),random(),random(),random())
	handle_1Qgates(U,2,1)
def test_H():
	handle_1Qgates(H,2,1)
def test_T():
	handle_1Qgates(T,2,1)

#test two_qubit gates 
#def test_HH(): # should return 00,11,10,01 
#	u = np.kron(H,H)	
#	synth(u,4,2)
#def test_XX(): # should return 11
#def test_HX(): # should return 01,11
#def test_TT():   

# test three_qubit gates 
def test_UUU():
	u = unitary(random(),random(),random(),random())
	u1=np.kron(u,u)
	u1=np.kron(u1,u)
	synth(u1,8,3)

def test_HHH(): # should return 000,101,100,001,... 
	u = np.kron(H,H)
	u = np.kron(u,H)
	synth(u,8,3)

def test_XXX(): # should return 111
	u = np.kron(X,X)
	u = np.kron(u,X)
	synth(u,8,3)
def test_HXX(): # should return 01,11
	u = np.kron(H,X)
	u = np.kron(u,X)
	synth(u,8,3)
def test_XHH(): # should return 100,111,101,110
	u = np.kron(X,H)
	u = np.kron(u,H)
	synth(u,8,3)	
# test four_qubit_gates 
def test_UUUU():
	u = unitary(random(),random(),random(),random())
	u1=np.kron(u,u)
	u1=np.kron(u1,u)
	u1=np.kron(u1,u)
	synth(u1,16,4)

def test_HHHH(): # yyyy
	u= np.kron(H,H)
	u= np.kron(u,H)
	u= np.kron(u,H)
	synth(u,16,4) 
#def test_XXXH(): # 111y
#def test_TXHS(): 
#def test_XHXX(): # 1y11 

def test_HXXHXH():# y11y1  20 mins
	u = np.kron(H,X)
	u = np.kron(u,X)
	u = np.kron(u,H)
	u = np.kron(u,X)
	u = np.kron(u,H)
	synth(u,64,6)

# test five_qubit_gates
def test_HXHHX():# y11y1  20 mins
	u = np.kron(H,X)
	u = np.kron(u,H)
	u = np.kron(u,H)
	u = np.kron(u,X)
	synth(u,32,5)

def test_HXXHXHH():# y11y1  20 mins
	u = np.kron(H,X)
	u = np.kron(u,X)
	u = np.kron(u,H)
	u = np.kron(u,X)
	u = np.kron(u,H)
	u = np.kron(u,H)
	synth(u,128,7)

#test_H()
#test_UUU()
#test_UUUU()
#test_HHHH()
test_HXHHX()
#test_HXXHXHH() 
 
