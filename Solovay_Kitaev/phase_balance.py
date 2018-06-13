#					GLOBAL PHASE BALANCER 
#	U  be any unitary matrix st UU*=I . 
#	U' be its transform in SU(2).
#	then U' = [1/det(U)]^0.5 

from util import determinant 
import numpy as np 
import cmath,math

def SU2(U):
	t  = complex(0,0) 
	t = t + determinant(U)
	globalPhase = (1/t)**0.5  
	return U/globalPhase 

#checked............................................................
