from util import *
from approx import *
from math import *

def commutator(V,W): # multiplicative group commutator
	return np.matmul(V,np.matmul(W,np.matmul(dagger(V),dagger(W))))


def Unit(a,b,c,d):
	return np.array([[complex(a,-1.*d),complex(-1.*c, -1.*b) ],[complex(c,-1.*b) ,complex(a,d)]])

def tolequal(r1,r2): # equality within a tolerance  
	if(abs(r1-r2) < 0.000001):
		return True
	else : 
		return False  
# UTILITY ........................................................................................................................

def sendSimilar(U) : # main funtion for group commutator decompose 
	I = U[0][0].real 
	X = (-1)*(U[1][0].imag)
	Y = U[1][0].real
	Z = (-1)*(U[0][0].imag)

	s   =  pow((1.- I)/2,0.25)
	c   =  sqrt(1.- pow(s,2))
	
	nn  =  sqrt(1 - pow(I,2))
	nx  =  X / nn 
	ny  =  Y / nn 
	nz  =  Z / nn 
	mn  =  sqrt(1. - pow(I,2))
	mx  =  2*pow(s,3)*c/mn
	my  =  -2*pow(s,3)*c/mn
	mz  =  -2*pow(s,2)*pow(c,2)/mn

	x = (nx + mx)/2 
	y = (ny + my)/2
	z = (nz + mz)/2
 
	n = sqrt(x*x+y*y+z*z)
	x=x/n
	y=y/n
	z=z/n
	V = Unit(c,s,0,0)
	W = Unit(c,0,s,0)
	S = Unit(0,x,y,z)
	Vt = np.matmul(S,np.matmul(V,dagger(S)))
	Wt = np.matmul(S,np.matmul(W,dagger(S)))
	
	
	w1,v1  = diagonalize(U)
	w2,v2  = diagonalize(commutator(Vt,Wt))
	return Vt,Wt
	
def decompose(U):
	Vt,Wt = sendSimilar(U)	
	V = commutator(Vt,Wt)
	w1,v1  = diagonalize(U)
	w2,v2  = diagonalize(V)
	if( abs(w1[0].imag-w2[0].imag) < 0.000001 ): # same permute
		# Eigen values same ...........................................................
		K = SU2(np.matmul(v2,np.matmul( np.array([[w2[0],0],[0,w2[1]]]) ,dagger(v2))))
		S=np.matmul(v1,dagger(v2))		
		Unew = np.matmul(S,np.matmul(K,dagger(S)))
		Vt = np.matmul(S,np.matmul(Vt,dagger(S)))
		Wt = np.matmul(S,np.matmul(Wt,dagger(S)))	
		return Vt,Wt
	else:
		# Eigen values swapped ....................................................
		v2 = np.array([[v2[0][1],v2[0][0]],[v2[1][1],v2[1][0]]])
		S=np.matmul(v1,dagger(v2))		
		Vt = np.matmul(S,np.matmul(Vt,dagger(S)))
		Wt = np.matmul(S,np.matmul(Wt,dagger(S)))
		Unew = np.matmul(S,np.matmul(V,dagger(S)))		
		return Vt,Wt
























