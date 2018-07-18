# This module contains all the utility codes 

import numpy as np 
from random import randint
import math,cmath
from numpy import linalg as LA

def dagger(Mat):       # adjoint of the matrix 
	return Mat.transpose().conjugate()

def determinant(Mat):  #determinant of a 2x2 matrix 
	return (Mat[0][0]*Mat[1][1]-Mat[0][1]*Mat[1][0])

def norm(vector):      #euclidean norm of a vector in R^n
	return  math.sqrt((vector[0]**2 + vector[1]**2 + vector[2]**2))


def isEqual(Mat1,Mat2):  
	tol = 0.00001 # check if two matrix are equal within the tolerance 
	r1 = abs(Mat1[0][0] - Mat2[0][0])
	r2 = abs(Mat1[1][0] - Mat2[1][0])
	r3 = abs(Mat1[0][1] - Mat2[0][1])
	r4 = abs(Mat1[1][1] - Mat2[1][1])
	if(r1 < tol and r2 < tol and r3 < tol and r4 < tol ) : return True
	return False,r1,r2,r3,r4
	
def is_SU2(Mat): #check if its in SU(2)
	I = np.array([[1,0],[0,1]])	
	tolerance = 0.01	
	if(abs(determinant(Mat)-1.0) < tolerance and isEqual(np.matmul(Mat,dagger(Mat)),I) == True): return True
	else : return False

def proj_Pauli(tx,ty,tz):	
	return np.array([[tz,complex(tx,-1.0*ty)],[complex(tx,ty),-1.0*tz]]) 

def generate_SU2():   # generate a random matrix in SU2() 
	theta = [randint(1,10) for i in range(3)]
	t = norm(theta)
	theta[0] = theta[0] / t 
	theta[1] = theta[1] / t
	theta[2] = theta[2] / t
	I = np.array([[1,0],[0,1]])
	Su2= math.cos(t/2)*I - 1j*math.sin(t/2)*proj_Pauli(theta[0],theta[1],theta[2])
	return Su2

def diagonalize(U):            # Diagonalize a matrix 
	w1,v1 = LA.eig(U)
	return w1,v1
	


# checked ...........................................................................................
