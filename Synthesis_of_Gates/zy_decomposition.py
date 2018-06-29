import math
import numpy as np 
from random import randint

def e_i(x): # complexexponent
	return (math.cos(x) + 1j*math.sin(x))  

def unitary(alpha,beta,gama,delta):               # Generate unitary matrix 
	mat = np.array([[1,0],[0,1]])
	a1 = e_i(-1.0*beta/2.0)
	a2 = e_i(beta/2.0)
	b1 = e_i(delta/2.0)
	b2 = e_i(-1.0*delta/2.0)
	c = e_i(alpha)
	A_mat = np.array([[a1,0],[0,a2]])
	B_mat = np.array([ [math.cos(gama/2.0),math.sin(gama/2.0)*-1], [math.sin(gama/2.0),math.cos(gama/2.0)] ])
	C_mat = np.array([[b2,0],[0,b1]])
	mat = np.matmul(A_mat,np.matmul(B_mat,np.matmul(C_mat,mat)))*c
	return mat

def determinant(Mat):  #determinant of a 2x2 matrix 
	return (Mat[0][0]*Mat[1][1]-Mat[0][1]*Mat[1][0])

def SU2(U): #................................................ checked.............................................. 
	t  = complex(0,0) 
	t = t + determinant(U)
	globalPhase = t**0.5  
	return U/globalPhase,globalPhase 
            # ...............................................checked .....................................
def rotateY(t):
	Y = np.array([[math.cos(t/2) ,-1.*math.sin(t/2)],[math.sin(t/2) , math.cos(t/2)]])
	return Y 	

def rotateZ(t):
	Z = np.array([[math.cos(t/2)-1j*math.sin(t/2) , 0],[0 ,math.cos(t/2) + 1j*math.sin(t/2)]])
	return Z 

def checkGC(U,a,b,g,d): # check the zy_decompostion 
	if(np.allclose((math.cos(a)+1j*math.sin(a))*np.matmul(rotateZ(b),np.matmul(rotateY(g),rotateZ(d))),U) == True): return True
	return False

def almost_equal(a,b):
	if(abs(a-b) < 0.00001):
		return True
	else : 
		return False  

def get_theta(cos,sin) : # solves inverse trigonometric equations 
	PI = math.acos(-1.0) 	
	
	# EXCEPTION CHECK .................................
	if(almost_equal(cos, 1.0)): cos =  1.0
	if(almost_equal(cos,-1.0)): cos = -1.0
	if(almost_equal(sin, 1.0)): sin =  1.0
	if(almost_equal(sin,-1.0)): sin = -1.0 
	# EXCEPTION CHECK .................................
	try:
		if(sin >= 0 and cos >= 0): # FIRST QUADRANT ! 
			return        math.asin(abs(sin)) 
		if(sin >= 0 and cos <= 0): # SECOND QUADRANT ! 
			return     PI-math.asin(abs(sin))
		if(sin <= 0 and cos <= 0): # THIRD QUADRANT ! 
			return   PI + math.asin(abs(sin))
		if(sin <= 0 and cos >= 0): # FOURTH QUADRANT ! 
			return 2*PI - math.asin(abs(sin)) 
	except:
		print(sin,cos)

def dist(U,a,b,g,d):
	V = (math.cos(a) + 1j*math.sin(a))*np.matmul(rotateZ(b),np.matmul(rotateY(g),rotateZ(d)))
	e = 0
	e = (U[0][0]-V[0][0])**2 + (U[1][0]-V[1][0])**2 + (U[0][1]-V[0][1])**2 + (U[1][1]-V[1][1])**2 
	return e


def zy_decompose(U): # mod = 0 , 1 representing positive and negative values  
	SU,gp = SU2(U) # extract global-phase 
	#if(mod>1) : return 
	
	# ........................SAVE ALPHA.........................................
	try:
		alpha = get_theta(gp.real,gp.imag)	
	except: # overflow 
		if(gp.real > 0): alpha = math.acos(1.0)
		if(gp.real < 0): alpha = math.acos(-1.0)
		if(gp.imag > 0): alpha = math.asin(1.0)
		if(gp.imag < 0): alpha = math.asin(-1.0)	
	# ..........................SAVED............................................
	
	u1 = SU[0][0]
	u2 = SU[1][0]

	#...........................SAVE GAMA........................................
	try: 
		gama = 2*get_theta(abs(u1),abs(u2))
	except:	
		if(abs(u1) > 0): gama=math.acos(1.0)
		if(abs(u2) > 0): gama=math.asin(1.0)
	# ..........................SAVED............................................
	#try:
	#	u1 = u1 / math.cos(gama/2.0)
	#	u2 = u2 / math.sin(gama/2.0)
	#except: # division by 0 
	#	if(np.allclose(math.cos(gama/2.0),0)==True): u1 = 0 
	#	if(np.allclose(math.sin(gama/2.0),0)==True): u2 = 0  
	
	if(almost_equal(math.cos(gama/2.0),0)): u1 = 0 
	if(almost_equal(math.sin(gama/2.0),0)): u2 = 0
	if(almost_equal(math.cos(gama/2.0),0) == False and almost_equal(math.sin(gama/2.0),0) == False):
		u1 = u1 / math.cos(gama/2.0)
		u2 = u2 / math.sin(gama/2.0)

	theta1 = get_theta(u1.real,-1.0*u1.imag)
	theta2 = get_theta(u2.real,u2.imag)		
	#print(theta1,theta2,u1,u2)
	# ..........................SAVE ALPHA,BETA..................................
	beta  = theta1 + theta2  
	delta = theta1 - theta2	
	#print(U) 
	# ...............................SAVED.......................................
	#if(checkGC(U,alpha,beta,gama,delta) == False):
	#	return "Failure!!",dist(U,alpha,beta,gama,delta)
	#if(checkGC(U,alpha,beta,gama,delta) == True ):
		#print("SUCCESS!")
	return alpha,beta,gama,delta
#.................................................................................
#................................................................................
def eAXBXC(U):
	a,b,g,d = zy_decompose(U)
	A = np.matmul(rotateZ(b),rotateY(g/2.0))
	B = np.matmul(rotateY(-1.0*g/2),rotateZ(-(d+b)/2.0))
	C = rotateZ((d-b)/2.0)
	X = np.array([[0,1],[1,0]])
	e = math.cos(a) + 1j*math.sin(a)
	return a,A,B,C,X 
# ................................................................................







