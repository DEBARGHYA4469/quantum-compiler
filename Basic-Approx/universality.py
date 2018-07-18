# Universality construction .... 

import numpy as np 
from math import *
from phase_balance import *
from Gates import *
from random import random

PI = acos(-1.0)
I = np.array([[1,0],[0,1]]) 
X = np.array([[0,1],[1,0]])
Y = np.array([[0,-1j],[1j,0]]) 
Z = np.array([[1,0],[0,-1]])

def norm(vec):
	return sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)

def rotate_n1(theta): # THTH
	n1 = np.array([cos(PI/8.),sin(PI/8.),cos(PI/8.)])
	n1_cap = n1/norm(n1)
	mat = cos(theta/2.)*I - 1j*sin(theta/2.)*(n1_cap[0]*X + n1_cap[1]*Y + n1_cap[2]*Z) 
	return mat

def rotate_n2(theta): # HTHT
	n2 = np.array([cos(PI/8.),-1.*sin(PI/8.),cos(PI/8.)])
	n2_cap = n2/norm(n2)
	mat = cos(theta/2.)*I - 1j*sin(theta/2.)*(n2_cap[0]*X + n2_cap[1]*Y + n2_cap[2]*Z) 
	return mat

# TESTING .............................

def sqnc(q1,q2,q3): # ...THTH....HTHT....THTH
	r1 = "THTH"
	r2 = "HTHT"
	ret = ""
	for __ in range(int(q1)):ret=ret + r1 
	for __ in range(int(q2)):ret=ret + r2 
	for __ in range(int(q3)):ret=ret + r1
	return ret 

from approx import *

file = open("gates.dat","a")

y = 0
average = 0
e=0.125
while(y < 30):
	t1 = random()
	t2 = random()
	t3 = random()
	mat = np.matmul(rotate_n1(t1),np.matmul(rotate_n2(t2),rotate_n1(t3)))
	print(mat)
	theta = 2*acos(cos(PI/8.0)**2)
	#epsilon = float(input()) # float ! 
	q1,q2,q3 = 0,0,0
	mod = 2*PI
	while(1):
		q1 = q1 + 1 
		#if(abs( (q1*theta)%mod - t1 ) < e): break;
		if(trace_norm(rotate_n1(t1),rotate_n1(q1*theta)) < e/3):break		
	while(1):
		q2 = q2 + 1 
		#if(abs( (q2*theta)%mod - t2 ) < e): break;
		if(trace_norm(rotate_n1(t2),rotate_n1(q2*theta)) < e/3):break	
	while(1):
		q3 = q3 + 1 
		#if(abs( (q3*theta)%mod - t3 ) < e): break;
		if(trace_norm(rotate_n1(t3),rotate_n1(q3*theta)) < e/3):break	
	#print(q1,q2,q3)	

	mat2 = np.matmul(rotate_n1(q1*theta),np.matmul(rotate_n2(q2*theta),rotate_n1(q3*theta)))
		
	print(trace_norm(mat,mat2))
	#print("No.of gates used : ",q1+q2+q3)
	print("Sequence:",sqnc(q1,q2,q3))
	#print(e,trace_norm(mat,mat2))
	# file.write(str(q1+q2+q3))\
	input()
	average = average + q1 +q2 + q3 
	print("---------------------------------------- ")
	y=y+1
average = average / 30
file.write(str(e)+" "+str(average)+"\n")
# TESTING .............................


