# ....................................................MAIN SOLOVAY KITAEV ROUTINE ....................................

import numpy as np
from approx import multiply 
from NNSearch import *
# from GC_Decompose import *
from util import dagger
from Similar import *

def dagger_seq(sequence): # HTHTSST ------> tssthth
	ll  = len(sequence)
	new_sequence = "" 
	for i in range(ll):
		if(sequence[i]=="H"): new_sequence = 'h' + new_sequence 
		if(sequence[i]=="h"): new_sequence = 'H' + new_sequence
		if(sequence[i]=="T"): new_sequence = 't' + new_sequence
		if(sequence[i]=="t"): new_sequence = 'T' + new_sequence
		if(sequence[i]=="I"): new_sequence = 'I' + new_sequence
		if(sequence[i]=="s"): new_sequence = 'S' + new_sequence
		if(sequence[i]=="S"): new_sequence = 's' + new_sequence
	return new_sequence 
e0 = 0.1

def solovay_kitaev(U,n): # input :: gate , depth 
	if(n==0): # base case 
		return basic_approx(U,e0) # returns the sequence 
	else :
		U_prev_seq = solovay_kitaev(U,n-1) # U(n-1)
		U_prev = multiply(U_prev_seq)
		V,W = decompose(np.matmul(U,dagger(U_prev)))
		V_prev_seq = solovay_kitaev(V,n-1)
		W_prev_seq = solovay_kitaev(W,n-1)
		U_now_seq  = V_prev_seq + W_prev_seq + dagger_seq(V_prev_seq) + dagger_seq(W_prev_seq) + U_prev_seq
		return U_now_seq 




