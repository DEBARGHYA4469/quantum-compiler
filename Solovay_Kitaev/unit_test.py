# This contains some of the unit tests 

from util import * 
from approx import *
from recur_sk import *


#  GROUP_COMMUTATOR_DECOMPOSE_TEST 
def test_gc_decompose():
	U = generate_SU2() 
	V,W = decompose(U)
	if(isEqual(U,commutator(V,W)) == True): print("Success")	
	else : print("Failure")

test_gc_decompose()
