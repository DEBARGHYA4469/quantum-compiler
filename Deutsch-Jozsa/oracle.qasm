# 
# File:   test2.qasm
# Date:   22-Mar-04
# Author: I. Chuang <ichuang@mit.edu>
#
# Sample qasm input file - simple teleportation circuit
#

 qubit q0
 qubit q1
 qubit q2
 qubit q3
 qubit q4
 qubit q5
 h q0
 h q1
 h q2
 h q3
 h q4
 X q5
 h q5
 cnot q0,q5
 cnot q1,q5
 cnot q2,q5
 cnot q4,q5
 h q0
 h q1
 h q2
 h q3
 h q4
 measure q0
 measure q1
 measure q2
 measure q3
 measure q4
