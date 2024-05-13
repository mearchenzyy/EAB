import sys
# sys.path.append(r"Y:\Users\Yingyue\Gates_Lab_Suite-master")
sys.path.append(r"/Volumes/funkflower/Users/Yingyue/Gates_Lab_Suite-master")
from Core_Definition import *
from Auto_Algorithm import *
from Visualization import *
import os
from SPAM import *
import numpy as np
# import xlsxwriter as xlsx
from scipy import optimize
import random
import math
from math import floor
import uuid

def state_init_bell_pairs(qc,n):
    for i in range (n):
        qc.Add_Gate(Quantum_Gate("HAD",i))
    for i in range (n):
        qc.Add_Gate(Quantum_Gate("CNOT",i,int(i+n)))
        
def state_init_bell_pairs_01(qc,n):
    for i in range (n):
        qc.Add_Gate(Quantum_Gate("HAD",i))
        qc.Add_Gate(Quantum_Gate("RX",i+n,angle=1))
    for i in range (n):
        qc.Add_Gate(Quantum_Gate("CNOT",i,int(i+n)))

def state_init_bell_pairs_explicitCNOT(qc,n):
    for i in range (n):
        qc.Add_Gate(Quantum_Gate("HAD",i))
    for i in range (n):
#         qc.Add_Gate(Quantum_Gate("CNOT",i,int(i+n)))
        qc.Add_Gate(Quantum_Gate("RY",i,angle=2))
        qc.Add_Gate(Quantum_Gate("FTXA",i,int(i+n),angle=1/4))
        qc.Add_Gate(Quantum_Gate("RY",i,angle=-2))
        qc.Add_Gate(Quantum_Gate("RX",int(i+n),angle=-2))
        qc.Add_Gate(Quantum_Gate("AZ",i,angle=(-1/2)))

def state_init_bell_pairs_explicitCNOT_01(qc,n):
    for i in range (n):
        qc.Add_Gate(Quantum_Gate("HAD",i))
        qc.Add_Gate(Quantum_Gate("RX",int(i+n),angle=1))
    for i in range (n):
#         qc.Add_Gate(Quantum_Gate("CNOT",i,int(i+n)))
        qc.Add_Gate(Quantum_Gate("RY",i,angle=2))
        qc.Add_Gate(Quantum_Gate("FTXA",i,int(i+n),angle=1/4))
        qc.Add_Gate(Quantum_Gate("RY",i,angle=-2))
        qc.Add_Gate(Quantum_Gate("RX",int(i+n),angle=-2))
        qc.Add_Gate(Quantum_Gate("AZ",i,angle=(-1/2)))
        
#should we set_mapping in this function? need to think about if it takes mapping into account   
def add_pauli_twirl(qc,n):
    pauliLayer = [random.choice(['I','X','Y','Z']) for j in range(n)]
    q_index=0
    for pauli in pauliLayer :
        if (pauli=="I"):
            pass
        elif (pauli=="X"):
            qc.Add_Gate(Quantum_Gate("RX",q_index,angle=1))
        elif (pauli=="Y"):
            qc.Add_Gate(Quantum_Gate("RY",q_index,angle=1))
        elif (pauli=="Z"):
            qc.Add_Gate(Quantum_Gate("AZ",q_index,angle=1))
        q_index+=1
    return pauliLayer
    
def add_clifford_layer(qc,n,clifford):
    if (clifford=="CNOT"):
        for i in range (floor(n/2)):
            qc.Add_Gate(Quantum_Gate("CNOT",2*i,2*i+1))
    if (clifford=="XX"):
        for i in range (floor(n/2)):
            qc.Add_Gate(Quantum_Gate("FTXA",2*i,2*i+1,angle=1/4))
            
def bell_measurement(qc,n):
    for i in range (n):
        qc.Add_Gate(Quantum_Gate("CNOT",i,int(i+n)))
    for i in range (n):
        qc.Add_Gate(Quantum_Gate("HAD",i))
        
def bell_measurement_explicitCNOT(qc,n):
    for i in range (n):
#         qc.Add_Gate(Quantum_Gate("CNOT",i,int(i+n)))
        qc.Add_Gate(Quantum_Gate("RY",i,angle=2))
        qc.Add_Gate(Quantum_Gate("FTXA",i,int(i+n),angle=1/4))
        qc.Add_Gate(Quantum_Gate("RY",i,angle=-2))
        qc.Add_Gate(Quantum_Gate("RX",int(i+n),angle=-2))
        qc.Add_Gate(Quantum_Gate("AZ",i,angle=(-1/2)))
    for i in range (n):
        qc.Add_Gate(Quantum_Gate("HAD",i))

def map_statepop_2_ibm_mapping(counts,n):
    counts_ibm_mapping=[0 for i in range (2**(2*n))]
    idx_ibm_mapping=[]
    counts_ibm_mapping_dic={}
    if len(str((2*n)))<2:
        f="00"+str(2*n)+"b"
    elif len(str((2*n)))<3:
        f="0"+str(2*n)+"b"
    else:
        raise ValueError("n is too big")
    for gates_idx in range (2**(2*n)):
        gates_idx_str=format(gates_idx,f)
        ibm_idx=0
        for i in range (2*n):
            ibm_idx+=2**(i)*int(gates_idx_str[i])
        ibm_idx_str=format(ibm_idx,f)
#         idx_ibm_mapping.append(ibm_idx_str)
        counts_ibm_mapping[ibm_idx]=counts[gates_idx]
        counts_ibm_mapping_dic[ibm_idx_str]=counts[gates_idx]
    return counts_ibm_mapping, counts_ibm_mapping_dic
        
        
#for qiskit
def prepare_bell_state_1q(circuit,index1,index2):
	circuit.h([index1])
	circuit.cx([index1],[index2])

def prepare_bell_state_1q_01(circuit,index1,index2):
	circuit.h([index1])
	circuit.x([index2])
	circuit.cx([index1],[index2])
        
def bell_measurement_1q(circuit,index1,index2): 
	# info qubit at index1
	circuit.cx([index1],[index2])
	circuit.h([index1])

def pauli_gate_1q(circuit,index,pauli=None): #For stabilizer simulator to work, cannot use Pauli class
	if pauli == 'I':
		circuit.id([index])
	elif pauli == 'Z':
		circuit.z([index])
	elif pauli == 'X':
		circuit.x([index])
	elif pauli == 'Y':
		circuit.y([index])
	else:
		assert 1==0

def add_XX_clifford(circ,index):
    circ.z(2*index)
    circ.h(2*index)
    circ.x(2*index)
    circ.z(2*index)
    circ.cx(2*index,2*index+1)
    circ.s(2*index)
    circ.h(2*index)
    circ.x(2*index)
    circ.h(2*index+1)
    circ.s(2*index+1)
    circ.h(2*index+1)