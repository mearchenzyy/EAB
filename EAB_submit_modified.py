import numpy as np
import sys, json, pickle, random
import qiskit
from qiskit import IBMQ, QuantumCircuit, execute
from qiskit.compiler.transpiler import transpile
from qiskit.providers.ibmq.managed import IBMQJobManager
from qiskit.extensions import UnitaryGate
from qiskit.quantum_info import Pauli, Clifford
from scipy.stats import sem, unitary_group
from scipy.linalg import sqrtm,expm


def prepare_bell_state_1q(circuit,index1,index2):
	circuit.h([index1])
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

def measure_pauli_1q(circuit,index,pauli=None):
	if pauli == 'I' or pauli == 'Z':
		circuit.id([index])
	elif pauli == 'X':
		circuit.h([index])
	elif pauli == 'Y':
		circuit.s([index])
		circuit.s([index])
		circuit.s([index])
		circuit.h([index])
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


def submit_eab(n,n_total,Lrange,C,batch,qubit_map,clifford_layer='Id',gset="Pauli",repeat=None,periodic=False,use_density_matrix=False):
	data_save = {}
	q = qubit_map
	eab_circ_all = []
	for b in range(batch):
		# data_batch = {}
		# num_jobs = 0

		circuits_batch = []
		data_batch = []

		for l in range(len(Lrange)):
			L = Lrange[l]
			for c in range(C):
				# run the circuit
				job_save = {}

				# n_total = # of qubits on devices. [1..n] qubits actually used
				state = QuantumCircuit(2*n,2*n)
				gates = QuantumCircuit(n)

				# state preparation
				# Bell state Psi_0
				for j in range(n):
					# may need to take device structure into consideration
					prepare_bell_state_1q(state,j,j+n)
					# prepare_bell_state_1q(state,q[j],q[j+n])
					# prepare_pauli_eigenstate_1q(state,q[j],pauli=pauliList[j])
					# apply_1q(circuit,q[j],gate = np.expm(pauli_sample[j]))
				state.barrier()

				for i in range(L):
					pauliLayer = [random.choice(['I','X','Y','Z']) for j in range(n)]
					#pauliTrans = Pauli(''.join(pauliLayer[::-1]))
					for j in range(n):
						# gates.pauli(pauliLayer[j],[q[j]])
						pauli_gate_1q(gates,j,pauli=pauliLayer[j])
						#pauli_gate_1q(gates,q[j],pauli=pauliLayer[j])
					
					'''Clifford Layer'''
					if clifford_layer == 'Id':
						pass
					elif clifford_layer == 'CNOT':
						ngates = int(n/2)
						for j in range(ngates):
							gates.cx(2*j,2*j+1)
					elif clifford_layer == 'XX':
						ngates = int(n/2)
						for j in range(ngates):
							add_XX_clifford(gates,j)
# 							gates.cz(2*j,2*j+1)
# 							gates.cx(2*j,2*j+1)
                           
					gates.barrier()

					# ngates = int(n/2)
					# for j in range(ngates):
					# 	gates.cx(q[2*j],q[2*j+1])
					# if n%2 == 1:
					# 	gates.id(q[n-1])
					#pauliOp = pauliOp.evolve(pauliTrans).evolve()		

					# start from qubit 0
					# ngates = int(n/2)
					# for j in range(ngates):
					# 	apply_1q_random(circuit,q[2*j],gset=gset,record_gates=gates)
					# 	apply_1q_random(circuit,q[2*j+1],gset=gset,record_gates=gates)
					# 	circuit.cx(q[2*j],q[2*j+1])
					# if n%2 == 1:
					# 	apply_1q_random(circuit,q[n-1],gset=gset,record_gates=gates)
					# 	circuit.id(q[n-1])

				# final layer:
				pauliLayer = [random.choice(['I','X','Y','Z']) for j in range(n)]
				# pauliTrans = Pauli(''.join(pauliLayer[::-1]))
				for j in range(n):
					# gates.pauli(pauliLayer[j],[q[j]])
					pauli_gate_1q(gates,j,pauli=pauliLayer[j])
					# pauli_gate_1q(gates,q[j],pauli=pauliLayer[j])


				# # calculate C(P), which decides our measurement setting
				# pauliOp = Pauli(''.join(pauliList[::-1])) # join in reverse order
				# pauliOp = pauliOp.evolve(Clifford(gates).adjoint()) # note: adjoint is necessary for Heisenberg evolution.

				'''The choice of depth guarantees 'gates' is Pauli'''	
				cliffordOp = Clifford(gates)



				circuit = state.compose(gates,range(n))

				#transpile circuit only here?
				

				if use_density_matrix:
					circuit.save_density_matrix()
				else:
					circuit.barrier()

					# # measurement for fidelity estimation
					# measurement_setting = pauliOp.to_label()
					# # there could be a '-' in the Pauli label
					# if measurement_setting[0].isupper() is False:
					# 	measurement_setting = measurement_setting[1:]
					# measurement_setting = measurement_setting[::-1]
					# #print(measurement_setting)
					# for j in range(n):
					# 	measure_pauli_1q(circuit,q[j],pauli=measurement_setting[j])

					for j in range(n):
						bell_measurement_1q(circuit,j,j+n)
						# bell_measurement_1q(circuit,q[j],q[j+n])

					circuit.barrier()
					#circuit.measure([q[i] for i in range(2*n)], [i for i in range(2*n)])
					circuit.measure([i for i in range(2*n)], [i for i in range(2*n)])
				
				#circuit.draw(filename="EAB_circuit")

				### use one of the following lines:
				# circuit = circuit.decompose().decompose()
				# circuit = qiskit.transpile(circuit,optimization_level=1,basis_gates=basis_gates)
				
				# circuit = transpile(circuit,initial_layout=q)

				# circuit.draw(output='mpl',filename='circuit3.png')
				# sys.exit(0)

				R = 1
				if repeat is not None:
					R = repeat[l]
				for r in range(R):
					circuits_batch.append(circuit)

				job_save["n"] = n
				job_save["L"] = L
				job_save["c"] = c
				# job_save["type"] = "cross_entropy_H"
				# job_save["gates"] = gates
				# job_save["pauli"] = pauliOp
				job_save['clifford'] = cliffordOp.to_dict() # actually Pauli
				job_save["repeat"] = R
				job_save["clifford_layer"] = clifford_layer

				# job_save["job_id"] = job_id

				data_batch.append(job_save)
				# num_jobs += 1
				# print(num_jobs)

		eab_circ_all.append(circuits_batch)
		data_save["batch_%d" % b] = data_batch
	return data_save, eab_circ_all

def submit_eab_ancilla_pauli_twirl(n,n_total,Lrange,C,batch,qubit_map,clifford_layer='Id',gset="Pauli",repeat=None,periodic=False,use_density_matrix=False):
	data_save = {}
	q = qubit_map
	eab_circ_all = []
	for b in range(batch):

		circuits_batch = []
		data_batch = []

		for l in range(len(Lrange)):
			L = Lrange[l]
			for c in range(C):
				# run the circuit
				job_save = {}

				# n_total = # of qubits on devices. [1..n] qubits actually used
				state = QuantumCircuit(2*n,2*n)
				gates_all= QuantumCircuit(2*n,2*n)
				gates = QuantumCircuit(n)
				gates_anc = QuantumCircuit(n)
				

				# state preparation
				# Bell state Psi_0
				for j in range(n):
					# may need to take device structure into consideration
					prepare_bell_state_1q(state,j,j+n)
				state.barrier()

				for i in range(L):
					pauliLayer = [random.choice(['I','X','Y','Z']) for j in range(n)]
					#pauliTrans = Pauli(''.join(pauliLayer[::-1]))
					for j in range(n):
						pauli_gate_1q(gates_all,j,pauli=pauliLayer[j])
						pauli_gate_1q(gates,j,pauli=pauliLayer[j])

					pauliLayer = [random.choice(['I','X','Y','Z']) for j in range(n)]
					for j in range(n,2*n):
						pauli_gate_1q(gates_all,j,pauli=pauliLayer[j])
						pauli_gate_1q(gates_anc,(j-n),pauli=pauliLayer[j])
					
					'''Clifford Layer'''
					if clifford_layer == 'Id':
						pass
					elif clifford_layer == 'CNOT':
						ngates = int(n/2)
						for j in range(ngates):
							gates.cx(2*j,2*j+1)
					elif clifford_layer == 'XX':
						ngates = int(n/2)
						for j in range(ngates):
							add_XX_clifford(gates,j)
                           
					gates.barrier()

				# final layer:
				pauliLayer = [random.choice(['I','X','Y','Z']) for j in range(n)]
				for j in range(n):
					pauli_gate_1q(gates_all,j,pauli=pauliLayer[j])
					pauli_gate_1q(gates,j,pauli=pauliLayer[j])
				
				gates_anc_inv=gates_anc.inverse()
				gates_all = gates_anc_inv.compose(gates_all,range(n,2*n))


				'''The choice of depth guarantees 'gates' is Pauli'''	
				cliffordOp = Clifford(gates)



				circuit = state.compose(gates_all,range(n))

				#transpile circuit only here?
				

				if use_density_matrix:
					circuit.save_density_matrix()
				else:
					circuit.barrier()

					for j in range(n):
						bell_measurement_1q(circuit,j,j+n)
						# bell_measurement_1q(circuit,q[j],q[j+n])

					circuit.barrier()
					#circuit.measure([q[i] for i in range(2*n)], [i for i in range(2*n)])
					circuit.measure([i for i in range(2*n)], [i for i in range(2*n)])

				R = 1
				if repeat is not None:
					R = repeat[l]
				for r in range(R):
					circuits_batch.append(circuit)

				job_save["n"] = n
				job_save["L"] = L
				job_save["c"] = c
				job_save['clifford'] = cliffordOp.to_dict() # actually Pauli
				job_save["repeat"] = R
				job_save["clifford_layer"] = clifford_layer


				data_batch.append(job_save)

		eab_circ_all.append(circuits_batch)
		data_save["batch_%d" % b] = data_batch
	return data_save, eab_circ_all
