import numpy as np
import sys, json, copy, pickle
import matplotlib.pyplot as plt
import qiskit
import random
from qiskit.result import Result
from qiskit import IBMQ, QuantumCircuit, execute
from qiskit.providers.aer import QasmSimulator, StatevectorSimulator
from qiskit.extensions import UnitaryGate
from qiskit.quantum_info.operators.symplectic import pauli
from scipy.stats import sem, entropy, linregress
from scipy.optimize import curve_fit
from qiskit.quantum_info import Pauli, Clifford
import pdb
# IBMQ.save_account('b3460dbc07ed93247ba3dd87b6619d71872d5d079f3f01bd5944678aa544b97203807ffcff040ca6d440ad990d907bbe59489179c190bd7b6670bf432e874940')

# IBMQ.load_account()
# provider = IBMQ.get_provider(hub='ibm-q-internal', group='deployed', project='default')
# backend = provider.get_backend('ibmq_athens')
# backend = StatevectorSimulator()

# Lrange = list(range(1,16))

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

def process_EAB(Lrange, cb_data, pauli_request_list, counts_batch=None, repeat=None, periodic=False, use_density_matrix = False, inject_pauli_sparse = None):
    
    fidelity_list = {}
    for pauli in pauli_request_list:
        fidelity_list[pauli] = {}
        for L in Lrange:
            fidelity_list[pauli][L] = []
    batch=1
    for b in range(batch):
#         data_batch = cb_data["batch_%d" % b]
#         result_batch = Result.from_dict(cb_data["result"][b])
        data_batch = cb_data
#         pdb.set_trace()

        circuit_count = 0  ### To extract data from result_batch

        for i in range(len(data_batch)):
            # print("batch", b, "circuit", i)
            job_data = data_batch[i]
            n = job_data["n"]
            L = job_data["L"]
            clifford_layer = job_data["clifford_layer"]
            #print("n=%d" % n, "batch", b, "circuit", i)
            # print(job_data["circuit"].item().shape)
            # sys.exit(0)

            

            ###
            # circuit = construct_circuit(job_data["n"], job_data["L"], job_data["circuit"], periodic=periodic)
            # job = execute(circuit, backend, backend_options={"max_parallel_threads": 0})
            # job = execute(circuit, backend, max_parallel_threads=0)
            # state_vector = np.array(job.result().get_statevector())
            
            # pauliOp = job_data["pauli"]
            cliffordOp = Clifford.from_dict(job_data['clifford'])

            '''TODO: add inject Pauli noise according to inject_pauli_sparse
            This should be done by calculated a new CliffordOp
            How efficient can we do that?
            '''
            if inject_pauli_sparse is not None: # let's inject some noise
                inject_circ = QuantumCircuit(n)
                for _k in range(L):
                    # clifford layer
                    if clifford_layer == "Id":
                        pass
                    elif clifford_layer == 'CNOT':
                        ngates = int(n/2)
                        for j in range(ngates):
                            inject_circ.cx(2*j,2*j+1)
                    else:
                        raise ValueError('Clifford layer undefined or not supported.')
                    # pauli layer
                    inject_pauli_sample = np.random.choice(inject_pauli_sparse[0],p=inject_pauli_sparse[1])
                    for j in range(n):
                        pauli_gate_1q(inject_circ,j,inject_pauli_sample[n-j-1]) # hope the ordering is correct

                injectOp = Clifford.from_circuit(inject_circ)
                cliffordOp = cliffordOp.compose(injectOp)

            if use_density_matrix is True:
                # rho = result_batch.data(circuit_count)['density_matrix']
                # F = np.real(np.trace(rho @ pauliOp.to_matrix()))
                # # print(F)
                pass
            else:
                # assert np.mod(pauliOp.phase,2) == 0
                # phase = (-1)**(pauliOp.phase>>1)
                # if phase == 1:
                #     label = pauliOp.to_label()
                # else:
                #     label = pauliOp.to_label()[1:]
                # assert len(label) == n

                # phase = pauliOp.phase
                # label = pauliOp.to_label()                
                # label = label[len(label)-n:] # remove phase

                #outcomes = result_batch.data(circuit_count)['counts']
#                 outcomes = result_batch.get_counts(circuit_count)
                outcomes = job_data["counts"]
                
                #print(outcomes)


                for pauli_label in pauli_request_list:

                    

                    pauli = Pauli(pauli_label)
                    # print("--------------")
                    # print("Request = ", pauli)
                    # print("Correction = ", pauliOp)
                    


                    assert pauli.to_label()[-n:] == pauli.evolve(cliffordOp.adjoint()).to_label()[-n:] # total clifford -> identity (except for intc_cb)
                    correction_phase = (-1)**(pauli.phase != pauli.evolve(cliffordOp.adjoint()).phase)


                    # print("correction_phase = %d" % correction_phase)
                    # label = pauli.to_label()
                    # label = label[len(label)-n:]

                    F = 0
                    tot = 0 
                    for key, counts in outcomes.items():
                        F_key = 1
                        for j in range(n):
                            # print(key)
                            # print([key[j]],[key[j+n]])
#                             pdb.set_trace()
                            pauli_bell = Pauli(([key[j+n]=='1'],[key[j]=='1'])) ### change order if not correct
                            # bell_label = Pauli(([key[j+n]],[key[j]])) ### change order if not correct
                            # print("%dth" % j)
                            # print("Bell state Pauli = ", pauli_bell)
                            # print("Request (1q) = ", Pauli(pauli_label[j]))

                            phase_1q = (-1)**(Pauli(pauli_label[j]).anticommutes(pauli_bell))
                            # print("phase_1q = %d" % phase_1q)
                            F_key *= phase_1q
                        F += F_key * counts
                        tot += counts
#                     if (i==20):
#                         pdb.set_trace()
                    
                    F = F*correction_phase/tot

                    fidelity_list[pauli_label][L].append(F)


                    # F = 0
                    # tot = 0
                    # for key, counts in outcomes.items():
                    #     F_key = 1
                    #     for j in range(n):
                    #         if label[j]!='I' and key[j] == '1':
                    #             F_key *= -1
                    #     F += F_key * counts
                    #     tot += counts

                    # F = F*phase/tot


                

            # counts = {}
            if repeat is None:
                R = 1
            else:
                R = job_data["repeat"]
#             print("circuit count:", circuit_count)

            circuit_count += R



    EAB_result = {
        "fidelity_list":                     fidelity_list,
    }
    return EAB_result


# use antiferromagnetic bell state 01+10, start from the 01 state
def process_EAB_01(Lrange, cb_data, pauli_request_list, counts_batch=None, repeat=None, periodic=False, use_density_matrix = False, inject_pauli_sparse = None):
    
    fidelity_list = {}
    for pauli in pauli_request_list:
        fidelity_list[pauli] = {}
        for L in Lrange:
            fidelity_list[pauli][L] = []
    batch=1
    for b in range(batch):
#         data_batch = cb_data["batch_%d" % b]
#         result_batch = Result.from_dict(cb_data["result"][b])
        data_batch = cb_data
#         pdb.set_trace()

        circuit_count = 0  ### To extract data from result_batch

        for i in range(len(data_batch)):
            # print("batch", b, "circuit", i)
            job_data = data_batch[i]
            n = job_data["n"]
            L = job_data["L"]
            clifford_layer = job_data["clifford_layer"]
            cliffordOp = Clifford.from_dict(job_data['clifford'])

            '''TODO: add inject Pauli noise according to inject_pauli_sparse
            This should be done by calculated a new CliffordOp
            How efficient can we do that?
            '''
            if inject_pauli_sparse is not None: # let's inject some noise
                inject_circ = QuantumCircuit(n)
                for _k in range(L):
                    # clifford layer
                    if clifford_layer == "Id":
                        pass
                    elif clifford_layer == 'CNOT':
                        ngates = int(n/2)
                        for j in range(ngates):
                            inject_circ.cx(2*j,2*j+1)
                    else:
                        raise ValueError('Clifford layer undefined or not supported.')
                    # pauli layer
                    inject_pauli_sample = np.random.choice(inject_pauli_sparse[0],p=inject_pauli_sparse[1])
                    for j in range(n):
                        pauli_gate_1q(inject_circ,j,inject_pauli_sample[n-j-1]) # hope the ordering is correct

                injectOp = Clifford.from_circuit(inject_circ)
                cliffordOp = cliffordOp.compose(injectOp)

            if use_density_matrix is True:
                pass
            else:
                if (any([L==d for d in Lrange])):
                    outcomes = job_data["counts"]
                    for pauli_label in pauli_request_list:
                        pauli = Pauli(pauli_label)
                        assert pauli.to_label()[-n:] == pauli.evolve(cliffordOp.adjoint()).to_label()[-n:] # total clifford -> identity (except for intc_cb)
                        correction_phase = (-1)**(pauli.phase != pauli.evolve(cliffordOp.adjoint()).phase)
    
                        F = 0
                        tot = 0 
                        for key, counts in outcomes.items():
                            F_key = 1
                            for j in range(n):
                                # print(key)
                                # print([key[j]],[key[j+n]])
    #                             pdb.set_trace()
                                pauli_bell = Pauli(([key[j+n]=='1'],[key[j]=='0'])) ### change order if not correct
    
                                phase_1q = (-1)**(Pauli(pauli_label[j]).anticommutes(pauli_bell))
                                # print("phase_1q = %d" % phase_1q)
                                F_key *= phase_1q
                            F += F_key * counts
                            tot += counts
    #                     if (i==20):
    #                         pdb.set_trace()
                        
                        F = F*correction_phase/tot
    #                     F = abs(F)/tot
    
                        fidelity_list[pauli_label][L].append(F)
                else:
                    continue

            if repeat is None:
                R = 1
            else:
                R = job_data["repeat"]

            circuit_count += R



    EAB_result = {
        "fidelity_list":                     fidelity_list,
    }
    return EAB_result


# use antiferromagnetic bell state 01+10, start from the 01 state
def process_EAB_01_identitytest(Lrange, cb_data, pauli_request_list, counts_batch=None, repeat=None, periodic=False, use_density_matrix = False, inject_pauli_sparse = None):
    
    fidelity_list = {}
    for pauli in pauli_request_list:
        fidelity_list[pauli] = {}
        for L in Lrange:
            fidelity_list[pauli][L] = []
    batch=1
    for b in range(batch):
        data_batch = cb_data
        circuit_count = 0  ### To extract data from result_batch

        for i in range(len(data_batch)):
            # print("batch", b, "circuit", i)
            job_data = data_batch[i]
            n = job_data["n"]
            L = job_data["L"]
            clifford_layer = job_data["clifford_layer"]
            cliffordOp = Clifford.from_dict(job_data['clifford'])

            '''TODO: add inject Pauli noise according to inject_pauli_sparse
            This should be done by calculated a new CliffordOp
            How efficient can we do that?
            '''
            if inject_pauli_sparse is not None: # let's inject some noise
                inject_circ = QuantumCircuit(n)
                for _k in range(L):
                    # clifford layer
                    if clifford_layer == "Id":
                        pass
                    elif clifford_layer == 'CNOT':
                        ngates = int(n/2)
                        for j in range(ngates):
                            inject_circ.cx(2*j,2*j+1)
                    else:
                        raise ValueError('Clifford layer undefined or not supported.')
                    # pauli layer
                    inject_pauli_sample = np.random.choice(inject_pauli_sparse[0],p=inject_pauli_sparse[1])
                    for j in range(n):
                        pauli_gate_1q(inject_circ,j,inject_pauli_sample[n-j-1]) # hope the ordering is correct

                injectOp = Clifford.from_circuit(inject_circ)
                cliffordOp = cliffordOp.compose(injectOp)

            if use_density_matrix is True:
                pass
            else:
                outcomes = job_data["counts"]

                for pauli_label in pauli_request_list:

                    

                    pauli = Pauli(pauli_label)
#                     assert pauli.to_label()[-n:] == pauli.evolve(cliffordOp.adjoint()).to_label()[-n:] # total clifford -> identity (except for intc_cb)
#                     correction_phase = (-1)**(pauli.phase != pauli.evolve(cliffordOp.adjoint()).phase)

                    F = 0
                    tot = 0 
                    for key, counts in outcomes.items():
                        F_key = 1
                        for j in range(n):
                            pauli_bell = Pauli(([key[j+n]=='1'],[key[j]=='0'])) ### change order if not correct
 
                            phase_1q = (-1)**(Pauli(pauli_label[j]).anticommutes(pauli_bell))
                            F_key *= phase_1q
                        F += F_key * counts
                        tot += counts
                    F = abs(F)/tot

                    fidelity_list[pauli_label][L].append(F)
                

            if repeat is None:
                R = 1
            else:
                R = job_data["repeat"]

            circuit_count += R



    EAB_result = {
        "fidelity_list":                     fidelity_list,
    }
    return EAB_result


def rcs_fit_fun(x, a, alpha):
        #return a * np.exp(-alpha * x)
        return a * (alpha ** x)

def fit_EAB(X, xeb_list):
    Y = [np.mean(xeb_list[L]) for L in X]
    Yerr = [sem(xeb_list[L]) for L in X]
    #print(linregress(X,np.log(Y)))
    
    
    try:
        params, pcov = curve_fit(rcs_fit_fun, X, Y, sigma=Yerr, absolute_sigma=True, p0=[1,1])
        alpha = params[1]
        params_err = np.sqrt(np.diag(pcov))
        alpha_err = params_err[1]

    except RuntimeError:
        alpha = 1.0
        alpha_err = 0.0

    # params, pcov = curve_fit(rcs_fit_fun, X, Y, sigma=Yerr, absolute_sigma=True, p0=[1,1])
    # #params, pcov = curve_fit(rcs_fit_fun, X, Y, absolute_sigma=True, p0=[1,1])


    # print(params)

    return alpha, alpha_err

    print(alpha, alpha_err)


