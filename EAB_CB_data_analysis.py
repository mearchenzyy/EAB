import itertools
import os
import random
import numpy as np
from scipy import optimize
from scipy.optimize import curve_fit
from scipy.stats import sem, entropy, linregress
import bootstrap
import pickle

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




def map_statepop_2_ibm_mapping_no_ancilla(counts,n):
    counts_ibm_mapping=[0 for i in range (2**(n))]
    idx_ibm_mapping=[]
    counts_ibm_mapping_dic={}
    if len(str(n))<2:
        f="00"+str(n)+"b"
    elif len(str(n))<3:
        f="0"+str(n)+"b"
    else:
        raise ValueError("n is too big")
    for gates_idx in range (2**(n)):
        gates_idx_str=format(gates_idx,f)
        ibm_idx=0
        for i in range (n):
            ibm_idx+=2**(i)*int(gates_idx_str[i])
        ibm_idx_str=format(ibm_idx,f)
#         idx_ibm_mapping.append(ibm_idx_str)
        counts_ibm_mapping[ibm_idx]=counts[gates_idx]
        counts_ibm_mapping_dic[ibm_idx_str]=counts[gates_idx]
    return counts_ibm_mapping, counts_ibm_mapping_dic
        


def CB_load_circuit (n,C, depth,pauli_sample_list,circ_path): 
    '''
    TO DO:
    this function is not working. Needs to be fixed.
    load CB circuits into a dictionary called "all_circuits"
    n: number of qubits
    C: sample for each depth
    d= list of depth, for example [2,4,8]
    pauli_sample_list: list of pauli strings to be sampled
    circ_path: path to circuit file

    return the dictionary all_circuits
    '''
    all_circuits={}
    for pauli in pauli_sample_list:
        all_circuits[pauli]={}
    for d in depth:
        all_circuits[pauli][d]=[]

    for f in os.listdir(circ_path):
        if (f.find(".txt")!=-1):
            n=f.find("_")
            n1=f.find(".")
            Plabel=f[n-2:n]
            # print(Plabel)
            dlabel=f[n+3:n1]
            # print (type(dlabel))
            if int(dlabel) in depth:
                file=open(circ_path+f)
                Lines= file.readlines()
        #         c_d8=[]
                count = 0
                for line in Lines:
                    # print(Plabel)
                    # print (dlabel)
                    all_circuits[Plabel][int(dlabel)]=[]
                    all_circuits[Plabel][int(dlabel)].append(line)
                    count += 1
        #         print (count)
                count=0
    
    return all_circuits

# def CB_load_datafiles()


 
def rcs_fit_fun(x, a, alpha):
        #return a * np.exp(-alpha * x)
        return a * (alpha ** x)

def fit_EAB_plot(X, xeb_list):
    Y = [np.mean(xeb_list[L]) for L in X]
    Yerr = [sem(xeb_list[L]) for L in X]
    #print(linregress(X,np.log(Y)))
       
    try:
        params, pcov = curve_fit(rcs_fit_fun, X, Y, sigma=Yerr, absolute_sigma=True, p0=[1,1])
        alpha = params[1]
        a=params[0]
        params_err = np.sqrt(np.diag(pcov))
        alpha_err = params_err[1]

    except RuntimeError:
        alpha = 1.0
        alpha_err = 0.0


    return alpha,a, alpha_err,Y, Yerr

    print(alpha, alpha_err)


def rcs_fit_fun_depth1(x, alpha):
        #return a * np.exp(-alpha * x)
        return alpha ** x

def fit_EAB_depth1(X, xeb_list):
    Y = [np.mean(xeb_list[L]) for L in X]
    Yerr = [sem(xeb_list[L]) for L in X]
    #print(linregress(X,np.log(Y)))
    
    
    try:
        params, pcov = curve_fit(rcs_fit_fun_depth1, X, Y, sigma=Yerr, absolute_sigma=True, p0=[1])
        alpha = params[0]
        params_err = np.sqrt(np.diag(pcov))
        alpha_err = params_err[0]

    except RuntimeError:
        alpha = 1.0
        alpha_err = 0.0


    return alpha, alpha_err

    print(alpha, alpha_err)

def int_to_pauli(i,n):
    p = np.base_repr(i,base=4)
    p = '0'*(n-len(p)) + p
    p = p.replace('0','I').replace('1', 'X').replace('2', 'Y').replace('3', 'Z')
    return p

def commute(p,q):
    c = 1
    n = len(p)
    for i in range(n):
        if p[i] != 'I' and q[i] != 'I':
            if p[i] != q[i]:
                c *= -1
    return c
def fidelity_to_error(pauli_fidelity,n):
    N = 4**n
    pauli_error = {}
    for i in range(N):
        p = int_to_pauli(i,n)
        pauli_error[p] = 0
        for j in range(N):
            q = int_to_pauli(j,n)
            pauli_error[p] += pauli_fidelity[q] * commute(p,q) / N
    return pauli_error

def fidelity_to_error_w_errorbar(pauli_fidelity,pauli_fidelity_errorbar,n):
    '''
    same as the fidelity to error function above but calculates the error bar on Puali err through simple error propogation for addition 
    and subtraction (use standard deviation of sum of squared errors)
    '''
    N = 4**n
    pauli_error = {}
    pauli_error_errorbar = {}
    for i in range(N):
        p = int_to_pauli(i,n)
        pauli_error[p] = 0
        pauli_error_errorbar[p] = 0
        for j in range(N):
            q = int_to_pauli(j,n)
            pauli_error[p] += pauli_fidelity[q] * commute(p,q) / N
            pauli_error_errorbar[p] += (pauli_fidelity_errorbar[q]/N)**2
        pauli_error_errorbar[p] = np.sqrt(pauli_error_errorbar[p])
    return pauli_error,pauli_error_errorbar

def bootstrap_N_save(s,depth,raw_fidelity_list,data_save_path,pauli_request_list,nqubit):
    '''
    s: sample size for bootstrapping
    depth: example [2,4,8]
    data_save_path= "/Volumes/funkflower/Users/Yingyue/Gates_Lab_Suite-master/PauliNoiseEstimation/data/for plotting/scan wait time test/150us/"
    '''
    
    fidelity_list = {}
    stdev_list = {}
    a_detail={}
    Y_detail={}
    Y_err_detail={}
    a_BS_dic={}
    Y_BS_dic={}
    Yerr_BS_dic={}
    alpha_detail={}
    alpha_error_detail={}    
    for pauli_label in pauli_request_list:
        if (pauli_label == 'I'*nqubit):
            alpha_detail[pauli_label] = [1.0]*10
            alpha_error_detail[pauli_label] = [0.0]*10
        else:
            alpha_bootstrap,alpha_err_bootstrap, alpha_rtn_frm_fit,a_rtn_frm_fit, alpha_err_rtn_frm_fit,Y_rtn_frm_fit,Yerr_rtn_frm_fit = bootstrap.bootstrap_fit_EAB_plot_rs(depth, raw_fidelity_list[pauli_label],s)
            alpha_detail[pauli_label]=alpha_rtn_frm_fit
            alpha_error_detail[pauli_label]=alpha_err_rtn_frm_fit
    filename="alpha_detail_d"
    filename_error="alpha_error_detail_d"
    for d in depth:
        filename+=str(d)
        filename_error+=str(d)
    with open(data_save_path+filename, "wb") as fp:
        pickle.dump(alpha_detail, fp)
    with open(data_save_path+filename_error, "wb") as fp:
        pickle.dump(alpha_error_detail, fp)
