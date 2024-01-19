import random   
import copy
from copy import copy, deepcopy
import numpy as np
from scipy.stats import sem, unitary_group
from scipy.optimize import curve_fit
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

    # params, pcov = curve_fit(rcs_fit_fun, X, Y, sigma=Yerr, absolute_sigma=True, p0=[1,1])
    # #params, pcov = curve_fit(rcs_fit_fun, X, Y, absolute_sigma=True, p0=[1,1])


    # print(params)

    return alpha,a, alpha_err,Y, Yerr

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
    
#only eliminate one element in bootstrapping
def bootstrap_fit_EAB_plot(X, xeb_list):
    full_list=deepcopy(xeb_list)
    alpha_rtn_frm_fit=[0]*10
    alpha_err_rtn_frm_fit=[0]*10
    a_rtn_frm_fit=[0]*10
    Y_rtn_frm_fit=[0]*10
    Yerr_rtn_frm_fit=[0]*10
    for i in range (10):
        xeb_list_temp=deepcopy(full_list)
        for j in X:
            r=random.randint(0, 19)
#             print ("before random removal", len(xeb_list_temp[j]))
            xeb_list_temp[j].remove(xeb_list_temp[j][r])
#             print ("after random removal", len(xeb_list_temp[j]))
            assert len(xeb_list_temp[j])==19
#         print (i)
        alpha_rtn_frm_fit[i], a_rtn_frm_fit[i],alpha_err_rtn_frm_fit[i],Y_rtn_frm_fit[i],Yerr_rtn_frm_fit[i]= fit_EAB_plot(X, xeb_list_temp)
#         print (alpha_rtn_frm_fit[i])
    alpha_bootstrap=np.mean(alpha_rtn_frm_fit)
    alpha_err_bootstrap=np.std(alpha_rtn_frm_fit)
    return alpha_bootstrap,alpha_err_bootstrap, alpha_rtn_frm_fit,a_rtn_frm_fit, alpha_err_rtn_frm_fit,Y_rtn_frm_fit,Yerr_rtn_frm_fit
    
#specify number of elements to resample : rs
def bootstrap_fit_EAB_plot_rs(X, xeb_list,rs):
    full_list=deepcopy(xeb_list)
    alpha_rtn_frm_fit=[0]*10
    alpha_err_rtn_frm_fit=[0]*10
    a_rtn_frm_fit=[0]*10
    Y_rtn_frm_fit=[0]*10
    Yerr_rtn_frm_fit=[0]*10
    for i in range (10):
        xeb_list_temp=deepcopy(full_list)
        for j in X:
#             print ("before random removal", len(xeb_list_temp[j]))
            for m in range (20-rs):
                r=random.randint(0,20-m-1)
#                 r=random.sample(range(20), 20-rs)
                xeb_list_temp[j].remove(xeb_list_temp[j][r])
#             print ("after random removal", len(xeb_list_temp[j]))
            assert len(xeb_list_temp[j])==rs
#         print (i)
        alpha_rtn_frm_fit[i], a_rtn_frm_fit[i],alpha_err_rtn_frm_fit[i],Y_rtn_frm_fit[i],Yerr_rtn_frm_fit[i]= fit_EAB_plot(X, xeb_list_temp)
#         print (alpha_rtn_frm_fit[i])
    alpha_bootstrap=np.mean(alpha_rtn_frm_fit)
    alpha_err_bootstrap=np.std(alpha_rtn_frm_fit)
    return alpha_bootstrap,alpha_err_bootstrap, alpha_rtn_frm_fit,a_rtn_frm_fit, alpha_err_rtn_frm_fit,Y_rtn_frm_fit,Yerr_rtn_frm_fit


def bootstrap_fit_EAB_plot_rs(X, xeb_list,rs):
    full_list=deepcopy(xeb_list)        
    alpha_rtn_frm_fit=[0]*10
    alpha_err_rtn_frm_fit=[0]*10
    a_rtn_frm_fit=[0]*10
    Y_rtn_frm_fit=[0]*10
    Yerr_rtn_frm_fit=[0]*10
    for i in range (10):
        xeb_list_temp=deepcopy(full_list)
        for j in X:
#             print ("before random removal", len(xeb_list_temp[j]))
            for m in range (20-rs):
                r=random.randint(0,20-m-1)
#                 r=random.sample(range(20), 20-rs)
                xeb_list_temp[j].remove(xeb_list_temp[j][r])
#             print ("after random removal", len(xeb_list_temp[j]))
            assert len(xeb_list_temp[j])==rs
#         print (i)
        alpha_rtn_frm_fit[i], a_rtn_frm_fit[i],alpha_err_rtn_frm_fit[i],Y_rtn_frm_fit[i],Yerr_rtn_frm_fit[i]= fit_EAB_plot(X, xeb_list_temp)
#         print (alpha_rtn_frm_fit[i])
    alpha_bootstrap=np.mean(alpha_rtn_frm_fit)
    alpha_err_bootstrap=np.std(alpha_rtn_frm_fit)
    return alpha_bootstrap,alpha_err_bootstrap, alpha_rtn_frm_fit,a_rtn_frm_fit, alpha_err_rtn_frm_fit,Y_rtn_frm_fit,Yerr_rtn_frm_fit

