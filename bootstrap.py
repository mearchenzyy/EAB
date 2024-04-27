import random   
import copy
from copy import copy, deepcopy
import numpy as np
from scipy.stats import sem, unitary_group
from scipy.optimize import curve_fit
import EAB_CB_data_analysis

    
#only eliminate one element in bootstrapping
def bootstrap_fit_EAB_plot(X, xeb_list):
    reps=30
    full_list=deepcopy(xeb_list)
    alpha_rtn_frm_fit=[0]*reps
    alpha_err_rtn_frm_fit=[0]*reps
    a_rtn_frm_fit=[0]*reps
    Y_rtn_frm_fit=[0]*reps
    Yerr_rtn_frm_fit=[0]*reps
    for i in range (reps):
        xeb_list_temp=deepcopy(full_list)
        for j in X:
            r=random.randint(0, 19)
#             print ("before random removal", len(xeb_list_temp[j]))
            xeb_list_temp[j].remove(xeb_list_temp[j][r])
#             print ("after random removal", len(xeb_list_temp[j]))
            assert len(xeb_list_temp[j])==19
#         print (i)
        alpha_rtn_frm_fit[i], a_rtn_frm_fit[i],alpha_err_rtn_frm_fit[i],Y_rtn_frm_fit[i],Yerr_rtn_frm_fit[i]= EAB_CB_data_analysis.fit_EAB_plot(X, xeb_list_temp)
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
        alpha_rtn_frm_fit[i], a_rtn_frm_fit[i],alpha_err_rtn_frm_fit[i],Y_rtn_frm_fit[i],Yerr_rtn_frm_fit[i]= EAB_CB_data_analysis.fit_EAB_plot(X, xeb_list_temp)
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
        alpha_rtn_frm_fit[i], a_rtn_frm_fit[i],alpha_err_rtn_frm_fit[i],Y_rtn_frm_fit[i],Yerr_rtn_frm_fit[i]= EAB_CB_data_analysis.fit_EAB_plot(X, xeb_list_temp)
#         print (alpha_rtn_frm_fit[i])
    alpha_bootstrap=np.mean(alpha_rtn_frm_fit)
    alpha_err_bootstrap=np.std(alpha_rtn_frm_fit)
    return alpha_bootstrap,alpha_err_bootstrap, alpha_rtn_frm_fit,a_rtn_frm_fit, alpha_err_rtn_frm_fit,Y_rtn_frm_fit,Yerr_rtn_frm_fit

