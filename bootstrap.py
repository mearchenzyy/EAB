import sys
import os
import random   
import copy
from copy import copy, deepcopy
import numpy as np
from scipy.stats import sem, unitary_group
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
sys.path.append(r"/Users/yzhu/yzhu_work/gates projects/EAB")
from EAB_process_modified import *
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

def bootstrap_and_generate_param_lists(pauli_request_list,nqubit,depth, raw_fidelity_list,sample_size):
    """ 
    only called by the two plotting functions below
    """
    fidelity_list = {}
    stdev_list = {}
    alpha_detail={}
    a_detail={}
    alpha_error_detail={}
    Y_detail={}
    Y_err_detail={}
    a_BS_dic={}
    Y_BS_dic={}
    Yerr_BS_dic={}
    for pauli_label in pauli_request_list:
        if(pauli_label == 'I'*nqubit):
            fidelity_list[pauli_label] = 1.0
            stdev_list[pauli_label] = 0.0
        else:
            alpha_bootstrap,alpha_err_bootstrap, alpha_rtn_frm_fit,a_rtn_frm_fit, alpha_err_rtn_frm_fit,Y_rtn_frm_fit,Yerr_rtn_frm_fit = bootstrap_fit_EAB_plot_rs(depth, raw_fidelity_list[pauli_label],sample_size)
            fidelity_list[pauli_label] = alpha_bootstrap #avg pauli fidelity from fit
            stdev_list[pauli_label] = alpha_err_bootstrap #error on pauli fidelity from fit
            a_BS_dic[pauli_label]=a_rtn_frm_fit
            Y_BS_dic[pauli_label]=Y_rtn_frm_fit
            Yerr_BS_dic[pauli_label]=Yerr_rtn_frm_fit 
            alpha_detail[pauli_label]=alpha_rtn_frm_fit
            alpha_error_detail[pauli_label]=alpha_err_rtn_frm_fit
    return fidelity_list, stdev_list,a_BS_dic,Y_BS_dic,Yerr_BS_dic,alpha_detail,alpha_error_detail

def bootstrap_and_plot_method_1(axs, pauli_request_list,nqubit,depth, raw_fidelity_list,sample_size=19):
    """
    Method I: use the average of the 10 sets of Pauli fidelities from 10 times of resampling in bootstrapping
    In the generated plot, data points are the average of 10 times of bootstrapping. Parameters for the fitted curve is average of the 10 fitted curves 
    """
    fidelity_list, stdev_list,a_BS_dic,Y_BS_dic,Yerr_BS_dic,alpha_detail,alpha_error_detail= bootstrap_and_generate_param_lists(pauli_request_list,nqubit,depth, raw_fidelity_list,sample_size)

    #figure 1
    # x_c=np.linspace(0,32,num=80)
    # fig, axs = plt.subplots(4, 4)
    # fig.set_figwidth(20)
    # fig.set_figheight(15)
    # fig.subplots_adjust(hspace=0.5,wspace=0.3) 
    # fig.text(0.5, 0.03, 'depth', ha='center',fontsize=20)
    # fig.text(0.05, 0.5, 'mean fidelty', va='center', rotation='vertical',fontsize=20)
    # for i in range (4):
    #     for j in range(4):
    #         pauli_label=pauli_request_list[4*i+j]
    #         if (pauli_label == 'I'*nqubit):
    # #                 fidelity_list[pauli_label] = 1.0
    # #                 stdev_list[pauli_label] = 0.0
    #             pass
    #         else:
    #             axs[i, j].set_xticks(depth)   
    #             for m in range (10):
    #                 axs[i, j].errorbar(depth,Y_BS_dic[pauli_label][m], yerr=Yerr_BS_dic[pauli_label][m], fmt='o',markersize=5)
    #                 axs[i, j].plot(x_c,rcs_fit_fun(x_c,a_BS_dic[pauli_label][m], alpha_detail[pauli_label][m]))
    #             axs[i, j].set_title(pauli_label[::-1])

    #figure 2
    x_c=np.linspace(0,32,num=80)
    # fig, axs = plt.subplots(4, 4)
    # fig.set_figwidth(20)
    # fig.set_figheight(15)
    # fig.subplots_adjust(hspace=0.5,wspace=0.3) 
    # fig.text(0.5, 0.03, 'depth', ha='center',fontsize=20)
    # fig.text(0.05, 0.5, 'mean fidelty', va='center', rotation='vertical',fontsize=20)
    for i in range (4):
        for j in range(4):
            pauli_label=pauli_request_list[4*i+j]
            axs[i, j].set_title(pauli_label[::-1])
            # axs[i, j].set_ylim([0.15,0.95])
            axs[i, j].set_xticks(depth) 
            if (pauli_label == 'I'*nqubit):
                # axs[i, j].set_ylim([0.1,0.9])
                if i!=3:
                    axs[i, j].set_xticks([])
                if j!=0:
                    axs[i, j].set_yticks([]) 
                else:
                    axs[i, j].set_yticks([0.2,0.4,0.6,0.8])
            else:  
                y_plot=[np.mean([Y_BS_dic[pauli_label][m][d] for m in range(10)]) for d in range(len(depth))]
                y_err_plot=[np.mean([Yerr_BS_dic[pauli_label][m][d] for m in range(10)]) for d in range(len(depth))]
                axs[i, j].errorbar(depth,y_plot, yerr=y_err_plot, fmt='o',markersize=5)
                axs[i, j].plot(x_c,rcs_fit_fun(x_c,np.mean([a_BS_dic[pauli_label][m] for m in range(10)]), np.mean([alpha_detail[pauli_label][m] for m in range(10)])))
                axs[i, j].set_title(pauli_label[::-1],y=1.0,pad=-18)
                # axs[i, j].set_ylim([0.1,0.95])
                if i!=3:
                    axs[i, j].set_xticks([])
                if j!=0:
                    axs[i, j].set_yticks([]) 
                else:
                    axs[i, j].set_yticks([0.2,0.4,0.6,0.8])