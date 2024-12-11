# -*- coding: utf-8 -*-
"""
Created on 11 / 5 / 2024

@author: AndrewTilman
"""
import numpy as np
import matplotlib.pyplot as plt
" MODEL PARAMETERS"

"Choose a ParamSets to import"
from ParamSet_middle import params_dict

import functions as fn
from scipy.integrate import solve_ivp

"SAVE FIGS???"
saveFigs = True



"SIMTYPE"
#! SimType defines the Subsidy and treatment probs for simulation
    #! 1: no private treatment or public treatment
    #! 4: no private treatment but optimal public treatment
    #! 2: no private subsidies and no public treatment
    #! 5: no private subsidies but optimal public treatment
    #! 3: optimal private subsidies but no public treatment
    #! 6: optimal private subsidies and optimal public treatment

"COLORS"
healthy = '#6DB33F'
infested = '#CFCA7A'
dead = '#8e7028'
opac = .75
linWidPub = 3
linWidPriv = 4
linWid = 4

t_max = 205
stepsize = .025

privFrac = 6/10
pubFrac = 1- privFrac
y0 = [99/100,1/100,0]
y0_mix = [y0[0]*privFrac,y0[1]*privFrac,y0[2]*privFrac,y0[0]*pubFrac,y0[1]*pubFrac,y0[2]*pubFrac]



figName = ''

beta = params_dict['beta']
gamma = params_dict['gamma']
alpha = params_dict['alpha']
t_eval = np.arange(0, t_max, stepsize)
args = (beta,gamma,alpha,params_dict)

"PLOT OPTIONS"
uMax = params_dict['uMax']
uMin=params_dict['uMin']
sMax=params_dict['sMax']
sMin = params_dict['sMin']
sMax2=params_dict['sMax2']
"OTHER PARAMS"

###! cost of treatment
c = params_dict['c']

# ###! survalence params 'h' 'i' and 'd' are the underlying states of the trees
# pr_h = params_dict['pr_h'] #frac healthy
# pr_i = params_dict['pr_i'] #frac infested
# pr_d = params_dict['pr_d'] #frac dead

###! assessment params: pr_Hh is the prob of assessing a tree as 'H' healthy given that it is truly healthy 'h'
pr_Hh = params_dict['pr_Hh']
pr_Ih = params_dict['pr_Ih']
pr_Dh = params_dict['pr_Dh']

pr_Hi = params_dict['pr_Hi']
pr_Ii = params_dict['pr_Ii']
pr_Di = params_dict['pr_Di']

pr_Hd = params_dict['pr_Hd']
pr_Id = params_dict['pr_Id']
pr_Dd = params_dict['pr_Dd']


t_horiz = params_dict['t_horiz']
eps_h = params_dict['eps_h']
eps_i = params_dict['eps_i']
spilloverCare = params_dict['spill_weight']

###! Tree values

a = params_dict['a'] # min value of Delta_o
b = params_dict['b'] # max value of Delta_o

V_m = params_dict['V_m'] # value of a healthy tree to society 
W_m = params_dict['W_m'] # Cost of a dead tree to society
W_m_alt = params_dict['W_m_alt'] # Cost of a dead tree to society, including removal costs

vFrac = params_dict['vFrac'] # fraction of social benefit of saving a tree derived from ecosystem service benefit
wFrac = params_dict['wFrac'] # fraction of social beenfit of saving a dree derived from avoiding dealing with dead tree
delt_min = params_dict['delt_min'] # min social value of saving a tree 
delt_max = params_dict['delt_max'] # max social value of saving a tree

Delta_m = V_m + W_m # value of saving a tree to society
Delta_m_alt = params_dict['Delta_m_alt']

#%%
sol_mix = [0,0,0,0,0]
simTypes = [1,3,4,6]
for i in range(len(simTypes)):
    params_dict["simType"] = simTypes[i]
    fun_mix = lambda t,y: fn.hid_pub_priv(t,y,*args)
    sol_mix[i] = solve_ivp(fun_mix, [0, t_max], y0_mix ,method='LSODA', t_eval=t_eval)
#%%
'Calc costs and benefits of scenarios (TOTAL COSTS UP UNTIL NOW)'

# fig0, axs = plt.subplots(1, 1, figsize=(6, 6))
# fig0 = plt.figure(figsize=(12,4))
# axs = fig0.add_subplot(131)
# bx = fig0.add_subplot(132)
# cx = fig0.add_subplot(133)
fig1 = plt.figure(figsize=(12,8))
axs1 = fig1.add_subplot(231)
axs2 = fig1.add_subplot(232)
axs3 = fig1.add_subplot(233)
axs4 = fig1.add_subplot(234)
axs5 = fig1.add_subplot(235)
axs6 = fig1.add_subplot(236)
# axs7 = fig1.add_subplot(337)
# axs8 = fig1.add_subplot(338)
# axs9 = fig1.add_subplot(339)

labels = ['  ','No treatment','No priv. subs / No pub. treat','Only private subsidies','Only public treatment','no priv. subs / opt. pub. treatment','Optimal pub. / priv. policies','Treat all assessed healthy and infested']
col=['k','k','k','k','k','k','k']
opac = [.15,.3,.45,.75,.85]
# axs.set_title('Net social per/tree value of public trees \n(total up to current year)')
# axs.set_xlabel('Time (years)')
# axs.set_xlim(0,t_max-5)
# bx.set_title('Net social per/tree value of private trees \n(total up to current year)')
# bx.set_xlabel('Time (years)')
# bx.set_xlim(0,t_max-5)
# cx.set_title('Net social per/tree value \n(total up to current year)')
# cx.set_xlabel('Time (years)')
# cx.set_xlim(0,t_max-5)
axs1.set_title('Net annual per/tree value of public trees')
axs1.set_xlabel('Time (years)')
axs1.set_xlim(3,t_max-5)
axs2.set_title('Net annual per/tree value of private trees')
axs2.set_xlabel('Time (years)')
axs2.set_xlim(3,t_max-5)
axs3.set_title('Net annual per/tree value')
axs3.set_xlabel('Time (years)')
axs3.set_xlim(3,t_max-5)
axs4.set_title('Annual per/tree treatment cost of public trees')
axs4.set_xlabel('Time (years)')
axs4.set_xlim(3,t_max-5)
axs5.set_title('Annual per/tree subsidy program cost')
axs5.set_xlabel('Time (years)')
axs5.set_xlim(3,t_max-5)
axs6.set_title('Annual per/tree benefits')
axs6.set_xlabel('Time (years)')
axs6.set_xlim(3,t_max-5)

axs1.set_ylim(0,350)
axs2.set_ylim(0,350)
axs3.set_ylim(0,350)
axs4.set_ylim(0,90)
axs5.set_ylim(0,90)
axs6.set_ylim(0,350)

for i in range(len(simTypes)):
    params_dict["simType"] = simTypes[i]
    t_series =  sol_mix[i].t # time series data
    h_o = sol_mix[i].y[0]   # healthy owned privately frac data ACTUAL 
    i_o = sol_mix[i].y[1]   # infested owned privately frac data ACTUAL
    d_o = sol_mix[i].y[2]   # dead owned privately frac data ACTUAL
    h_m = sol_mix[i].y[3]   # healthy muni frac data ACTUAL
    i_m = sol_mix[i].y[4]   # infested muni frac data ACTUAL 
    d_m = sol_mix[i].y[5]   # dead muni frac data ACTUAL
    
    H_o = pr_Hh*h_o + pr_Hi*i_o + pr_Hd*d_o    # healthy owned privately frac data ASSESSED
    I_o = pr_Ih*h_o + pr_Ii*i_o + pr_Id*d_o # infested owned privately frac data ASSESSED
    D_o = pr_Dh*h_o + pr_Di*i_o + pr_Dd*d_o   # dead owned privately frac data ASSESSED
    H_m = pr_Hh*h_m + pr_Hi*i_m + pr_Hd*d_m   # healthy muni frac data ASSESSED
    I_m = pr_Ih*h_m + pr_Ii*i_m + pr_Id*d_m   # infested muni frac data ASSESSED
    D_m = pr_Dh*h_m + pr_Di*i_m + pr_Dd*d_m   # dead muni frac data ASSESSED
    
    s_H = np.zeros(len(t_series))
    s_I = np.zeros(len(t_series))
    s_D = np.zeros(len(t_series))
    
    B_annual_m = (V_m/t_horiz * (h_m + i_m))/pubFrac 
    B_annual_o = (V_m/t_horiz * (h_o + i_o))/privFrac
    C_annual_m = np.zeros(len(t_series))
    C_annual_o = np.zeros(len(t_series))
    
    B_tot_m = np.zeros(len(t_series))
    B_tot_o = np.zeros(len(t_series))
    C_tot_m = np.zeros(len(t_series))
    C_tot_o = np.zeros(len(t_series))
    
    P_thm = np.zeros(len(t_series))
    P_tim = np.zeros(len(t_series))
    P_tdm = np.zeros(len(t_series))
    
    P_tho = np.zeros(len(t_series))
    P_tio = np.zeros(len(t_series))
    P_tdo = np.zeros(len(t_series))
    
    P_tHm = np.zeros(len(t_series))
    P_tIm = np.zeros(len(t_series))
    P_tDm = np.zeros(len(t_series))
    
    P_tHo = np.zeros(len(t_series))
    P_tIo = np.zeros(len(t_series))
    P_tDo = np.zeros(len(t_series))
    
    for j in range(len(t_series)):
        B_tot_m[j] = np.trapezoid(B_annual_m[0:j+1],t_series[0:j+1])
        B_tot_o[j] = np.trapezoid(B_annual_o[0:j+1],t_series[0:j+1])
        [P_tho[j],P_tio[j],P_tdo[j],P_thm[j],P_tim[j],P_tdm[j],P_tHo[j],P_tIo[j],P_tDo[j],P_tHm[j],P_tIm[j],P_tDm[j],s_H[j],s_I[j],s_D[j]] = fn.pTreat_pub_priv([h_o[j],i_o[j],d_o[j],h_m[j],i_m[j],d_m[j]], beta, gamma, alpha, params_dict)
    
    C_annual_treat_m = (c/t_horiz * (P_thm*h_m + P_tim*i_m + P_tdm*d_m))/pubFrac
    C_annual_treat_o = 1/t_horiz *(s_H*P_tHo*H_o + s_I*P_tIo*I_o + s_D*P_tDo*D_o)/privFrac
    
    for j in range(len(t_series)):
        C_tot_m[j] = np.trapezoid(C_annual_treat_m[0:j+1],t_series[0:j+1]) + d_m[j]*W_m_alt/pubFrac
        C_tot_o[j] = np.trapezoid(C_annual_treat_o[0:j+1],t_series[0:j+1]) + d_o[j]*W_m/privFrac
        if j != len(t_series)-1:
            C_annual_m[j] = C_annual_treat_m[j] + (d_m[j+1] - d_m[j])/(t_series[j+1]-t_series[j])*W_m_alt/pubFrac
            C_annual_o[j] = C_annual_treat_o[j] + (d_o[j+1] - d_o[j])/(t_series[j+1]-t_series[j])*W_m/pubFrac
    # ax.plot(t_series,B_tot_m,label='Total benefits so far')
    # ax.plot(t_series,C_tot_m,label='Total costs so far')
    # axs.plot(t_series,B_tot_m-C_tot_m,label=labels[simTypes[i]],color=col[i],linewidth=4,alpha=opac[i])
    # bx.plot(t_series,B_tot_o-C_tot_o,label=labels[simTypes[i]],color=col[i],linewidth=4,alpha=opac[i])
    # cx.plot(t_series, (B_tot_m-C_tot_m)*pubFrac + (B_tot_o-C_tot_o)*privFrac,label=labels[simTypes[i]],color=col[i],linewidth=4,alpha=opac[i])
    axs1.plot(t_series,B_annual_m-C_annual_m,label=labels[simTypes[i]],color=col[i],linewidth=4,alpha=opac[i])
    axs2.plot(t_series,B_annual_o-C_annual_o,label=labels[simTypes[i]],color=col[i],linewidth=4,alpha=opac[i])
    axs3.plot(t_series,(B_annual_m-C_annual_m)*pubFrac + (B_annual_o-C_annual_o)*privFrac,label=labels[simTypes[i]],color=col[i],linewidth=4,alpha=opac[i])
    if simTypes[i] in [4,5,6,7]:
        # axs4.plot(t_series,C_annual_treat_m,label=labels[simTypes[i]],color=col[i],linewidth=4,alpha=opac[i])
        # axs4.plot(t_series,C_annual_m-C_annual_treat_m ,label=labels[i],color=healthy,linewidth=4,alpha=opac[i]+.22)
        axs4.plot(t_series,C_annual_treat_m,label=labels[simTypes[i]],color=col[i],linewidth=4,alpha=opac[i])
        
    if simTypes[i] in [3,6,7]:
        axs5.plot(t_series,C_annual_treat_o,label=labels[simTypes[i]],color=col[i],linewidth=4,alpha=opac[i])
        # axs9.plot(t_series,s_I ,label=labels[i],color=infested,linewidth=4,alpha=opac[i]+.22)
        # axs7.plot(t_series,P_tHo ,label=labels[i],color=healthy,linewidth=4,alpha=opac[i]+.22)
        # axs8.plot(t_series,P_tIo ,label=labels[i],color=infested,linewidth=4,alpha=opac[i]+.22)
        # axs9.plot(t_series,P_tDo ,label=labels[i],color=dead,linewidth=4,alpha=opac[i]+.22)
        
        # axs6.plot(t_series,s_D ,label=labels[i],color=dead,linewidth=4,alpha=opac[i])
    # axs6.plot(t_series,B_annual_m ,label=labels[i],color=col[i],linewidth=4,linestyle='dashed',alpha=opac[i])
    axs6.plot(t_series,B_annual_m*pubFrac + B_annual_o*privFrac ,label=labels[simTypes[i]],color=col[i],linewidth=4,alpha=opac[i])
    # axs6.plot(t_series,s_H ,label=labels[i],color=col[i],linewidth=2,alpha=opac[i])
    # axs6.plot(t_series,s_I ,label=labels[i],color=col[i],linewidth=4,alpha=opac[i])
    # axs6.plot(t_series,s_D ,label=labels[i],color=col[i],linewidth=7,alpha=opac[i])

# axs.legend()
# fig0.tight_layout()
fig1.tight_layout()
axs3.legend()
axs4.legend(loc='upper right')
axs5.legend()
# axs6.legend(loc='upper right')

axs1.text(-0.05, 1.14, '(a)', transform=axs1.transAxes, va='top', fontsize=14, fontweight='bold')
axs2.text(-0.05, 1.14, '(b)', transform=axs2.transAxes, va='top', fontsize=14, fontweight='bold')
axs3.text(-0.05, 1.14, '(c)', transform=axs3.transAxes, va='top', fontsize=14, fontweight='bold')
axs4.text(-0.05, 1.14, '(d)', transform=axs4.transAxes, va='top', fontsize=14, fontweight='bold')
axs5.text(-0.05, 1.14, '(e)', transform=axs5.transAxes, va='top', fontsize=14, fontweight='bold')
axs6.text(-0.05, 1.14, '(f)', transform=axs6.transAxes, va='top', fontsize=14, fontweight='bold')

"SAVE FIGS"
if saveFigs == True:
    fig1.savefig("FIGS/FigSM9"+figName+".png", bbox_inches='tight',dpi=250)



