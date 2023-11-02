# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 16:22:14 2023

@author: AndrewTilman
"""

import numpy as np
import matplotlib.pyplot as plt
" MODEL PARAMETERS"
"Choose a ParamSet 1 2 3 or 4"
from ParamSet1 import params_dict
import functions as fn
from scipy.integrate import solve_ivp

"SAVE FIGS???"
saveFigs = False

"SIMTYPE"
#! SimType defines the Subsidy and treatment probs for simulation
    #! 1: no private treatment or public treatment
    #! 2: no private treatment but optimal public treatment
    #! 3: no private subsidies and no public treatment
    #! 4: no private subsidies but optimal public treatment
    #! 5: optimal private subsidies but no public treatment
    #! 6: optimal private subsidies and optimal public treatment

"COLORS"
healthy = '#6DB33F'
infested = '#CFCA7A'
dead = '#8e7028'
opac = .75

"PLOT OPTIONS"
uMax = params_dict['uMax']
uMin=params_dict['uMin']
sMax=params_dict['sMax']
sMin = params_dict['sMin']
sMax2=params_dict['sMax2']
linWidPub = 6
linWidPriv = 3
linWid = 4

"SIM OPTIONS"
t_max=500
stepsize=.05
timeGap=3.5
sims= 8
privFrac = .6
pubFrac = 1- privFrac
y0 = [.98,.01,0.01]
y0_mix = [y0[0]*privFrac,y0[1]*privFrac,y0[2]*privFrac,y0[0]*pubFrac,y0[1]*pubFrac,y0[2]*pubFrac]

"EPI PARAMS"
beta = 1
gamma = .3
alpha = 1
t_eval = np.arange(0, t_max, stepsize)
args = (beta,gamma,alpha,params_dict)

"OTHER PARAMS"
###! cost of treatment
c = params_dict['c']

###! survalence params 'h' 'i' and 'd' are the underlying states of the trees
pr_h = params_dict['pr_h'] #frac healthy
pr_i = params_dict['pr_i'] #frac infested
pr_d = params_dict['pr_d'] #frac dead

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

###! outcomes of treatment
hth = params_dict['hth'] # prob of a tree staying healty 'h', given that it is treated 't' and healthy 'h'
huh = params_dict['huh']  # prob of a tree staying healty 'h', given that it is untreated 'u' and healthy 'h'
hti = params_dict['hti'] # prob of a tree becoming healty 'h', given that it is treated 't' and infested 'i'
hui = params_dict['hui'] # prob of a tree becoming healty 'h', given that it is untreated 'u' and infested 'i'
htd = params_dict['htd']   # prob of saving a dead tree -> 0
hud = params_dict['hud']   # prob of saving a dead tree -> 0

###! Tree values
# V_o = value of a healthy tree to its owner 
# W_o = cost of a dead tree to its owner
# Delta_o = V_o + W_o = net benefit of saving a tree to its owner

a = params_dict['a'] # min value of Delta_o
b = params_dict['b'] # max value of Delta_o

V_m = params_dict['V_m'] # value of a healthy tree to society 
W_m = params_dict['W_m'] # Cost of a dead tree to society

vFrac = params_dict['vFrac'] # fraction of social benefit of saving a tree derived from ecosystem service benefit
wFrac = params_dict['wFrac'] # fraction of social beenfit of saving a dree derived from avoiding dealing with dead tree
delt_min = params_dict['delt_min'] # min social value of saving a tree 
delt_max = params_dict['delt_max'] # max social value of saving a tree

Delta_m = V_m + W_m # value of saving a tree to society
Delta_m_alt = params_dict['Delta_m_alt']

figName = params_dict['paramSet']

if vFrac + wFrac != 1:
    print('fix your value fractions! vFrac / wFrac so they add to one!')
if pr_h + pr_i +pr_d != 1:
    print('fix your survailence probs!')
if pr_Hh + pr_Ih + pr_Dh != 1 or pr_Hi + pr_Ii + pr_Di != 1 or pr_Hd + pr_Id + pr_Dd != 1 :
    print('fix your assessment accuracy probs so they sum to 1!')
if np.max((hth,huh,hti,hui,htd,hud))>1 or np.min((hth,huh,hti,hui,htd,hud))<0:
    print('sort out your treatment sucess probablities!')    
    
###! Bayes Theorem: pr_hH is the prob that a tree assessed as 'H' is truly healthy 'h'
pr_hH = fn.prob_xH('h',pr_h,pr_i,pr_d,pr_Hh,pr_Hi,pr_Hd)
pr_iH = fn.prob_xH('i',pr_h,pr_i,pr_d,pr_Hh,pr_Hi,pr_Hd)
pr_dH = fn.prob_xH('d',pr_h,pr_i,pr_d,pr_Hh,pr_Hi,pr_Hd)

pr_hI = fn.prob_xI('h',pr_h,pr_i,pr_d,pr_Ih,pr_Ii,pr_Id)
pr_iI = fn.prob_xI('i',pr_h,pr_i,pr_d,pr_Ih,pr_Ii,pr_Id)
pr_dI = fn.prob_xI('d',pr_h,pr_i,pr_d,pr_Ih,pr_Ii,pr_Id)

pr_hD = fn.prob_xD('h',pr_h,pr_i,pr_d,pr_Dh,pr_Di,pr_Dd)
pr_iD = fn.prob_xD('i',pr_h,pr_i,pr_d,pr_Dh,pr_Di,pr_Dd)
pr_dD = fn.prob_xD('d',pr_h,pr_i,pr_d,pr_Dh,pr_Di,pr_Dd)

###! fraction of trees with each classification
pr_H = pr_Hh * pr_h + pr_Hi * pr_i + pr_Hd * pr_d
pr_I = pr_Ih * pr_h + pr_Ii * pr_i + pr_Id * pr_d
pr_D = pr_Dh * pr_h + pr_Di * pr_i + pr_Dd * pr_d

###! marginal impact of treatment
t_h = hth - huh
t_i = hti - hui
t_d = htd - hud

###! exected impact of treatemt on saving a tree, given its assessed state
k_H = fn.k_X(t_h,t_i,t_d,pr_hH,pr_iH,pr_dH)
k_I = fn.k_X(t_h,t_i,t_d,pr_hI,pr_iI,pr_dI)
k_D = fn.k_X(t_h,t_i,t_d,pr_hD,pr_iD,pr_dD)


"""
DYNAMICS OF THE INVASION UNDER NO subsidy or pub treatment VS OPTIMAL POLICY
"""

"no treatmeents"
params_dict["simType"] = 1
fun_mix_none = lambda t,y: fn.hid_pub_priv(t,y,*args)
sol = solve_ivp(fun_mix_none, [0, t_max], y0_mix ,method='LSODA', t_eval=t_eval)

"no subs / no pub treatment"
params_dict["simType"] = 3
fun_mix_noSub = lambda t,y: fn.hid_pub_priv(t,y,*args)
sol_noSub = solve_ivp(fun_mix_noSub, [0, t_max], y0_mix ,method='LSODA', t_eval=t_eval)

'opt subs / opt pub treatments'
params_dict["simType"] = 6
fun_mix_opt = lambda t,y: fn.hid_pub_priv(t,y,*args)
sol_opt = solve_ivp(fun_mix_opt, [0, t_max], y0_mix ,method='LSODA', t_eval=t_eval)

t_ev = [0,0,0,0,0,0,0,0,0,0]
sol_opt1 = [0,0,0,0,0,0,0,0,0]

for j in range(sims):
    t_ev[j] = np.arange(j*timeGap, t_max, stepsize)
    sol_opt1[j] = solve_ivp(fun_mix_opt, [j*timeGap, t_max], sol_noSub.y[:,int(j*timeGap/stepsize)] ,method='LSODA', t_eval=t_ev[j]) 
    
#%%
# Plotting stream plot
fig4 = plt.figure(figsize=(12,3.73))

rx = fig4.add_subplot(131)
plt.text(.36,1.05,'Infested',size='14')
plt.text(-.05,-.1,'Healthy',size='14')
plt.text(.59,-.1,'Dying / dead',size='14')
plt.box(on=None)

rx.set_xticks([])
rx.set_yticks([])

x = np.linspace(0.01,0.99, 200)
y = np.linspace(0.01, 0.99, 200)
X, Y = np.meshgrid(x, y)
"STREAMPLOT OPT"
dxdt = np.zeros((len(x),len(y)))
dydt = np.zeros((len(x),len(y)))
"STREAMPLOT NO SUB / no pub treatment"
dxdt_a = np.zeros((len(x),len(y)))
dydt_a = np.zeros((len(x),len(y)))
"STREAMPLOT NO treatment"
dxdt_nt = np.zeros((len(x),len(y)))
dydt_nt = np.zeros((len(x),len(y)))

for i in range(len(x)):
    for j in range(len(y)):
        heal = 1 - X[i,j] - 1/2*Y[i,j]
        infe = Y[i,j]
        dea = X[i,j] - 1/2*Y[i,j]
        if heal >= 0 and infe >=0 and dea >=0:
            params_dict["simType"] = 6
            [dhPrivdt,diPrivdt,ddPrivdt,dhPubdt,diPubdt,ddPubdt] = fn.hid_pub_priv(0,[heal/2,infe/2,dea/2,heal/2,infe/2,dea/2],beta,gamma,alpha,params_dict)
            didt = diPrivdt + diPubdt
            dddt = ddPrivdt + ddPubdt
            params_dict["simType"] = 3
            [dhPrivdt,diPrivdt,ddPrivdt,dhPubdt,diPubdt,ddPubdt] = fn.hid_pub_priv(0,[heal/2,infe/2,dea/2,heal/2,infe/2,dea/2],beta,gamma,alpha,params_dict)
            didt_a = diPrivdt + diPubdt
            dddt_a = ddPrivdt + ddPubdt
            params_dict["simType"] = 1
            [dhPrivdt,diPrivdt,ddPrivdt,dhPubdt,diPubdt,ddPubdt] = fn.hid_pub_priv(0,[heal/2,infe/2,dea/2,heal/2,infe/2,dea/2],beta,gamma,alpha,params_dict)
            didt_nt = diPrivdt + diPubdt
            dddt_nt = ddPrivdt + ddPubdt
            dxdt[i,j] = dddt + (1/np.sqrt(3))*didt 
            dydt[i,j] = didt
            dxdt_a[i,j] = dddt_a + (1/np.sqrt(3))*didt_a 
            dydt_a[i,j] = didt_a
            dxdt_nt[i,j] = dddt_nt + (1/np.sqrt(3))*didt_nt 
            dydt_nt[i,j] = didt_nt

            

rx.streamplot(X, Y, dxdt_nt, dydt_nt,color='k',density=4/14,arrowsize=1.5)

dShow = (sol_noSub.y[2]+sol_noSub.y[5])+1/2*(sol_noSub.y[1]+sol_noSub.y[4])
iShow = sol_noSub.y[1]+sol_noSub.y[4]

dShow_untreated = sol.y[2]+sol.y[5]+1/2*(sol.y[1]+sol.y[4])
iShow_untreated = sol.y[1]+sol.y[4]

rx.plot(dShow_untreated,iShow_untreated,color='k',alpha=.85,linewidth=3,label='No treatment')
rx.plot(dShow,iShow,color='dimgrey',alpha=1,linewidth=4,label='No public action')

for j in range(sims):
    dShow_opt1 = sol_opt1[j].y[2]+sol_opt1[j].y[5]+1/2*(sol_opt1[j].y[1]+sol_opt1[j].y[4])
    iShow_opt1 = sol_opt1[j].y[1]+sol_opt1[j].y[4]
    if j==0:
        
        rx.plot(dShow_opt1,iShow_opt1,'silver',alpha=1,linewidth=4,label='Optimal policies')
    else:
        
        rx.plot(dShow_opt1,iShow_opt1,'silver',alpha=1,linewidth=4)


rx.plot(dShow,iShow,color='dimgrey',alpha=1,linewidth=4)

 


for j in range(sims):
    dShow_opt1 = sol_opt1[j].y[2]+sol_opt1[j].y[5]+1/2*(sol_opt1[j].y[1]+sol_opt1[j].y[4])
    iShow_opt1 = sol_opt1[j].y[1]+sol_opt1[j].y[4]
    

rx.legend(loc='upper right', bbox_to_anchor=(0.4, 1))
rx.plot((0,1/2,1,0),(0,1,0,0),'k',linewidth=2)


ex = fig4.add_subplot(132)
ex.plot(sol_noSub.t, sol_noSub.y[0]+sol_noSub.y[3],'--',linewidth=2,color=healthy)
ex.plot(sol_noSub.t, sol_noSub.y[1]+sol_noSub.y[4],'--',linewidth=2,color=infested)
ex.plot(sol_noSub.t, sol_noSub.y[2]+sol_noSub.y[5],'--',linewidth=2,color=dead)
ex.plot(sol_opt1[1].t, sol_opt1[1].y[0]+sol_opt1[1].y[3],label='Healthy',linewidth=4,color=healthy)
ex.plot(sol_opt1[1].t, sol_opt1[1].y[1]+sol_opt1[1].y[4],label='Infested',linewidth=4,color=infested)
ex.plot(sol_opt1[1].t, sol_opt1[1].y[2]+sol_opt1[1].y[5],label='Dying / dead',linewidth=4,color=dead)
ex.set_xlim([0,50])
ex.set_xlabel('time (years)',size='14')
ex.set_ylabel('Fraction of trees',size='14')
ex.set_title('Early intervention',size='14')
ex.legend(loc='right')


fx = fig4.add_subplot(133)
fx.plot(sol_noSub.t, sol_noSub.y[0]+sol_noSub.y[3],'--',linewidth=2,color=healthy)
fx.plot(sol_noSub.t, sol_noSub.y[1]+sol_noSub.y[4],'--',linewidth=2,color=infested)
fx.plot(sol_noSub.t, sol_noSub.y[2]+sol_noSub.y[5],'--',linewidth=2,color=dead)
fx.plot(sol_opt1[3].t, sol_opt1[3].y[0]+sol_opt1[3].y[3],label='Healthy',linewidth=4,color=healthy)
fx.plot(sol_opt1[3].t, sol_opt1[3].y[1]+sol_opt1[3].y[4],label='Infested',linewidth=4,color=infested)
fx.plot(sol_opt1[3].t, sol_opt1[3].y[2]+sol_opt1[3].y[5],label='Dying / dead',linewidth=4,color=dead)
fx.set_xlim([0,50])
fx.set_xlabel('time (years)',size='14')
fx.set_ylabel('Fraction of trees',size='14')
fx.set_title('Late intervention',size='14')
plt.tight_layout()

if saveFigs == True:
    fig4.savefig("FIGS/dynam_simplex_mixed_"+figName+".pdf",bbox_inches='tight', dpi=50)
    fig4.savefig("FIGS/dynam_simplex_mixed_"+figName+".png", bbox_inches='tight',dpi=250)