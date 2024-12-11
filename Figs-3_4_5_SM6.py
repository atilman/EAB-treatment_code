# -*- coding: utf-8 -*-
"""
Created on Fri May  5 15:44:15 2023

@author: AndrewTilman
"""
import numpy as np
import matplotlib.pyplot as plt

" MODEL PARAMETERS"
"Choose a ParamSet to import"
from ParamSet_middle import params_dict  # may have to run these files first for successful import
import functions as fn                   # may have to run these files first for successful import
from scipy.integrate import solve_ivp

"SAVE FIGS???"
saveFigs = True



"SIMTYPES"
#! SimType defines the Subsidy and treatment probs for simulation
    #! 1: no private treatment or public treatment
    #! 2: no private treatment but optimal public treatment
    #! 3: no private subsidies and no public treatment
    #! 4: no private subsidies but optimal public treatment
    #! 5: optimal private subsidies but no public treatment
    #! 6: optimal private subsidies and optimal public treatment


#COLORS
healthy = '#6DB33F'
infested = '#CFCA7A'
dead = '#8e7028'
opac = .75
linWidPub = 3
linWidPriv = 4
linWid = 4

t_max = 50
stepsize = .05

privFrac = 6/10
pubFrac = 1- privFrac
y0 = [99/100,1/100,0]
y0_mix = [y0[0]*privFrac,y0[1]*privFrac,y0[2]*privFrac,y0[0]*pubFrac,y0[1]*pubFrac,y0[2]*pubFrac]

beta = params_dict['beta']
gamma = params_dict['gamma']
alpha = params_dict['alpha']

t_eval = np.arange(0, t_max, stepsize)
args = (beta,gamma,alpha,params_dict)



sMax2=params_dict['sMax2']
"OTHER PARAMS"

###! cost of treatment
c = params_dict['c']


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

a = params_dict['a'] # min value of Delta_o
b = params_dict['b'] # max value of Delta_o

V_m = params_dict['V_m'] # value of a healthy tree to society 
W_m = params_dict['W_m'] # Cost of a dead tree to society


delt_min = params_dict['delt_min'] # min social value of avoiding tree mortaltiy 
delt_max = params_dict['delt_max'] # max social value of avoiding tree mortaltiy 

Delta_m = V_m + W_m # value to society
Delta_m_alt = params_dict['Delta_m_alt']

figName = params_dict['paramSet']





#%%
sol_mix = [0,0,0,0,0,0]
for i in range(6):
    params_dict["simType"] = i+1
    fun_mix = lambda t,y: fn.hid_pub_priv(t,y,*args)
    sol_mix[i] = solve_ivp(fun_mix, [0, t_max], y0_mix ,method='LSODA', t_eval=t_eval)


#%%
# Labels for figs, plotting options, etc.

titles = ['No private treatment, \nNo public treatment','No private subsidies, \nNo public treatment','Optimal private subsidies, \nNo public treatment','No private treatment, \nOptimal public treatment','No private subsidies, \nOptimal public treatment','Optimal private subsidies, \nOptimal public treatment']
labels = ['Healthy','Infested','Dying / dead']
panelAlphas = ['(a)','(b)','(c)','(d)','(e)','(f)']
linWid = [linWidPriv,linWidPriv,linWidPriv,linWidPub,linWidPub,linWidPub]
col = [healthy,infested,dead,healthy,infested,dead]
linSty = ['solid','solid','solid','dashed','dashed','dashed']
opacity = [opac,opac,opac,1,1,1]
textSize = 14


# Dynamics: Figure 2
fig2, axs2 = plt.subplots(2, 3, figsize=(12, 7.66))

for i, ax in enumerate(axs2.flat):
    for j in range(len(labels)):
        ax.plot(sol_mix[i].t, sol_mix[i].y[j]+sol_mix[i].y[j+3],label=labels[j], linewidth=linWid[j], color=col[j])
    ax.set_title(titles[i],fontsize=textSize)  
    ax.set_xlabel('Time (years)',fontsize=textSize)
    ax.set_ylabel('Fraction of trees',fontsize=textSize)
    ax.text(0.02, 1.125, panelAlphas[i], transform=ax.transAxes, va='top', fontsize=textSize, fontweight='bold')
    if i==0:
        ax.legend(bbox_to_anchor=(1, .65),loc='upper right',fontsize=12)


# Dynamics: Figure 3
fig3, axs3 = plt.subplots(2, 3, figsize=(12, 7.66))
labels = ['Healthy-Private','Infested-Private','Dying-Private','Healthy-Public','Infested-Public','Dying-Public']

for i, ax in enumerate(axs3.flat):
    for j in range(len(labels)):
        ax.plot(sol_mix[i].t, sol_mix[i].y[j], label=labels[j],linestyle=linSty[j], linewidth=linWid[j], color=col[j], alpha=opacity[j])
    ax.set_title(titles[i],fontsize=textSize)  
    ax.set_xlabel('Time (years)',fontsize=textSize)
    ax.set_ylabel('Fraction of trees',fontsize=textSize)
    ax.text(0.02, 1.125, panelAlphas[i], transform=ax.transAxes, va='top', fontsize=textSize, fontweight='bold')
    if i==0:
        ax.legend(bbox_to_anchor=(1, .68),loc='upper right',fontsize=12)


fig2.tight_layout()
fig3.tight_layout()

"SAVE FIGS"
if saveFigs == True:
    fig2.savefig("FIGS/FigSM6_"+figName+".png", bbox_inches='tight',dpi=250)
    fig3.savefig("FIGS/Fig4_"+figName+".png", bbox_inches='tight',dpi=250)


#%%
"CALCS OF OPTIMAL SUBS for HEATMAPS"

delta = 0.0021
levs = np.linspace(0, 250, 300)
levsPR = np.linspace(0, 1, 300)
colorMap = 'viridis'



dFrac = np.linspace(0.001, .99, 300) #frac dying / dead trees
iFrac = np.linspace(0.001, .99, 300) # frac infested trees
dFRAC, iFRAC = np.meshgrid(dFrac, iFrac)

sH = np.zeros((len(dFrac),len(iFrac))) #subsidy levels
sI = np.zeros((len(dFrac),len(iFrac)))
sD = np.zeros((len(dFrac),len(iFrac)))
ptH = np.zeros((len(dFrac),len(iFrac))) # Treatment probability for private tree
ptI = np.zeros((len(dFrac),len(iFrac)))
ptD = np.zeros((len(dFrac),len(iFrac)))
ptH_noSub = np.zeros((len(dFrac),len(iFrac))) # Treatment probability W/O subsidies for private tree
ptI_noSub = np.zeros((len(dFrac),len(iFrac)))
ptD_noSub = np.zeros((len(dFrac),len(iFrac)))
ptH_pub = np.zeros((len(dFrac),len(iFrac))) # treatment probability for public tree
ptI_pub = np.zeros((len(dFrac),len(iFrac)))
ptD_pub = np.zeros((len(dFrac),len(iFrac)))
dfSHOW = np.zeros((len(dFrac),len(iFrac))) # a rescaling to make equilateral triangles rather than right triangles for simplex


for i in range(len(dFrac)):
    for j in range(len(iFrac)):
        if dFrac[i] + iFrac[j] > 1:
            sH[i,j]=np.nan
            sI[i,j]=np.nan
            sD[i,j]=np.nan
            ptH[i,j]=np.nan
            ptI[i,j]=np.nan
            ptD[i,j]=np.nan
            ptH_noSub[i,j]=np.nan
            ptI_noSub[i,j]=np.nan
            ptD_noSub[i,j]=np.nan
            ptH_pub[i,j]=np.nan
            ptI_pub[i,j]=np.nan
            ptD_pub[i,j]=np.nan

        else:
            
            muh = fn.m_uh(t_horiz,iFRAC[i,j],beta,gamma)
            mth = fn.m_th(t_horiz,iFRAC[i,j],beta,gamma,eps_h,eps_i)
            mui = fn.m_ui(t_horiz,gamma)
            mti = fn.m_ti(t_horiz,gamma,eps_i)
            mud = 1
            mtd = 1
            
            t_h = muh - mth 
            t_i = mui - mti
            t_d = mud - mtd
            
            hFrac = 1-(dFrac[i] + iFrac[j])
            smuh = fn.sm_uh(t_horiz,iFRAC[i,j],hFrac,beta,gamma)
            smth = fn.sm_th(t_horiz,iFRAC[i,j],hFrac,beta,gamma,alpha,eps_h,eps_i)
            smui = fn.sm_ui(hFrac,beta,gamma)
            smti = fn.sm_ti(hFrac,beta,gamma,alpha,eps_i)
            smud = 0
            smtd = 0

            s_h = (smuh - smth) * spilloverCare
            s_i = (smui - smti) * spilloverCare
            s_d = (smud - smtd) * spilloverCare
            
            
            kHtemp = fn.k_X(t_h,t_i,t_d,\
                                fn.prob_xH('h',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Hh,pr_Hi,pr_Hd),\
                                fn.prob_xH('i',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Hh,pr_Hi,pr_Hd),\
                                fn.prob_xH('d',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Hh,pr_Hi,pr_Hd))
            kItemp = fn.k_X(t_h,t_i,t_d,\
                                fn.prob_xI('h',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Ih,pr_Ii,pr_Id),\
                                fn.prob_xI('i',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Ih,pr_Ii,pr_Id),\
                                fn.prob_xI('d',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Ih,pr_Ii,pr_Id))
            kDtemp = fn.k_X(t_h,t_i,t_d,\
                                fn.prob_xD('h',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Dh,pr_Di,pr_Dd),\
                                fn.prob_xD('i',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Dh,pr_Di,pr_Dd),\
                                fn.prob_xD('d',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Dh,pr_Di,pr_Dd))
                
            LHtemp = fn.L_X(s_h,s_i,s_d,\
                                fn.prob_xH('h',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Hh,pr_Hi,pr_Hd),\
                                fn.prob_xH('i',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Hh,pr_Hi,pr_Hd),\
                                fn.prob_xH('d',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Hh,pr_Hi,pr_Hd))
            LItemp = fn.L_X(s_h,s_i,s_d,\
                                fn.prob_xI('h',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Ih,pr_Ii,pr_Id),\
                                fn.prob_xI('i',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Ih,pr_Ii,pr_Id),\
                                fn.prob_xI('d',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Ih,pr_Ii,pr_Id))
            LDtemp = fn.L_X(s_h,s_i,s_d,\
                                fn.prob_xD('h',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Dh,pr_Di,pr_Dd),\
                                fn.prob_xD('i',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Dh,pr_Di,pr_Dd),\
                                fn.prob_xD('d',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Dh,pr_Di,pr_Dd))
                                

            sH[i,j] = fn.s_opt(c,kHtemp,LHtemp,Delta_m,a,b)
            sI[i,j] = fn.s_opt(c,kItemp,LItemp,Delta_m,a,b)
            sD[i,j] = fn.s_opt(c,kDtemp,LDtemp,Delta_m,a,b)

            ptH[i,j] = fn.pr_treated(c,kHtemp,LHtemp,Delta_m,a,b)
            ptI[i,j] = fn.pr_treated(c,kItemp,LItemp,Delta_m,a,b)
            ptD[i,j] = fn.pr_treated(c,kDtemp,LDtemp,Delta_m,a,b)
            
            ptH_noSub[i,j] = fn.pr_treated(c,kHtemp,LHtemp,0,a,b)
            ptI_noSub[i,j] = fn.pr_treated(c,kItemp,LItemp,0,a,b)
            ptD_noSub[i,j] = fn.pr_treated(c,kDtemp,LDtemp,0,a,b)
            
            ptH_pub[i,j] = fn.pr_treated_pub(c,kHtemp,LHtemp,Delta_m_alt,a,b)
            ptI_pub[i,j] = fn.pr_treated_pub(c,kItemp,LItemp,Delta_m_alt,a,b)
            ptD_pub[i,j] = fn.pr_treated_pub(c,kDtemp,LDtemp,Delta_m_alt,a,b)
            
            dfSHOW[i,j]= dFRAC[i,j] + 1/2*iFRAC[i,j]
#%%
"PLOTS OF OPTIMAL SUBS WITH HEATMAPS"

"LABELS AND SUCH"
rowTitles = ['Assessed \nhealthy','','','','Assessed \ninfested','','','','Assessed \ndying / dead','','','']
colTitles = [' Optimal subsidy \n','Treatment probability \n   w/ optimal subsidy', 'Treatment probability \n   w/o subsidy','Treatment probability \n   for public trees','','','','','','','','','','','','','','','','','','']
heatMapsData = [sH,ptH,ptH_noSub,ptH_pub,sI,ptI,ptI_noSub,ptI_pub,sD,ptD,ptD_noSub,ptD_pub]
cbMax = [sMax2+1,2,2,2,sMax2+1,2,2,2,sMax2+1,2,2,2]
cbSteps = [50,1,1,1,50,1,1,1,50,1,1,1]
levsList = [levs,levsPR,levsPR,levsPR,levs,levsPR,levsPR,levsPR,levs,levsPR,levsPR,levsPR]
panelAlphas = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)']

"PLOTTING"
fig4, bxs = plt.subplots(3, 4, figsize=(15.3, 9))

for i, bx in enumerate(bxs.flat):
    bx.text(.35,1.05,'All infested',size='12')
    bx.text(-.05,-.1,'All healthy',size='12')
    bx.text(.76,-.1,'All dying',size='12')
    bx.text(-.5,.5,rowTitles[i],size='16')
    bx.text(.12,1.25,colTitles[i],size='16')
    bx.set_frame_on(None)
    BX = bx.contourf( dfSHOW,iFRAC,heatMapsData[i],levels=levsList[i],cmap=colorMap)
    bx.set_ylim([0,1])
    bx.set_xlim([0,1])
    fig4.colorbar(BX,ticks=range(0,cbMax[i], cbSteps[i]),ax=bx)
    bx.set_xticks(())
    bx.set_yticks(())
    bx.text(-0.05, 1.08, panelAlphas[i], transform=bx.transAxes, va='top', fontsize=12, fontweight='bold')

fig4.tight_layout()

"SAVE PLOTS?"
if saveFigs == True:
    fig4.savefig("FIGS/Fig3_"+figName+".png", bbox_inches='tight',dpi=250)


#%%
"""
DYNAMICS OF THE INVASION UNDER NO subsidy or pub treatment VS OPTIMAL POLICY
"""


    
t_max=500
stepsize=.05
timeGap=3.5
sims= 8
privFrac = .6
pubFrac = 1- privFrac
y0 = [.99,.01,0]
y0_mix = [y0[0]*privFrac,y0[1]*privFrac,y0[2]*privFrac,y0[0]*pubFrac,y0[1]*pubFrac,y0[2]*pubFrac]
beta = 1
gamma = .3
alpha = 1
t_eval = np.arange(0, t_max, stepsize)
args = (beta,gamma,alpha,params_dict)

"no treatmeents"
params_dict["simType"] = 1
fun_mix_none = lambda t,y: fn.hid_pub_priv(t,y,*args)
sol = solve_ivp(fun_mix_none, [0, t_max], y0_mix ,method='LSODA', t_eval=t_eval)


"no subs / no pub treatment"
params_dict["simType"] = 2
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
fig5 = plt.figure(figsize=(12,3.73))

rx = fig5.add_subplot(131)
plt.text(.315,1.05,'All infested',size='14')
plt.text(-.07,-.1,'All healthy',size='14')
plt.text(.75,-.1,'All dying',size='14')
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
            params_dict["simType"] = 2
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

i_early = 1
i_late = 4


ex = fig5.add_subplot(132)
ex.plot(sol_noSub.t, sol_noSub.y[0]+sol_noSub.y[3],linestyle='dotted',linewidth=2,color=healthy)
ex.plot(sol_noSub.t, sol_noSub.y[1]+sol_noSub.y[4],linestyle='dotted',linewidth=2,color=infested)
ex.plot(sol_noSub.t, sol_noSub.y[2]+sol_noSub.y[5],linestyle='dotted',linewidth=2,color=dead)
ex.plot(sol_opt1[i_early].t, sol_opt1[i_early].y[0]+sol_opt1[i_early].y[3],label='Healthy',linewidth=4,color=healthy)
ex.plot(sol_opt1[i_early].t, sol_opt1[i_early].y[1]+sol_opt1[i_early].y[4],label='Infested',linewidth=4,color=infested)
ex.plot(sol_opt1[i_early].t, sol_opt1[i_early].y[2]+sol_opt1[i_early].y[5],label='Dying / dead',linewidth=4,color=dead)
ex.set_xlim([0,50])
ex.set_xlabel('time (years)',size='14')
ex.set_ylabel('Fraction of trees',size='14')
ex.set_title('Early intervention',size='14')
ex.legend(loc='right')


fx = fig5.add_subplot(133)
fx.plot(sol_noSub.t, sol_noSub.y[0]+sol_noSub.y[3],linestyle='dotted',linewidth=2,color=healthy)
fx.plot(sol_noSub.t, sol_noSub.y[1]+sol_noSub.y[4],linestyle='dotted',linewidth=2,color=infested)
fx.plot(sol_noSub.t, sol_noSub.y[2]+sol_noSub.y[5],linestyle='dotted',linewidth=2,color=dead)
fx.plot(sol_opt1[i_late].t, sol_opt1[i_late].y[0]+sol_opt1[i_late].y[3],label='Healthy',linewidth=4,color=healthy)
fx.plot(sol_opt1[i_late].t, sol_opt1[i_late].y[1]+sol_opt1[i_late].y[4],label='Infested',linewidth=4,color=infested)
fx.plot(sol_opt1[i_late].t, sol_opt1[i_late].y[2]+sol_opt1[i_late].y[5],label='Dying / dead',linewidth=4,color=dead)
fx.set_xlim([0,50])
fx.set_xlabel('time (years)',size='14')
fx.set_ylabel('Fraction of trees',size='14')
fx.set_title('Late intervention',size='14')


rx.text(-0.05, 1.08, '(a)', transform=rx.transAxes, va='top', fontsize=12, fontweight='bold')
ex.text(-0.05, 1.08, '(b)', transform=ex.transAxes, va='top', fontsize=12, fontweight='bold')
fx.text(-0.05, 1.08, '(c)', transform=fx.transAxes, va='top', fontsize=12, fontweight='bold')


plt.tight_layout()

if saveFigs == True:
    fig5.savefig("FIGS/Fig5_"+figName+".png", bbox_inches='tight',dpi=250)