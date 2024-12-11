# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 07:57:52 2023

@author: AndrewTilman
"""
import numpy as np
"""
BAYES THEOREM FOR ASSESSMENT: "H", "I", and "D" are the assessed states
                                    "x \in {h,i,d} is the underlying state"
                                    probs depend on prevalence of states and assessment accuracies
""" 

def prob_xH(x,pr_h,pr_i,pr_d,pr_Hh,pr_Hi,pr_Hd):
    if x =='h':
        bayes = (pr_Hh * pr_h) / (pr_Hh * pr_h + pr_Hi * pr_i + pr_Hd * pr_d)
    elif x=='i':
        bayes = (pr_Hi * pr_i) / (pr_Hh * pr_h + pr_Hi * pr_i + pr_Hd * pr_d)
    elif x=='d':
        bayes = (pr_Hd * pr_d) / (pr_Hh * pr_h + pr_Hi * pr_i + pr_Hd * pr_d)
    else:
        print('Choose a valid string, h, i ,d, to correspond to the underlying health state')
    return bayes

def prob_xI(x,pr_h,pr_i,pr_d,pr_Ih,pr_Ii,pr_Id):
    if x =='h':
        bayes = (pr_Ih * pr_h) / (pr_Ih * pr_h + pr_Ii * pr_i + pr_Id * pr_d)
    elif x=='i':
        bayes = (pr_Ii * pr_i) / (pr_Ih * pr_h + pr_Ii * pr_i + pr_Id * pr_d)
    elif x=='d':
        bayes = (pr_Id * pr_d) / (pr_Ih * pr_h + pr_Ii * pr_i + pr_Id * pr_d)
    else:
        print('Choose a valid string, h, i ,d, to correspond to the underlying health state')
    return bayes

def prob_xD(x,pr_h,pr_i,pr_d,pr_Dh,pr_Di,pr_Dd):
    if x =='h':
        bayes = (pr_Dh * pr_h) / (pr_Dh * pr_h + pr_Di * pr_i + pr_Dd * pr_d)
    elif x=='i':
        bayes = (pr_Di * pr_i) / (pr_Dh * pr_h + pr_Di * pr_i + pr_Dd * pr_d)
    elif x=='d':
        bayes = (pr_Dd * pr_d) / (pr_Dh * pr_h + pr_Di * pr_i + pr_Dd * pr_d)
    else:
        print('Choose a valid string, h, i ,d, to correspond to the underlying health state')
    return bayes

"OPTIMAL SUBSIDY"
def s_opt(c,k_Phi,L_Phi,Delta_m,a,b): #This is for ecosystem service optimization
    if c <= a*k_Phi:
        s = 0
    elif c + b*k_Phi - 2*a*k_Phi <= Delta_m*(k_Phi+L_Phi):
        s = c - a*k_Phi
    elif np.abs(c - b*k_Phi) < Delta_m*(k_Phi+L_Phi):
        s = (1/2)*(Delta_m*(k_Phi+L_Phi) + c - b*k_Phi)
    else:
        s = 0
    return s


"TREATMET PROBABILITY UNDER NO SUBSIDY"
def pr_treated_nosub(c,k_Phi,a,b):
    if c <= a*k_Phi:
        pr_T = 1
    elif c < b*k_Phi :
        pr_T = (b*k_Phi - c)/(k_Phi*(b-a))
    elif b*k_Phi <= c:
        pr_T = 0
    else:
        # print('theres an error somewhere!')
        pr_T = 0
    return pr_T

"TREATMET PROBABILITY UNDER OPTIMAL SUBSIDY, basaed on assessed states"
def pr_treated(c,k_Phi,L_Phi,Delta_m,a,b):
    if c <= a*k_Phi:
        pr_T = 1
    elif c + b*k_Phi - 2*a*k_Phi <= Delta_m*(k_Phi+L_Phi):
        pr_T = 1
    elif np.abs(c - b*k_Phi) < Delta_m*(k_Phi+L_Phi):
        pr_T = (Delta_m*(k_Phi+L_Phi) - c + b*k_Phi)/(2*k_Phi*(b-a))
    elif c < b*k_Phi and Delta_m*(k_Phi+L_Phi) <= b*k_Phi - c:
        pr_T = (b*k_Phi - c)/(k_Phi*(b-a))
    elif b*k_Phi <= c and Delta_m*(k_Phi+L_Phi) <= c - b*k_Phi:
        pr_T = 0
    else:
        # print('theres an error somewhere!')
        pr_T = 0
    return pr_T


"TREATMENT OF PUBLICLY OWNED TREES"
def pr_treated_pub(c,k_Phi,L_Phi,Delta_m,a,b):
    if c <= Delta_m*(k_Phi+L_Phi):
        pr_T = 1
    else:
        pr_T = 0
    return pr_T


"SURVIVAL PROBABILITY OF A focal TREE ASSESSED IN STATE X GIVEN treatment probability"
def pr_survX(pr_treat,pr_hX,pr_iX,pr_dX,hth,huh,hti,hui,htd,hud):
    prSurv = pr_treat*(hth*pr_hX + hti*pr_iX + htd*pr_dX) + (1-pr_treat)*(huh*pr_hX + hui*pr_iX + hud*pr_dX)
    return prSurv



"MORTALITY RISK PERCEPTION BY A TREE OWNER / MUNI"
"for a healthy untreated tree"
def m_uh(t_horiz,I0,beta,gamma):
    muh = (beta*I0*(1-np.exp(-gamma*t_horiz))-gamma*(1-np.exp(-beta*I0*t_horiz)))/(beta*I0-gamma)
    return muh

"for an infested untreated tree"
def m_ui(t_horiz,gamma):
    mui = 1-np.exp(-gamma*t_horiz)
    return mui

"for a healthy treated tree"
def m_th(t_horiz,I0,beta,gamma,eps_h,eps_i):
    ri = 1 - eps_i
    rh = 1 - eps_h
    mth = (beta*I0*rh*(1-np.exp(-gamma*ri*t_horiz))-gamma*ri*(1-np.exp(-beta*I0*rh*t_horiz)))/(beta*rh*I0-gamma*ri)
    return mth

'for an infested treated tree'
def m_ti(t_horiz,gamma,eps_i):
    ri = 1 - eps_i
    mti = 1 - np.exp(-gamma*ri*t_horiz)
    return mti

"SPILLOVER MORTALITY FROM INFESTATION ON TREE COMMUNITY"
"from an untreated healthy tree"
def sm_uh(t_horiz,I0,H0,beta,gamma):
    Pr_i = 1-np.exp(-beta*I0*t_horiz)
    R_eff = beta*H0/gamma
    return Pr_i*R_eff

"from a treated healthy tree"
def sm_th(t_horiz,I0,H0,beta,gamma,alpha,eps_h,eps_i):
    Pr_i = 1-np.exp(-beta*(1-eps_h)*I0*t_horiz)
    R_eff = beta*H0/(gamma*(1-eps_i)+alpha*eps_i)
    return Pr_i*R_eff

"from an untreated infested tree (R_effective)"
def sm_ui(H0,beta,gamma):
    R_eff = beta*H0/gamma
    return R_eff

"from a treated infested tree"
def sm_ti(H0,beta,gamma,alpha,eps_i):
    R_eff = beta*H0/(gamma*(1-eps_i)+alpha*eps_i)
    return R_eff

"Increase in the survival prob of a *focal* TREE in response to TREATMENT GIVEN IT'S ASSESSED STATE IS X"
def k_X(t_h,t_i,t_d,pr_hX,pr_iX,pr_dX):
    treatmentImpact = pr_hX * t_h + pr_iX * t_i + pr_dX * t_d 
    return treatmentImpact

"SPILLOVER SURVIVAL BENEFIT: Expected increse (# of trees) in tree survival of non-focal trees in response to treatment GIVEN the treated tree's ASSESSED STATE IS X"
def L_X(s_h,s_i,s_d,pr_hX,pr_iX,pr_dX):
    treatmentImpact = pr_hX * s_h + pr_iX * s_i + pr_dX * s_d 
    return treatmentImpact

"""
FOR DYNAMICAL SIMULATIONS:
    DISEASE DYNAMICS EQUATIONS
"""

"treatment probs calculator"
def pTreat_pub_priv(y,beta,gamma,alpha,P):
    
    pr_h_pub = y[3]
    pr_i_pub = y[4]
    pr_d_pub = y[5]
    
    pr_h_priv = y[0]
    pr_i_priv = y[1]
    pr_d_priv = y[2]
    

    pr_h = pr_h_pub + pr_h_priv
    pr_i = pr_i_pub + pr_i_priv
    pr_d = pr_d_pub + pr_d_priv
    
    ###! Subsidy and treatment probs for simulation
        ###! 1: no private treatment or public treatment
        ###! 4: no private treatment but optimal public treatment
        ###! 2: no private subsidies and no public treatment
        ###! 5: no private subsidies but optimal public treatment
        ###! 3: optimal private subsidies but no public treatment
        ###! 6: optimal private subsidies and optimal public treatment
        ###! 7: treat every tree assessed as healthy or infested
        
    simType = P['simType']
    ###! cost of treatment
    c = P['c']
    ###! VALUE FO SAVING A TREE TO OWNER (uniformly distributed from a to b) AND MUNI
    a = P['a']
    b = P['b']
    Delta_m = P['Delta_m']
    Delta_m_alt = P['Delta_m_alt']
    
    ###! assessment params: pr_Hh is the prob of assessing a tree as 'H' healthy given that it is truly healthy 'h'
    pr_Hh = P['pr_Hh']
    pr_Ih = P['pr_Ih']
    pr_Dh = P['pr_Dh']

    pr_Hi = P['pr_Hi']
    pr_Ii = P['pr_Ii']
    pr_Di = P['pr_Di']

    pr_Hd = P['pr_Hd']
    pr_Id = P['pr_Id']
    pr_Dd = P['pr_Dd']


    t_horiz = P['t_horiz']
    eps_h = P['eps_h']
    eps_i = P['eps_i']
    spilloverCare = P['spill_weight']
    
    muh = m_uh(t_horiz,pr_i,beta,gamma)
    mth = m_th(t_horiz,pr_i,beta,gamma,eps_h,eps_i)
    mui = m_ui(t_horiz,gamma)
    mti = m_ti(t_horiz,gamma,eps_i)
    mud = 1
    mtd = 1
    
    t_h = muh - mth 
    t_i = mui - mti
    t_d = mud - mtd
    
    smuh = sm_uh(t_horiz,pr_i,pr_h,beta,gamma)
    smth = sm_th(t_horiz,pr_i,pr_h,beta,gamma,alpha,eps_h,eps_i)
    smui = sm_ui(pr_h,beta,gamma)
    smti = sm_ti(pr_h,beta,gamma,alpha,eps_i)
    smud = 0
    smtd = 0
    
    s_h = (smuh - smth) * spilloverCare
    s_i = (smui - smti) * spilloverCare
    s_d = (smud - smtd) * spilloverCare
    

    # Probabiblity of being in actual state (h,i,d) given assessed state (H,I,D).
    pr_hH = prob_xH('h',pr_h,pr_i,pr_d,pr_Hh,pr_Hi,pr_Hd)
    pr_iH = prob_xH('i',pr_h,pr_i,pr_d,pr_Hh,pr_Hi,pr_Hd)
    pr_dH = prob_xH('d',pr_h,pr_i,pr_d,pr_Hh,pr_Hi,pr_Hd)
    
    pr_hI = prob_xI('h',pr_h,pr_i,pr_d,pr_Ih,pr_Ii,pr_Id)
    pr_iI = prob_xI('i',pr_h,pr_i,pr_d,pr_Ih,pr_Ii,pr_Id)
    pr_dI = prob_xI('d',pr_h,pr_i,pr_d,pr_Ih,pr_Ii,pr_Id)
    
    pr_hD = prob_xD('h',pr_h,pr_i,pr_d,pr_Dh,pr_Di,pr_Dd)
    pr_iD = prob_xD('i',pr_h,pr_i,pr_d,pr_Dh,pr_Di,pr_Dd)
    pr_dD = prob_xD('d',pr_h,pr_i,pr_d,pr_Dh,pr_Di,pr_Dd)
    
    # Change in PROBABILITY OF TREE survival given TREATMENT and IT'S ASSESSED STATE
    k_H = k_X(t_h,t_i,t_d,pr_hH,pr_iH,pr_dH)
    k_I = k_X(t_h,t_i,t_d,pr_hI,pr_iI,pr_dI)
    k_D = k_X(t_h,t_i,t_d,pr_hD,pr_iD,pr_dD)
    
    L_H = L_X(s_h,s_i,s_d,pr_hH,pr_iH,pr_dH)
    L_I = L_X(s_h,s_i,s_d,pr_hI,pr_iI,pr_dI)
    L_D = L_X(s_h,s_i,s_d,pr_hD,pr_iD,pr_dD)
    
    if simType == 1:
        pr_th_priv = 0
        pr_ti_priv = 0
        pr_td_priv = 0
        pr_th_pub = 0
        pr_ti_pub = 0
        pr_td_pub = 0
        pt_tH_priv = 0
        pt_tI_priv = 0
        pt_tD_priv = 0
        pt_tH_pub = 0
        pt_tI_pub = 0
        pt_tD_pub = 0
        s_H = 0
        s_I = 0
        s_D = 0
    elif simType == 4:
        pr_th_priv = 0
        pr_ti_priv = 0
        pr_td_priv = 0
        pr_th_pub = pr_Hh*pr_treated_pub(c,k_H,L_H,Delta_m_alt,a,b) +pr_Ih*pr_treated_pub(c,k_I,L_I,Delta_m_alt,a,b) + pr_Dh*pr_treated_pub(c,k_D,L_D,Delta_m_alt,a,b)
        pr_ti_pub = pr_Hi*pr_treated_pub(c,k_H,L_H,Delta_m_alt,a,b) +pr_Ii*pr_treated_pub(c,k_I,L_I,Delta_m_alt,a,b) + pr_Di*pr_treated_pub(c,k_D,L_D,Delta_m_alt,a,b)
        pr_td_pub = pr_Hd*pr_treated_pub(c,k_H,L_H,Delta_m_alt,a,b) +pr_Id*pr_treated_pub(c,k_I,L_I,Delta_m_alt,a,b) + pr_Dd*pr_treated_pub(c,k_D,L_D,Delta_m_alt,a,b)
        pt_tH_priv = 0
        pt_tI_priv = 0
        pt_tD_priv = 0
        pt_tH_pub = pr_treated_pub(c,k_H,L_H,Delta_m_alt,a,b)
        pt_tI_pub = pr_treated_pub(c,k_I,L_I,Delta_m_alt,a,b)
        pt_tD_pub = pr_treated_pub(c,k_D,L_D,Delta_m_alt,a,b)
        s_H = 0
        s_I = 0
        s_D = 0
    elif simType == 2:
        pr_th_priv = pr_Hh*pr_treated_nosub(c,k_H,a,b) +pr_Ih*pr_treated_nosub(c,k_I,a,b) + pr_Dh*pr_treated_nosub(c,k_D,a,b)
        pr_ti_priv = pr_Hi*pr_treated_nosub(c,k_H,a,b) +pr_Ii*pr_treated_nosub(c,k_I,a,b) + pr_Di*pr_treated_nosub(c,k_D,a,b)
        pr_td_priv = pr_Hd*pr_treated_nosub(c,k_H,a,b) +pr_Id*pr_treated_nosub(c,k_I,a,b) + pr_Dd*pr_treated_nosub(c,k_D,a,b)
        pr_th_pub = 0
        pr_ti_pub = 0
        pr_td_pub = 0
        pt_tH_priv = pr_treated_nosub(c,k_H,a,b)
        pt_tI_priv = pr_treated_nosub(c,k_I,a,b)
        pt_tD_priv = pr_treated_nosub(c,k_D,a,b)
        pt_tH_pub = 0
        pt_tI_pub = 0
        pt_tD_pub = 0
        s_H = 0
        s_I = 0
        s_D = 0
    elif simType == 5:
        pr_th_priv = pr_Hh*pr_treated_nosub(c,k_H,a,b) +pr_Ih*pr_treated_nosub(c,k_I,a,b) + pr_Dh*pr_treated_nosub(c,k_D,a,b)
        pr_ti_priv = pr_Hi*pr_treated_nosub(c,k_H,a,b) +pr_Ii*pr_treated_nosub(c,k_I,a,b) + pr_Di*pr_treated_nosub(c,k_D,a,b)
        pr_td_priv = pr_Hd*pr_treated_nosub(c,k_H,a,b) +pr_Id*pr_treated_nosub(c,k_I,a,b) + pr_Dd*pr_treated_nosub(c,k_D,a,b)
        pr_th_pub = pr_Hh*pr_treated_pub(c,k_H,L_H,Delta_m_alt,a,b) +pr_Ih*pr_treated_pub(c,k_I,L_I,Delta_m_alt,a,b) + pr_Dh*pr_treated_pub(c,k_D,L_D,Delta_m_alt,a,b)
        pr_ti_pub = pr_Hi*pr_treated_pub(c,k_H,L_H,Delta_m_alt,a,b) +pr_Ii*pr_treated_pub(c,k_I,L_I,Delta_m_alt,a,b) + pr_Di*pr_treated_pub(c,k_D,L_D,Delta_m_alt,a,b)
        pr_td_pub = pr_Hd*pr_treated_pub(c,k_H,L_H,Delta_m_alt,a,b) +pr_Id*pr_treated_pub(c,k_I,L_I,Delta_m_alt,a,b) + pr_Dd*pr_treated_pub(c,k_D,L_D,Delta_m_alt,a,b)
        pt_tH_priv = pr_treated_nosub(c,k_H,a,b)
        pt_tI_priv = pr_treated_nosub(c,k_I,a,b)
        pt_tD_priv = pr_treated_nosub(c,k_D,a,b)
        pt_tH_pub = pr_treated_pub(c,k_H,L_H,Delta_m_alt,a,b)
        pt_tI_pub = pr_treated_pub(c,k_I,L_I,Delta_m_alt,a,b)
        pt_tD_pub = pr_treated_pub(c,k_D,L_D,Delta_m_alt,a,b)
        s_H = 0
        s_I = 0
        s_D = 0
    elif simType == 3:
        pr_th_priv = pr_Hh*pr_treated(c,k_H,L_H,Delta_m,a,b) +pr_Ih*pr_treated(c,k_I,L_I,Delta_m,a,b) + pr_Dh*pr_treated(c,k_D,L_D,Delta_m,a,b)
        pr_ti_priv = pr_Hi*pr_treated(c,k_H,L_H,Delta_m,a,b) +pr_Ii*pr_treated(c,k_I,L_I,Delta_m,a,b) + pr_Di*pr_treated(c,k_D,L_D,Delta_m,a,b)
        pr_td_priv = pr_Hd*pr_treated(c,k_H,L_H,Delta_m,a,b) +pr_Id*pr_treated(c,k_I,L_I,Delta_m,a,b) + pr_Dd*pr_treated(c,k_D,L_D,Delta_m,a,b)
        pr_th_pub = 0
        pr_ti_pub = 0
        pr_td_pub = 0
        pt_tH_priv = pr_treated(c,k_H,L_H,Delta_m,a,b)
        pt_tI_priv = pr_treated(c,k_I,L_I,Delta_m,a,b)
        pt_tD_priv = pr_treated(c,k_D,L_D,Delta_m,a,b)
        pt_tH_pub = 0
        pt_tI_pub = 0
        pt_tD_pub = 0
        s_H = s_opt(c,k_H,L_H,Delta_m,a,b)
        s_I = s_opt(c,k_I,L_I,Delta_m,a,b)
        s_D = s_opt(c,k_D,L_D,Delta_m,a,b)
    elif simType == 6:
        pr_th_priv = pr_Hh*pr_treated(c,k_H,L_H,Delta_m,a,b) +pr_Ih*pr_treated(c,k_I,L_I,Delta_m,a,b) + pr_Dh*pr_treated(c,k_D,L_D,Delta_m,a,b)
        pr_ti_priv = pr_Hi*pr_treated(c,k_H,L_H,Delta_m,a,b) +pr_Ii*pr_treated(c,k_I,L_I,Delta_m,a,b) + pr_Di*pr_treated(c,k_D,L_D,Delta_m,a,b)
        pr_td_priv = pr_Hd*pr_treated(c,k_H,L_H,Delta_m,a,b) +pr_Id*pr_treated(c,k_I,L_I,Delta_m,a,b) + pr_Dd*pr_treated(c,k_D,L_D,Delta_m,a,b)
        pr_th_pub = pr_Hh*pr_treated_pub(c,k_H,L_H,Delta_m_alt,a,b) +pr_Ih*pr_treated_pub(c,k_I,L_I,Delta_m_alt,a,b) + pr_Dh*pr_treated_pub(c,k_D,L_D,Delta_m_alt,a,b)
        pr_ti_pub = pr_Hi*pr_treated_pub(c,k_H,L_H,Delta_m_alt,a,b) +pr_Ii*pr_treated_pub(c,k_I,L_I,Delta_m_alt,a,b) + pr_Di*pr_treated_pub(c,k_D,L_D,Delta_m_alt,a,b)
        pr_td_pub = pr_Hd*pr_treated_pub(c,k_H,L_H,Delta_m_alt,a,b) +pr_Id*pr_treated_pub(c,k_I,L_I,Delta_m_alt,a,b) + pr_Dd*pr_treated_pub(c,k_D,L_D,Delta_m_alt,a,b)
        pt_tH_priv = pr_treated(c,k_H,L_H,Delta_m,a,b)
        pt_tI_priv = pr_treated(c,k_I,L_I,Delta_m,a,b)
        pt_tD_priv = pr_treated(c,k_D,L_D,Delta_m,a,b)
        pt_tH_pub = pr_treated_pub(c,k_H,L_H,Delta_m_alt,a,b)
        pt_tI_pub = pr_treated_pub(c,k_I,L_I,Delta_m_alt,a,b)
        pt_tD_pub = pr_treated_pub(c,k_D,L_D,Delta_m_alt,a,b)
        s_H = s_opt(c,k_H,L_H,Delta_m,a,b)
        s_I = s_opt(c,k_I,L_I,Delta_m,a,b)
        s_D = s_opt(c,k_D,L_D,Delta_m,a,b)
    else:
        pr_th_priv = pr_Hh + pr_Ih 
        pr_ti_priv = pr_Hi + pr_Ii
        pr_td_priv = pr_Hd + pr_Id
        pr_th_pub = pr_Hh + pr_Ih 
        pr_ti_pub = pr_Hi + pr_Ii
        pr_td_pub = pr_Hd + pr_Id
        pt_tH_priv = 1
        pt_tI_priv = 1
        pt_tD_priv = 0
        pt_tH_pub = 1
        pt_tI_pub = 1
        pt_tD_pub = 0
        s_H = c - a * k_H
        s_I = c - a * k_I
        s_D = 0
        
    return (pr_th_priv,pr_ti_priv,pr_td_priv,pr_th_pub,pr_ti_pub,pr_td_pub,pt_tH_priv,pt_tI_priv,pt_tD_priv,pt_tH_pub,pt_tI_pub,pt_tD_pub,s_H,s_I,s_D)

"DYNAMICS in a landscape of mixed pub and priv ownership"
def hid_pub_priv(t,y,beta,gamma,alpha,P):
    h_pub = y[3]
    i_pub = y[4]
    # d_pub = y[5]
    
    h_priv = y[0]
    i_priv = y[1]
    # d_priv = y[2]
    
    pr_h_pub = y[3]
    pr_i_pub = y[4]
    pr_d_pub = y[5]
    
    pr_h_priv = y[0]
    pr_i_priv = y[1]
    pr_d_priv = y[2]
    
    # h = h_pub + h_priv
    # i = i_pub + i_priv
    # d = d_pub + d_priv

    pr_h = pr_h_pub + pr_h_priv
    pr_i = pr_i_pub + pr_i_priv
    pr_d = pr_d_pub + pr_d_priv
    
    ###! Subsidy and treatment probs for simulation
        ###! 1: no private treatment or public treatment
        ###! 4: no private treatment but optimal public treatment
        ###! 2: no private subsidies and no public treatment
        ###! 5: no private subsidies but optimal public treatment
        ###! 3: optimal private subsidies but no public treatment
        ###! 6: optimal private subsidies and optimal public treatment
        ###! 7: treat every tree assessed as healthy or infested
        
    simType = P['simType']
    ###! cost of treatment
    c = P['c']
    ###! VALUE FO SAVING A TREE TO OWNER (uniformly distributed from a to b) AND MUNI
    a = P['a']
    b = P['b']
    Delta_m = P['Delta_m']
    Delta_m_alt = P['Delta_m_alt']
    
    ###! assessment params: pr_Hh is the prob of assessing a tree as 'H' healthy given that it is truly healthy 'h'
    pr_Hh = P['pr_Hh']
    pr_Ih = P['pr_Ih']
    pr_Dh = P['pr_Dh']

    pr_Hi = P['pr_Hi']
    pr_Ii = P['pr_Ii']
    pr_Di = P['pr_Di']

    pr_Hd = P['pr_Hd']
    pr_Id = P['pr_Id']
    pr_Dd = P['pr_Dd']


    t_horiz = P['t_horiz']
    eps_h = P['eps_h']
    eps_i = P['eps_i']
    spilloverCare = P['spill_weight']
    
    muh = m_uh(t_horiz,pr_i,beta,gamma)
    mth = m_th(t_horiz,pr_i,beta,gamma,eps_h,eps_i)
    mui = m_ui(t_horiz,gamma)
    mti = m_ti(t_horiz,gamma,eps_i)
    mud = 1
    mtd = 1
    
    t_h = muh - mth 
    t_i = mui - mti
    t_d = mud - mtd
    
    smuh = sm_uh(t_horiz,pr_i,pr_h,beta,gamma)
    smth = sm_th(t_horiz,pr_i,pr_h,beta,gamma,alpha,eps_h,eps_i)
    smui = sm_ui(pr_h,beta,gamma)
    smti = sm_ti(pr_h,beta,gamma,alpha,eps_i)
    smud = 0
    smtd = 0
    
    s_h = (smuh - smth) * spilloverCare
    s_i = (smui - smti) * spilloverCare
    s_d = (smud - smtd) * spilloverCare
    

    
    pr_hH = prob_xH('h',pr_h,pr_i,pr_d,pr_Hh,pr_Hi,pr_Hd)
    pr_iH = prob_xH('i',pr_h,pr_i,pr_d,pr_Hh,pr_Hi,pr_Hd)
    pr_dH = prob_xH('d',pr_h,pr_i,pr_d,pr_Hh,pr_Hi,pr_Hd)
    
    pr_hI = prob_xI('h',pr_h,pr_i,pr_d,pr_Ih,pr_Ii,pr_Id)
    pr_iI = prob_xI('i',pr_h,pr_i,pr_d,pr_Ih,pr_Ii,pr_Id)
    pr_dI = prob_xI('d',pr_h,pr_i,pr_d,pr_Ih,pr_Ii,pr_Id)
    
    pr_hD = prob_xD('h',pr_h,pr_i,pr_d,pr_Dh,pr_Di,pr_Dd)
    pr_iD = prob_xD('i',pr_h,pr_i,pr_d,pr_Dh,pr_Di,pr_Dd)
    pr_dD = prob_xD('d',pr_h,pr_i,pr_d,pr_Dh,pr_Di,pr_Dd)
    
    # Increase in PROBABILITY OF TREE survival given TREATMENT and IT'S ASSESSED STATE
    k_H = k_X(t_h,t_i,t_d,pr_hH,pr_iH,pr_dH)
    k_I = k_X(t_h,t_i,t_d,pr_hI,pr_iI,pr_dI)
    k_D = k_X(t_h,t_i,t_d,pr_hD,pr_iD,pr_dD)
    
    L_H = L_X(s_h,s_i,s_d,pr_hH,pr_iH,pr_dH)
    L_I = L_X(s_h,s_i,s_d,pr_hI,pr_iI,pr_dI)
    L_D = L_X(s_h,s_i,s_d,pr_hD,pr_iD,pr_dD)
    
    if simType == 1:
        pr_th_priv = 0
        pr_ti_priv = 0
        pr_th_pub = 0
        pr_ti_pub = 0
    elif simType == 4:
        pr_th_priv = 0
        pr_ti_priv = 0
        pr_th_pub = pr_Hh*pr_treated_pub(c,k_H,L_H,Delta_m_alt,a,b) +pr_Ih*pr_treated_pub(c,k_I,L_I,Delta_m_alt,a,b) + pr_Dh*pr_treated_pub(c,k_D,L_D,Delta_m_alt,a,b)
        pr_ti_pub = pr_Hi*pr_treated_pub(c,k_H,L_H,Delta_m_alt,a,b) +pr_Ii*pr_treated_pub(c,k_I,L_I,Delta_m_alt,a,b) + pr_Di*pr_treated_pub(c,k_D,L_D,Delta_m_alt,a,b)
    elif simType == 2:
        pr_th_priv = pr_Hh*pr_treated_nosub(c,k_H,a,b) +pr_Ih*pr_treated_nosub(c,k_I,a,b) + pr_Dh*pr_treated_nosub(c,k_D,a,b)
        pr_ti_priv = pr_Hi*pr_treated_nosub(c,k_H,a,b) +pr_Ii*pr_treated_nosub(c,k_I,a,b) + pr_Di*pr_treated_nosub(c,k_D,a,b)
        pr_th_pub = 0
        pr_ti_pub = 0
    elif simType == 5:
        pr_th_priv = pr_Hh*pr_treated_nosub(c,k_H,a,b) +pr_Ih*pr_treated_nosub(c,k_I,a,b) + pr_Dh*pr_treated_nosub(c,k_D,a,b)
        pr_ti_priv = pr_Hi*pr_treated_nosub(c,k_H,a,b) +pr_Ii*pr_treated_nosub(c,k_I,a,b) + pr_Di*pr_treated_nosub(c,k_D,a,b)
        pr_th_pub = pr_Hh*pr_treated_pub(c,k_H,L_H,Delta_m_alt,a,b) +pr_Ih*pr_treated_pub(c,k_I,L_I,Delta_m_alt,a,b) + pr_Dh*pr_treated_pub(c,k_D,L_D,Delta_m_alt,a,b)
        pr_ti_pub = pr_Hi*pr_treated_pub(c,k_H,L_H,Delta_m_alt,a,b) +pr_Ii*pr_treated_pub(c,k_I,L_I,Delta_m_alt,a,b) + pr_Di*pr_treated_pub(c,k_D,L_D,Delta_m_alt,a,b)
    elif simType == 3:
        pr_th_priv = pr_Hh*pr_treated(c,k_H,L_H,Delta_m,a,b) +pr_Ih*pr_treated(c,k_I,L_I,Delta_m,a,b) + pr_Dh*pr_treated(c,k_D,L_D,Delta_m,a,b)
        pr_ti_priv = pr_Hi*pr_treated(c,k_H,L_H,Delta_m,a,b) +pr_Ii*pr_treated(c,k_I,L_I,Delta_m,a,b) + pr_Di*pr_treated(c,k_D,L_D,Delta_m,a,b)
        pr_th_pub = 0
        pr_ti_pub = 0
    elif simType == 6:
        pr_th_priv = pr_Hh*pr_treated(c,k_H,L_H,Delta_m,a,b) +pr_Ih*pr_treated(c,k_I,L_I,Delta_m,a,b) + pr_Dh*pr_treated(c,k_D,L_D,Delta_m,a,b)
        pr_ti_priv = pr_Hi*pr_treated(c,k_H,L_H,Delta_m,a,b) +pr_Ii*pr_treated(c,k_I,L_I,Delta_m,a,b) + pr_Di*pr_treated(c,k_D,L_D,Delta_m,a,b)
        pr_th_pub = pr_Hh*pr_treated_pub(c,k_H,L_H,Delta_m_alt,a,b) +pr_Ih*pr_treated_pub(c,k_I,L_I,Delta_m_alt,a,b) + pr_Dh*pr_treated_pub(c,k_D,L_D,Delta_m_alt,a,b)
        pr_ti_pub = pr_Hi*pr_treated_pub(c,k_H,L_H,Delta_m_alt,a,b) +pr_Ii*pr_treated_pub(c,k_I,L_I,Delta_m_alt,a,b) + pr_Di*pr_treated_pub(c,k_D,L_D,Delta_m_alt,a,b)
    else:
        pr_th_priv = pr_Hh + pr_Ih 
        pr_ti_priv = pr_Hi + pr_Ii
        pr_th_pub = pr_Hh + pr_Ih 
        pr_ti_pub = pr_Hi + pr_Ii

    

    
    infestRate_priv = beta*(1-eps_h*pr_th_priv)*h_priv*(i_priv+i_pub)
    deathRate_priv = gamma*(1-eps_i*pr_ti_priv)*i_priv
    recoveryRate_priv = alpha*eps_i*pr_ti_priv*i_priv
   
    

    
    infestRate_pub = beta*(1-eps_h*pr_th_pub)*h_pub*(i_priv+i_pub)
    deathRate_pub = gamma*(1-eps_i*pr_ti_pub)*i_pub
    recoveryRate_pub = alpha*eps_i*pr_ti_pub*i_pub
    
    dhPrivdt =  - infestRate_priv + recoveryRate_priv
    diPrivdt = infestRate_priv - recoveryRate_priv - deathRate_priv
    ddPrivdt = deathRate_priv
    dhPubdt =  - infestRate_pub + recoveryRate_pub
    diPubdt = infestRate_pub - recoveryRate_pub - deathRate_pub
    ddPubdt = deathRate_pub
    
    return (dhPrivdt,diPrivdt,ddPrivdt,dhPubdt,diPubdt,ddPubdt)


