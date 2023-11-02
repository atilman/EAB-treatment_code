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
"This calculates the likelihood of an underlying true health state, given assessed as healthy 'H' "
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

"This calculates the likelihood of an underlying true health state, given assessed as infested 'I' "
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

"This calculates the likelihood of an underlying true health state, given assessed as dying / dead 'D' "
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
def s_opt(c,k_Phi,Delta_m,a,b): #This is for ecosystem service optimization
    if c <= a*k_Phi:
        s = 0
    elif c + b*k_Phi - 2*a*k_Phi <= Delta_m*k_Phi:
        s = c - a*k_Phi
    elif np.abs(c - b*k_Phi) < Delta_m*k_Phi:
        s = (1/2)*(Delta_m*k_Phi + c - b*k_Phi)
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

"TREATMET PROBABILITY UNDER OPTIMAL SUBSIDY"
def pr_treated(c,k_Phi,Delta_m,a,b):
    if c <= a*k_Phi:
        pr_T = 1
    elif c + b*k_Phi - 2*a*k_Phi <= Delta_m*k_Phi:
        pr_T = 1
    elif np.abs(c - b*k_Phi) < Delta_m*k_Phi:
        pr_T = (Delta_m*k_Phi - c + b*k_Phi)/(2*k_Phi*(b-a))
    elif c < b*k_Phi and Delta_m*k_Phi <= b*k_Phi - c:
        pr_T = (b*k_Phi - c)/(k_Phi*(b-a))
    elif b*k_Phi <= c and Delta_m*k_Phi <= c - b*k_Phi:
        pr_T = 0
    else:
        # print('theres an error somewhere!')
        pr_T = 0
    return pr_T


"TREATMENT OF PUBLICLY OWNED TREES"
def pr_treated_pub(c,k_Phi,Delta_m,a,b):
    if c <= Delta_m*k_Phi:
        pr_T = 1
    else:
        pr_T = 0
    return pr_T


"SURVIVAL PROBABILITY OF A TREE ASSESSED IN STATE X GIVEN OPTIMAL POLICY FOLLOWED"
def pr_survX(pr_treat,pr_hX,pr_iX,pr_dX,hth,huh,hti,hui,htd,hud):
    prSurv = pr_treat*(hth*pr_hX + hti*pr_iX + htd*pr_dX) + (1-pr_treat)*(huh*pr_hX + hui*pr_iX + hud*pr_dX)
    return prSurv

"MUNICIPAL PAYOFF FUNCTION UNDER TREATMENT, X=t, or non-treatment, X=u"
def Pi(V_m,W_m,hXh,hXi,pr_hPhi,pr_iPhi,pr_dPhi):
    Pi = pr_hPhi*(hXh*V_m-(1-hXh)*W_m) + pr_iPhi*(hXi*V_m-(1-hXi)*W_m)-pr_dPhi*W_m
    return Pi

"EXPECTED MUNI UTILITY UNDER OPTIMAL POLICY"
def EU_m(c,k_Phi,Delta_m,a,b,PiU):
    if c <= a*k_Phi:
        EUm = PiU + Delta_m*k_Phi
    elif c + b*k_Phi - 2*a*k_Phi <= Delta_m*k_Phi:
        EUm = PiU + Delta_m*k_Phi - c + a*k_Phi
    elif np.abs(c - b*k_Phi) < Delta_m*k_Phi:
        EUm = PiU + (Delta_m*k_Phi - c + b*k_Phi)**2/(4*k_Phi*(b-a))
    elif c < b*k_Phi and Delta_m*k_Phi <= b*k_Phi - c:
        EUm = PiU + Delta_m*k_Phi * (b*k_Phi - c)/(k_Phi*(b-a))
    elif b*k_Phi <= c and Delta_m*k_Phi <= c - b*k_Phi:
        EUm = PiU
    else:
        print('theres an error somewhere!')
    return EUm


"PROBABILITY OF SAVING A TREEE WITH TREATMENT GIVEN IT'S ASSESSED STATE IS X"
def k_X(t_h,t_i,t_d,pr_hX,pr_iX,pr_dX):
    treatmentImpact = pr_hX * t_h + pr_iX * t_i + pr_dX * t_d 
    return treatmentImpact

"""
FOR DYNAMICAL SIMULATIONS:
    DISEASE DYNAMICS EQUATIONS
"""
"DYNAMICAL EQUATIONS UNDER NO TREATMENT"
def hid(t,y,beta,gamma,alpha,P):
    h=y[0]
    i=y[1]
    #d=y[2]
    dhdt = -beta*h*i
    didt = beta*h*i - gamma*i
    dddt = gamma*i
    return (dhdt,didt,dddt)



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
        ###! 2: no private treatment but optimal public treatment
        ###! 3: no private subsidies and no public treatment
        ###! 4: no private subsidies but optimal public treatment
        ###! 5: optimal private subsidies but no public treatment
        ###! 6: optimal private subsidies and optimal public treatment
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

    ###! outcomes of treatment
    hth = P['hth'] # prob of a tree staying healty 'h', given that it is treated 't' and healthy 'h'
    huh = P['huh']  # prob of a tree staying healty 'h', given that it is untreated 'u' and healthy 'h'
    hti = P['hti'] # prob of a tree becoming healty 'h', given that it is treated 't' and infested 'i'
    hui = P['hui'] # prob of a tree becoming healty 'h', given that it is untreated 'u' and infested 'i'
    htd = P['htd']   # prob of saving a dead tree -> 0
    hud = P['hud']   # prob of saving a dead tree -> 0
    
    t_h = hth - huh
    t_i = hti - hui
    t_d = htd - hud
    
    pr_hH = prob_xH('h',pr_h,pr_i,pr_d,pr_Hh,pr_Hi,pr_Hd)
    pr_iH = prob_xH('i',pr_h,pr_i,pr_d,pr_Hh,pr_Hi,pr_Hd)
    pr_dH = prob_xH('d',pr_h,pr_i,pr_d,pr_Hh,pr_Hi,pr_Hd)
    
    pr_hI = prob_xI('h',pr_h,pr_i,pr_d,pr_Ih,pr_Ii,pr_Id)
    pr_iI = prob_xI('i',pr_h,pr_i,pr_d,pr_Ih,pr_Ii,pr_Id)
    pr_dI = prob_xI('d',pr_h,pr_i,pr_d,pr_Ih,pr_Ii,pr_Id)
    
    pr_hD = prob_xD('h',pr_h,pr_i,pr_d,pr_Dh,pr_Di,pr_Dd)
    pr_iD = prob_xD('i',pr_h,pr_i,pr_d,pr_Dh,pr_Di,pr_Dd)
    pr_dD = prob_xD('d',pr_h,pr_i,pr_d,pr_Dh,pr_Di,pr_Dd)
    
    # PROBABILITY OF SAVING A TREEE WITH TREATMENT GIVEN IT'S ASSESSED STATE IS X
    k_H = k_X(t_h,t_i,t_d,pr_hH,pr_iH,pr_dH)
    k_I = k_X(t_h,t_i,t_d,pr_hI,pr_iI,pr_dI)
    k_D = k_X(t_h,t_i,t_d,pr_hD,pr_iD,pr_dD)
    if simType == 1:
        pr_th_priv = 0
        pr_ti_priv = 0
        pr_th_pub = 0
        pr_ti_pub = 0
    elif simType == 2:
        pr_th_priv = 0
        pr_ti_priv = 0
        pr_th_pub = pr_Hh*pr_treated_pub(c,k_H,Delta_m_alt,a,b) +pr_Ih*pr_treated_pub(c,k_I,Delta_m_alt,a,b) + pr_Dh*pr_treated_pub(c,k_D,Delta_m_alt,a,b)
        pr_ti_pub = pr_Hi*pr_treated_pub(c,k_H,Delta_m_alt,a,b) +pr_Ii*pr_treated_pub(c,k_I,Delta_m_alt,a,b) + pr_Di*pr_treated_pub(c,k_D,Delta_m_alt,a,b)
    elif simType == 3:
        pr_th_priv = pr_Hh*pr_treated_nosub(c,k_H,a,b) +pr_Ih*pr_treated_nosub(c,k_I,a,b) + pr_Dh*pr_treated_nosub(c,k_D,a,b)
        pr_ti_priv = pr_Hi*pr_treated_nosub(c,k_H,a,b) +pr_Ii*pr_treated_nosub(c,k_I,a,b) + pr_Di*pr_treated_nosub(c,k_D,a,b)
        pr_th_pub = 0
        pr_ti_pub = 0
    elif simType == 4:
        pr_th_priv = pr_Hh*pr_treated_nosub(c,k_H,a,b) +pr_Ih*pr_treated_nosub(c,k_I,a,b) + pr_Dh*pr_treated_nosub(c,k_D,a,b)
        pr_ti_priv = pr_Hi*pr_treated_nosub(c,k_H,a,b) +pr_Ii*pr_treated_nosub(c,k_I,a,b) + pr_Di*pr_treated_nosub(c,k_D,a,b)
        pr_th_pub = pr_Hh*pr_treated_pub(c,k_H,Delta_m_alt,a,b) +pr_Ih*pr_treated_pub(c,k_I,Delta_m_alt,a,b) + pr_Dh*pr_treated_pub(c,k_D,Delta_m_alt,a,b)
        pr_ti_pub = pr_Hi*pr_treated_pub(c,k_H,Delta_m_alt,a,b) +pr_Ii*pr_treated_pub(c,k_I,Delta_m_alt,a,b) + pr_Di*pr_treated_pub(c,k_D,Delta_m_alt,a,b)
    elif simType == 5:
        pr_th_priv = pr_Hh*pr_treated(c,k_H,Delta_m,a,b) +pr_Ih*pr_treated(c,k_I,Delta_m,a,b) + pr_Dh*pr_treated(c,k_D,Delta_m,a,b)
        pr_ti_priv = pr_Hi*pr_treated(c,k_H,Delta_m,a,b) +pr_Ii*pr_treated(c,k_I,Delta_m,a,b) + pr_Di*pr_treated(c,k_D,Delta_m,a,b)
        pr_th_pub = 0
        pr_ti_pub = 0
    else:
        pr_th_priv = pr_Hh*pr_treated(c,k_H,Delta_m,a,b) +pr_Ih*pr_treated(c,k_I,Delta_m,a,b) + pr_Dh*pr_treated(c,k_D,Delta_m,a,b)
        pr_ti_priv = pr_Hi*pr_treated(c,k_H,Delta_m,a,b) +pr_Ii*pr_treated(c,k_I,Delta_m,a,b) + pr_Di*pr_treated(c,k_D,Delta_m,a,b)
        pr_th_pub = pr_Hh*pr_treated_pub(c,k_H,Delta_m_alt,a,b) +pr_Ih*pr_treated_pub(c,k_I,Delta_m_alt,a,b) + pr_Dh*pr_treated_pub(c,k_D,Delta_m_alt,a,b)
        pr_ti_pub = pr_Hi*pr_treated_pub(c,k_H,Delta_m_alt,a,b) +pr_Ii*pr_treated_pub(c,k_I,Delta_m_alt,a,b) + pr_Di*pr_treated_pub(c,k_D,Delta_m_alt,a,b)
    

    
    infestRate_priv = beta*(1-hth*pr_th_priv)*h_priv*(i_priv+i_pub)
    deathRate_priv = gamma*(1-hui-t_i*pr_ti_priv)*i_priv
    recoveryRate_priv = alpha*(hui+t_i*pr_ti_priv)*i_priv

    
    infestRate_pub = beta*(1-hth*pr_th_pub)*h_pub*(i_priv+i_pub)
    deathRate_pub = gamma*(1-hui-t_i*pr_ti_pub)*i_pub
    recoveryRate_pub = alpha*(hui+t_i*pr_ti_pub)*i_pub
    
    dhPrivdt =  - infestRate_priv + recoveryRate_priv
    diPrivdt = infestRate_priv - recoveryRate_priv - deathRate_priv
    ddPrivdt = deathRate_priv
    dhPubdt =  - infestRate_pub + recoveryRate_pub
    diPubdt = infestRate_pub - recoveryRate_pub - deathRate_pub
    ddPubdt = deathRate_pub
    
    return (dhPrivdt,diPrivdt,ddPrivdt,dhPubdt,diPubdt,ddPubdt)


"RETURNS TREATMENT PROBS FOR TREES WITH A TRULY HEALTHY OR INFESTED UNDERLYING STATE"
def pr_treat_hi(y,P):
    pr_h=y[0]
    pr_i=y[1]
    pr_d=y[2]
    
    ###! cost of treatment
    c = P['c']
    ###! VALUE FO SAVING A TREE TO OWNER AND MUNI
    a = P['a']
    b = P['b']
    Delta_m = P['Delta_m']

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

    ###! outcomes of treatment
    hth = P['hth'] # prob of a tree staying healty 'h', given that it is treated 't' and healthy 'h'
    huh = P['huh']  # prob of a tree staying healty 'h', given that it is untreated 'u' and healthy 'h'
    hti = P['hti'] # prob of a tree becoming healty 'h', given that it is treated 't' and infested 'i'
    hui = P['hui'] # prob of a tree becoming healty 'h', given that it is untreated 'u' and infested 'i'
    htd = P['htd']   # prob of saving a dead tree -> 0
    hud = P['hud']   # prob of saving a dead tree -> 0
    
    t_h = hth - huh
    t_i = hti - hui
    t_d = htd - hud
    
    pr_hH = prob_xH('h',pr_h,pr_i,pr_d,pr_Hh,pr_Hi,pr_Hd)
    pr_iH = prob_xH('i',pr_h,pr_i,pr_d,pr_Hh,pr_Hi,pr_Hd)
    pr_dH = prob_xH('d',pr_h,pr_i,pr_d,pr_Hh,pr_Hi,pr_Hd)
    
    pr_hI = prob_xI('h',pr_h,pr_i,pr_d,pr_Ih,pr_Ii,pr_Id)
    pr_iI = prob_xI('i',pr_h,pr_i,pr_d,pr_Ih,pr_Ii,pr_Id)
    pr_dI = prob_xI('d',pr_h,pr_i,pr_d,pr_Ih,pr_Ii,pr_Id)
    
    pr_hD = prob_xD('h',pr_h,pr_i,pr_d,pr_Dh,pr_Di,pr_Dd)
    pr_iD = prob_xD('i',pr_h,pr_i,pr_d,pr_Dh,pr_Di,pr_Dd)
    pr_dD = prob_xD('d',pr_h,pr_i,pr_d,pr_Dh,pr_Di,pr_Dd)
    
    k_H = k_X(t_h,t_i,t_d,pr_hH,pr_iH,pr_dH)
    k_I = k_X(t_h,t_i,t_d,pr_hI,pr_iI,pr_dI)
    k_D = k_X(t_h,t_i,t_d,pr_hD,pr_iD,pr_dD)
    
    pr_th = pr_Hh*pr_treated(c,k_H,Delta_m,a,b) +pr_Ih*pr_treated(c,k_I,Delta_m,a,b) + pr_Dh*pr_treated(c,k_D,Delta_m,a,b)
    pr_ti = pr_Hi*pr_treated(c,k_H,Delta_m,a,b) +pr_Ii*pr_treated(c,k_I,Delta_m,a,b) + pr_Di*pr_treated(c,k_D,Delta_m,a,b)
    
    return pr_th , pr_ti


