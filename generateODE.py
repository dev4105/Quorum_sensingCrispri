"""
This script generates the ODEs for different quorum sensing systems 
"""



def ode_genQuorum(y, t, alpha, n, kd, deg, misc ):
    #unpack y values. 
    inAcP, AcP, N, LuxR, A, LuxR_A, GFP = y
    
    #unpack grouped constants:
    alpha0, alphaI, alphaR = alpha   #Fractional and max transcription rates
    nA, nG = n
    k1R, k2R, kd1, kd2 = kd  #Dissociation constants and reverse rates
    dR, dA, dG, dc = deg  #Degradation rates
    r, K, Vc, pT, mu, m, c = misc # misc parameters
    
    #return values as a list
    values = [
              #1. inAcP: Inactive promoter (basal/leaky expression)
              -(k2R/kd2)*LuxR_A*inAcP + k2R*AcP,
              
              #2. AcP: Active promoter
              (k2R/kd2)*LuxR_A*inAcP - k2R*AcP, 
              
              #3. N (cell density CFU/ml)
              r*N*(1- N/K),
              
              #4. LuxR - constitutive exp of luxR - luxR converted to luxR_Ac + 
              #   deactivation of luxR_Ac - (degradation+dilution) of luxR
              alphaR*pT - 2*(k1R/kd1)*(LuxR**2)*(A**2) + 2*k1R*LuxR_A - (mu+ dR)*LuxR,
              
              #5. AHL - creation: expression from inactivated and activated LuxI promoter
              #             + dissociation of activated LuxR_A
              #   destruction: association of AHL with LuxR + degradation of AHL
              (alphaI*AcP + alpha0*inAcP)*nA*Vc*N + 2*k1R*LuxR_A - 2*(k1R/kd1)*(LuxR**2)*(A**2)  - (mu + dA)*A,
              
              #6. LuxR_A - creation: association of luxR and AHL
              #   destruction: dissociation of activated LuxR + binding to inactive luxI  + dissociation from luxI 
              #              degradation of activated luxR
              (k1R/kd1)*(LuxR**2)*(A**2) - k1R*LuxR_A - (k2R/kd2)*LuxR_A*inAcP + k2R*AcP -(mu + dc)*LuxR_A,
              
              #7. GFP - creation: expression from active and inactive luxI promoter
              (alphaI*AcP + alpha0*inAcP)*nG - (mu+ dG)*GFP 
             ]  
              
    
    return(values)


#Function to generate ODEs for the qCRISPRi system

def ode_gen_qCRISPRi(y, t, alpha, n, kd, deg, misc ):
    #unpack y values. 
    inAcP, AcP, N, LuxR, A, LuxR_A, GFPc, GFP, dCas9, pC, pD = y
    
    #unpack grouped constants:
    alpha0, alphaI, alphaR, alphaG, alphaG_deac, alphaT = alpha   #Fractional and max transcription rates
    nA = n
    k1R, k2R, kcasR, kd1, kd2, kdCas = kd  #Dissociation constants and reverse rates
    dR, dA, dG, dCas, dc = deg  #Degradation rates
    r, K, Vc, pT, mu, m, c = misc # misc parameters
    
    #return values as a list
    values = [
              #1. inAcP Inactive promoter (basal/leaky expression)
              -(k2R/kd2)*LuxR_A*inAcP + k2R*AcP, 
              #2. AcP Active promoter (LuxR:AHL activated)
              (k2R/kd2)*LuxR_A*inAcP - k2R*AcP,  
              #3. N (cell density CFU/l)
              r*N*(1- N/K),  
              #4. LuxR
              alphaR*pT - 2*(k1R/kd1)*(LuxR)*(A) + 2*k1R*LuxR_A - (mu + dR)*LuxR,   
              #5. AHL
              (alphaI*AcP + alpha0*inAcP)*nA*Vc*N + k1R*LuxR_A - (k1R/kd1)*(LuxR)*(A)  - (mu + dA)*A,  
              #6. LuxR_A   
              (k1R/kd1)*(LuxR)*(A) - k1R*LuxR_A - (k2R/kd2)*LuxR_A*inAcP + k2R*AcP - (dc+mu)*LuxR_A ,      
              #7. GFP just constitutive (control)
              m*(alphaG*pT - (mu + dG)*GFPc) + c,
              #8. GFP Cas9 repressed
              m*(alphaG*pC + alphaG_deac*pD - (mu + dG)*GFP) + c,  
              #9. dCas9
              #alphaT is a placeholder for alpha0. The leakiness is only reduced for the dCas9 expressing luxI promoter and not the promoter expressing luxI itself in order for the system
              #to remain autonomous. when varying leakiness of the system only alphaT is varied. 
              alphaI*AcP + alphaT*inAcP - ( mu+ dCas)*dCas9 -(kcasR/kdCas)*dCas9*pC + kcasR*pD,   
              #10. pC - Constitutive promoter
              -(kcasR/kdCas)*dCas9*pC + kcasR*pD,  
              #11. pD - dCas9 bound repressed promoter 
              (kcasR/kdCas)*dCas9*pC - kcasR*pD
             ]
    return(values)




# simplified qCRIPSRi system to model LuxR stringency

def ode_gen_qCRISPRiSimple(y, t, alpha, n, kd, deg, misc ):
    #unpack y values. 
    AHL, dCas9, GFP, N  = y
    
    #unpack grouped constants:
    alpha_pl, alpha_rl, alpha_ra, alpha_g = alpha   #Fractional and max transcription rates
    nHill_A, nHill_c, nAHL = n
    Kd, KdCas = kd  #Dissociation constants and reverse rates
    Kdeg_A, Kdeg_Cas, Kdeg_G = deg  #Degradation rates
    r1, K1, Vc, m1, c1 = misc # misc parameters8
    
    #return values as a list
    values = [
             (alpha_pl + alpha_rl + (alpha_ra)*(AHL**nHill_A)/(Kd**nHill_A + AHL**nHill_A))*Vc*nAHL*N - Kdeg_A*AHL, #AHL
              alpha_pl + alpha_rl + (alpha_ra)*(AHL**nHill_A)/(Kd**nHill_A + AHL**nHill_A) - Kdeg_Cas*dCas9, #dCas9
              m1*(alpha_g*0.1 + alpha_g*0.9*(KdCas**nHill_c)/(KdCas**nHill_c + dCas9**nHill_c)- Kdeg_G*GFP) + c1, #GFP
              r1*N*(1-N/K1) #N
             ] 
    return(values)



#Function to generate ODEs for the qCRISPRi sponge system

def ode_gen_qCRISPRi_sponge(y, t, alpha, n, kd, deg, misc ):
    #unpack y values. 
    inAcP, AcP, N, LuxR, A, LuxR_A, GFPc, GFP, dCas9, pC, pD, pS, pSB = y
    
    #unpack grouped constants:
    alpha0, alphaI, alphaR, alphaG, alphaG_deac, alphaT = alpha   #Fractional and max transcription rates
    nA = n
    k1R, k2R, kcasR, kd1, kd2, kdCas = kd  #Dissociation constants and reverse rates
    dR, dA, dG, dCas, dc = deg  #Degradation rates
    r, K, Vc, pT, mu, m, c = misc # misc parameters
    
    #return values as a list
    values = [
            #1. inAcP Inactive promoter (basal/leaky expression)
            -(k2R/kd2)*LuxR_A*inAcP + k2R*AcP,  
            #2. AcP Active promoter (LuxR:AHL activated)
            (k2R/kd2)*LuxR_A*inAcP - k2R*AcP,  
            #3. N (cell density CFU/l)
            r*N*(1- N/K),  
            #4. LuxR
            alphaR*pT - 2*(k1R/kd1)*(LuxR)*(A) + 2*k1R*LuxR_A - (mu + dR)*LuxR,
            #5. AHL
            (alphaI*AcP + alpha0*inAcP)*nA*Vc*N + k1R*LuxR_A - (k1R/kd1)*(LuxR)*(A)  - dA*A,  
            #6. LuxR_A      
            (k1R/kd1)*(LuxR)*(A) - k1R*LuxR_A - (k2R/kd2)*LuxR_A*inAcP + k2R*AcP - (dc+mu)*LuxR_A , 
            #7. GFP just constitutive (control)
            m*(alphaG*pT - (mu + dG)*GFPc) + c, 
            #8. GFP Cas9 repressed
            m*(alphaG*pC + alphaG_deac*pD - (mu + dG)*GFP) + c,  
            #9. dCas9
            alphaI*AcP + alphaT*inAcP - ( mu+ dCas)*dCas9 -(kcasR/kdCas)*dCas9*pC + kcasR*pD -(kcasR/kdCas)*dCas9*pS + kcasR*pSB ,   
            #10. pC - Constitutive promoter     
            -(kcasR/kdCas)*dCas9*pC + kcasR*pD, 
            #11. pD - dCas9 bound repressed promoter       
            (kcasR/kdCas)*dCas9*pC - kcasR*pD,  
            #12. pS - Free sponge site
            -(kcasR/kdCas)*dCas9*pS + kcasR*pSB,     
            #13. pSB - bound sponge site 
            (kcasR/kdCas)*dCas9*pS - kcasR*pSB
            ]      
    return(values)