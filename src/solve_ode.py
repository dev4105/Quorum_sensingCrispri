"""This script contains functions that solves different models based on quorum sensing and returns the solution array. """

from configparser import ConfigParser
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from .generate_ode import ode_gen_quorum, ode_gen_qcrispri, ode_gen_qcrispri_simple, ode_gen_qcrispri_sponge

#load constants form config file
config_object = ConfigParser()
config_object.read("./src/config.ini")

productionRates = config_object["PRODUCTION_RATES"]
valency = config_object["VALENCY"]
rxnRates = config_object['REACTION_RATES']
degRates = config_object['DEGRADATION_RATES']
misc = config_object['MISCELLANEOUS']
stringent_sys_params = config_object['STRINGENT_PARAMETERS']

alpha0 = float(productionRates["alpha0"])    #basal expression of luxI
alphaI = float(productionRates["alphaI"])     #activated expression of luxI
alphaR = float(productionRates["alphaR"])      #constitutive expresssion of luxR promoter
alphaG = float(productionRates["alphaG"])     #constitutive expression of GFP
alphaIr = float(productionRates['alphaIr'])
alphaIrt = float(productionRates['alphaIrt'])

nA = float(valency["nA"])   #valency of AHL
nG = float(valency["nG"])   #valency of GFP
n = float(valency['n'])     #hill coefficient for luxR-AHL
ncas = float(valency['ncas']) #hill coefficient for dCas9-promoter DNA

k1R = float(rxnRates["k1R"])    #hr-1
k2R = float(rxnRates["k2R"])       #hr-1
kd1 = float(rxnRates["kd1"])       #nM LuxR_AHL complex dissociation constant
kd2 = float(rxnRates["kd2"])   #nM LuxR_AHL DNA binding constant
kcasR = float(rxnRates["kcasR"])     #hr-1
kdCas = float(rxnRates["kdCas"])   #nM

dR = float(degRates["dR"])
dA = float(degRates["dA"])
dCas = float(degRates["dCas"])
dG = float(degRates["dG"])
dc = float(degRates["dc"])

r = float(misc["r"]) #hr-1
K = eval(misc["K"])  #CFU/mL
pT = float(misc["pT"]) #nM
pTsponge = float(misc["pTsponge"])  #nM
Vc = eval(misc["Vc"])  #mL Volume of one cell
mu = float(misc["mu"]) #hr-1
ec50 = float(misc["ec50"])
m = float(misc["m"]) 
c = float(misc["c"])


# parameters for stringent system
alpha_pl = float(stringent_sys_params["alpha_pl"])    #basal expression of luxI
alpha_rl = float(stringent_sys_params["alpha_rl"])     #basal expression due to leaky LuxR activation
alpha_ra = float(stringent_sys_params["alpha_ra"])      #activated expression of pluxI
alpha_g = float(stringent_sys_params["alpha_g"])     #constitutive expression of GFP
nHill_A = float(stringent_sys_params["nHill_A"])    #hill coeff for AHL activation of pLuxI
nHill_c = float(stringent_sys_params["nHill_c"])     #hill coeff for dCas9 binding to dna
nAHL = float(stringent_sys_params["nAHL"]) 
Kd = float(stringent_sys_params["Kd"])      #saturation coeff for AHL activation of pLuxI
KdCas = float(stringent_sys_params["KdCas"])     #saturation coeff for dCas9 binding to dna
Kdeg_A = float(stringent_sys_params["Kdeg_A"])    #deg rate for AHL
Kdeg_Cas = float(stringent_sys_params["Kdeg_Cas"])     #deg rate for dCas9
Kdeg_G = float(stringent_sys_params["Kdeg_G"])      #deg rate for GFP
m1 = float(stringent_sys_params["m1"])     #
c1 = float(stringent_sys_params["c1"])     #
r1 = float(stringent_sys_params["r1"])     #
K1 = eval(stringent_sys_params["K1"])     #


t = np.linspace(0,48, 3000)


def convert2df(sol_array:list[list], column_list: list):
    sol_df = pd.DataFrame(columns = column_list)
    for i, column in enumerate(sol_df.columns):
        sol_df[column] = sol_array[:,i]
    return sol_df
        

def solve_quorum():
    alpha = [alpha0, alphaI, alphaR]
    n = [nA, nG]
    kd = [k1R, k2R, kd1, kd2]
    deg = [dR, dA, dG, dc]
    miscc = [r, K, Vc, pT, mu, m, c]
    y0 = [pT, 0, 10**7, 0, 0, 0, 0 ]
    # Solve ODEs for "BASE MODEL

    sol = odeint(ode_gen_quorum, y0, t,args=(alpha, n, kd, deg, miscc))
    cols = ['inAcP', 'AcP', 'N', 'LuxR', 'AHL', 'LuxR_A', 'GFP']
    return(convert2df(sol, cols))


def solve_qcrispri(alphaT = alpha0, alphaIac = alphaI,  activity_ratio=0.1, Kd1 = kd1, Kd2 = kd2, gfp0 = 1000):

    #use 50% as default to study + and - variation
    alphaG_deac = alphaG*activity_ratio
    kd1 = Kd1
    kd2 = Kd2
    y0=[pT, 0, 10**7, 0, 0, 0, 1000, gfp0, 0, pT, 0]
    alpha = [alpha0, alphaIac, alphaR, alphaG, alphaG_deac, alphaT]
    n = nA
    kd = [k1R, k2R,kcasR, kd1, kd2, kdCas]
    deg = [dR, dA, dG, dCas, dc]
    miscc = [r, K, Vc, pT, mu, m, c]

    # Solve ODEs for "qCRISPRi" 
    sol = odeint(ode_gen_qcrispri, y0, t, args=(alpha, n, kd, deg, miscc))
    cols = ['inAcP', 'AcP', 'N', 'LuxR', 'AHL', 'LuxR_A', 'GFPc', 'GFP', 'dCas9', 'pC', 'pD']
    return(convert2df(sol, cols))



def solve_qcrispri_simple(alphaRL = alpha_rl, alphaRA = alpha_ra, nH_AHL = nHill_A, KD = Kd, rs = r1):
    
    y0=[0, 0, 1000, 10**7]
    alpha = [alpha_pl, alphaRL, alphaRA, alpha_g]
    n = [nH_AHL, nHill_c, nAHL]
    kd = [KD, KdCas]
    deg = [Kdeg_A, Kdeg_Cas, Kdeg_G]
    miscc = [rs, K1, Vc, m1, c1]

    # Solve ODEs for "qCRISPRi" 
    sol = odeint(ode_gen_qcrispri_simple, y0, t, args=(alpha, n, kd, deg, miscc))
    cols = ['AHL', 'dCas9', 'GFP', 'N']
    return(convert2df(sol, cols))



def solve_qcrispri_sponge(alphaT = alpha0, activity_ratio = 0.1, Kd1 = kd1, Kd2 = kd2, decoy_sites = 10):
    
    alphaG_deac = alphaG*activity_ratio
    kd1 = Kd1
    kd2 = Kd2
    pTsponge = pT*decoy_sites
  
    y0=[pT, 0, 10**7, 0, 0, 0, 1000, 1000, 0, pT, 0, pTsponge, 0]
    alpha = [alpha0, alphaI, alphaR, alphaG, alphaG_deac, alphaT]
    n = nA
    kd = [k1R, k2R,kcasR, kd1, kd2, kdCas]
    deg = [dR, dA, dG, dCas, dc]
    miscc = [r, K, Vc, pT, mu, m, c]

    # Solve ODEs for "qCRISPRi sponge" 
    sol = odeint(ode_gen_qcrispri_sponge, y0, t,args=(alpha, n, kd, deg, miscc))
    cols = ['inAcP', 'AcP', 'N', 'LuxR', 'AHL', 'LuxR_A', 'GFPc', 'GFP', 'dCas9', 'pC', 'pD', 'pS', 'pSB']
    return(convert2df(sol, cols))