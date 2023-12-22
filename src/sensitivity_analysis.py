"""This script contains functions that performs sensitivity analysis of different parameters in the quorum sensing system 
    and returns a dataframe containing the range of parameter values and dynamics metrics
"""

import pandas as pd
from .solve_ode import solve_qcrispri, solve_qcrispri_sponge, alpha0, kd1, kd2, t
from .calculate_metrics import CalculateDynamicsMetrics


def single_param_vary(param, param_range):
    """solve qcrispri model and calculate switch density, foldrepression, transition time for a range of values of any parameter: 
    alphaT, activity_ratio, kDNA_bind, kcomplex.
    
    INPUTS:
    param: name of the parameter to vary
    param_range: list of values to use for the mentioned parameter <param>
    """
    metric_df = pd.DataFrame(columns=[f"{param}_value", 'switching_density', 'fold_repression', 'transition_time'])
    gfp_df = pd.DataFrame(columns = [f"{param}_value",'control', 'GFP', 'N' ])
    temp = pd.DataFrame()
    params = {'alphaT': alpha0,
              'activity_ratio': 0.1,
              'Kd1': kd1,
              'Kd2': kd2}
    for i in param_range:
        params[param] = i
        qcrispri_df = solve_qcrispri(**params)
        
        temp[f"{param}_value"] = [i]*len(t)
        temp['control'] = qcrispri_df.GFPc
        temp['GFP'] = qcrispri_df.GFP
        temp['N'] = qcrispri_df.N
        gfp_df = pd.concat([gfp_df, temp])

        switch_density, fold_repression, transition_time = CalculateDynamicsMetrics(qcrispri_df.GFP, t, qcrispri_df.N).get_dynamic_metrics()

        new_row = {f"{param}_value": i, 'switching_density': switch_density, 'fold_repression': fold_repression, 'transition_time':transition_time}
        metric_df = metric_df.append(new_row, ignore_index = True)
    return(gfp_df, metric_df)


# solve qcrispri model while varying 2 parameters and calcuate dynamics parameters and compile gfp data into a df
def twoParam_vary(param1, param1_range, param2, param2_range):
    param_df = pd.DataFrame(columns=[param1, param2, 'act_density', 'FoldRepression', 'Transition_time', 'MaxGFP'])
    gfp_df = pd.DataFrame(columns = [param1, param2, 'control', 'GFP', 'N' ])
    temp = pd.DataFrame()
    params = {'alphaT': alpha0,
              'activity_ratio': 0.1,
              'Kd1': kd1,
              'Kd2': kd2}
    qcrispridf = convert2df()
    for i in param1_range:
        params[param1] = i
        for j in param2_range:
            params[param2] = j
            qcrispridf.qcrispri(solve_qCRISPRi(**params))
            qcrispri_df = qcrispridf.returndf()
            #calculate dynamic parameters
            s_density = switch_density(qcrispri_df, t)
            fold_rep = foldRepression(qcrispri_df)
            t_time = transition_time(qcrispri_df, t)
            max_gfp = max(qcrispri_df.GFP)
            new_row = {param1:i, param2:j, 'act_density': s_density, 'FoldRepression': fold_rep, 'Transition_time':t_time, 'MaxGFP':max_gfp}
            param_df = param_df.append(new_row, ignore_index = True)
            
            #create gfp dataframe
            temp[param1] = [i]*len(t)
            temp[param2] = [j]*len(t)
            temp['control'] = qcrispri_df.GFPc
            temp['GFP'] = qcrispri_df.GFP
            temp['N'] = qcrispri_df.N
            gfp_df = pd.concat([gfp_df, temp])
    return(param_df, gfp_df)

#vary decoy sites and calculate parameters and create a df with gfp and related values
def decoysite_vary(param, param_range):
    param_df = pd.DataFrame(columns=[param, 'act_density', 'FoldRepression', 'Transition_time', 'MaxGFP'])
    gfp_df = pd.DataFrame(columns = [param,'control', 'GFP', 'N' ])
    temp = pd.DataFrame()
    params = {'alphaT': alpha0,
              'activity_ratio': 0.1,
              'Kd1': kd1,
              'Kd2': kd2,
              'decoy_sites': 10}
    qcrispridf = convert2df()
    for i in param_range:
        params[param] = i
        qcrispridf.qcrispri_sponge(solve_qcrispri_sponge(**params))
        qcrispri_df = qcrispridf.returndf()

        temp[param] = [i]*len(t)
        temp['control'] = qcrispri_df.GFPc
        temp['GFP'] = qcrispri_df.GFP
        temp['N'] = qcrispri_df.N
        gfp_df = pd.concat([gfp_df, temp])

        s_density = switch_density(qcrispri_df, t)
        fold_rep = foldRepression(qcrispri_df)
        t_time = transition_time(qcrispri_df, t)
        new_row = {param:i, 'act_density': s_density, 'FoldRepression': fold_rep, 'Transition_time':t_time, 'MaxGFP': max(qcrispri_df.GFP)}
        param_df = param_df.append(new_row, ignore_index = True)
    return(gfp_df, param_df)

#gets data required to make bar plots for 2 param vary
def newRow_forBarplot_df(params, alphaT, param2, param2_value):
    idx = params[(params['alphaT'] == alphaT) & (params[param2] == param2_value)].index[0]
    newrow = {'alphaT': alphaT, param2:param2_value, 'act_density':params['act_density'].iloc[idx], 
              'Transition_time':params['Transition_time'].iloc[idx], 'MaxGFP':params['MaxGFP'].iloc[idx], 'FoldRepression':params['FoldRepression'].iloc[idx]}
    return(newrow)

#returns the default parameter value of some parameters

def get_default_param_value(param_name):
    default_params = {'alphaT': alpha0,
                      'activity_ratio': 0.1,
                      'Kd1': kd1,
                      'Kd2': kd2,
                      'Stringency': 10}
    return(default_params[param_name])

#gets the index of a default parameter in a dataframe

def get_default_param_index(param_range, param_name):
    key = get_default_param_value(param_name)
    return(param_range[param_range == key].index[0])
    