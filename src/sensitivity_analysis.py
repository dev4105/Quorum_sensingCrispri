"""This script contains functions that performs sensitivity analysis of different parameters in the quorum sensing system 
    and returns a dataframe containing the range of parameter values and dynamics metrics
"""

import pandas as pd
from .solve_ode import solve_qcrispri, solve_qcrispri_sponge, alpha0, kd1, kd2, t
from .calculate_metrics import CalculateDynamicsMetrics


def vary_single_parameter(param: str, param_range: list)-> (pd.DataFrame,pd.DataFrame):
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
        metric_df = pd.concat([metric_df, pd.DataFrame([new_row])], ignore_index = True)
    return(gfp_df, metric_df)


# solve qcrispri model while varying 2 parameters and calcuate dynamics parameters and compile gfp data into a df
def vary_two_parameter(param1: str, param1_range: list, param2: str, param2_range: list) -> (pd.DataFrame, pd.DataFrame):
    """solve qcrispri model and calculate switch density, foldrepression, transition time while varying two parameters at once

    Args:
        param1: name of first parameter to vary
        param1_range: list of values of first parameter to use 
        param2: name of second parameter to vary
        param2_range: list of values of second parameter to use
    """
    metric_df = pd.DataFrame(columns=[param1, param2, 'switching_density', 'fold_repression', 'transition_time', 'max_gfp'])
    gfp_df = pd.DataFrame(columns = [param1, param2, 'control', 'GFP', 'N' ])
    temp = pd.DataFrame()
    params = {'alphaT': alpha0,
              'activity_ratio': 0.1,
              'Kd1': kd1,
              'Kd2': kd2}
    
    for i in param1_range:
        params[param1] = i
        for j in param2_range:
            params[param2] = j
            qcrispri_df = solve_qcrispri(**params)
            #calculate dynamic parameters
            switch_density, fold_repression, transition_time = CalculateDynamicsMetrics(qcrispri_df.GFP, t, qcrispri_df.N).get_dynamic_metrics()

            max_gfp = max(qcrispri_df.GFP)
            new_row = {param1:i, param2:j, 'switching_density': switch_density, 'fold_repression': fold_repression, 'transition_time':transition_time, 'max_gfp':max_gfp}
            metric_df = pd.concat([metric_df, pd.DataFrame([new_row])], ignore_index = True)
            
            #create gfp dataframe
            temp[param1] = [i]*len(t)
            temp[param2] = [j]*len(t)
            temp['control'] = qcrispri_df.GFPc
            temp['GFP'] = qcrispri_df.GFP
            temp['N'] = qcrispri_df.N
            gfp_df = pd.concat([gfp_df, temp])
    return(gfp_df, metric_df)

#vary decoy sites and calculate parameters and create a df with gfp and related values
def vary_decoysite(param_name:str, param_range: list):
    metric_df = pd.DataFrame(columns=[param_name, 'act_density', 'FoldRepression', 'Transition_time', 'MaxGFP'])
    gfp_df = pd.DataFrame(columns = [param_name,'control', 'GFP', 'N' ])
    temp = pd.DataFrame()
    params = {'alphaT': alpha0,
              'activity_ratio': 0.1,
              'Kd1': kd1,
              'Kd2': kd2,
              'decoy_sites': 10}
    for i in param_range:
        params[param_name] = i
        qcrispri_df = solve_qcrispri_sponge(**params)
    
        temp[param_name] = [i]*len(t)
        temp['control'] = qcrispri_df.GFPc
        temp['GFP'] = qcrispri_df.GFP
        temp['N'] = qcrispri_df.N
        gfp_df = pd.concat([gfp_df, temp])

        switch_density, fold_repression, transition_time = CalculateDynamicsMetrics(qcrispri_df.GFP, t, qcrispri_df.N).get_dynamic_metrics()

        new_row = {param_name:i, 'switching_density': switch_density, 'fold_repression': fold_repression, 'transition_time':transition_time, 'max_gfp': max(qcrispri_df.GFP)}
        metric_df = pd.concat([metric_df, pd.DataFrame([new_row])], ignore_index = True)
    return(gfp_df, metric_df)

#returns the default parameter value of some parameters
def get_default_param_value(param_name:str)-> float:
    default_params = {'alphaT': alpha0,
                      'activity_ratio': 0.1,
                      'Kd1': kd1,
                      'Kd2': kd2,
                      'Stringency': 10}
    return(default_params[param_name])

#gets the index of a default parameter in a dataframe

def get_default_param_index(param_range, param_name):
    default_val = get_default_param_value(param_name)
    return(param_range.index(default_val))
    