import pandas as pd
from np2df import *
from solve_ode import *

# Cell Density at the end of growth phase. 

def switch_density(data,t):
    max_gfp = max(data.GFP)
    max_gfp_index = data[data['GFP'] == max_gfp ].index[0] 
    max_90_id= min(range(max_gfp_index,len(data.GFP)), key=lambda i: abs(data.GFP[i]- (max_gfp*0.9)))
    act_density = data.N[max_90_id]
    return(act_density)

def switch_density_generic(gfp, N):
    max_gfp = max(gfp)
    max_gfp_index = list(gfp).index(max_gfp)
    max_90_id= min(range(max_gfp_index,len(gfp)), key=lambda i: abs(gfp[i]- (max_gfp*0.9)))
    act_density = N[max_90_id]
    return(act_density)



# To calculate fold repression Max GFP/ min GFP

def foldRepression(data):
    max_gfp = max(data.GFP)
    max_gfp_index = data[data['GFP'] == max_gfp].index[0]
    min_gfp = min(data.GFP[max_gfp_index:])
    fold  = max_gfp/min_gfp
    return(fold)

def foldRepression_generic(gfp):
    max_gfp = max(gfp)
    max_gfp_index = list(gfp).index(max_gfp)
    min_gfp = min(gfp[max_gfp_index:])
    fold  = max_gfp/min_gfp
    return(fold)



#Calculate transition time from max_cell density to min_cell density

def transition_time(data,t):
    max_gfp = max(data.GFP)
    max_gfp_index = data[data['GFP'] == max_gfp ].index[0]
    min_gfp = min(data.GFP[max_gfp_index:])
    min_10_id= min(range(max_gfp_index,len(data.GFP)), key=lambda i: abs(data.GFP[i]- (min_gfp*1.02)))
    trans_time = t[min_10_id] - t[max_gfp_index]
    return(trans_time)

# In experimental the slow drop in much longer and does not come to a steady state as smoothly as simulated. So the threshold for minimum is a little higher. 
def transition_time_generic(gfp,t):
    max_gfp = max(gfp)
    max_gfp_index = list(gfp).index(max_gfp)
    min_gfp = min(gfp[max_gfp_index:])  
    min_10_id= min(range(max_gfp_index,len(gfp)), key=lambda i: abs(gfp[i]- (min_gfp*1.15)))
    trans_time = t[min_10_id] - t[max_gfp_index]
    return(trans_time)



#solve qcrispri model and calculate switch density, foldrepression, transition time
# for a range of values of : alphaT, activity_ratio, kDNA_bind, kcomplex

def single_param_vary(param, param_range):
    param_df = pd.DataFrame(columns=[param, 'act_density', 'FoldRepression', 'Transition_time'])
    gfp_df = pd.DataFrame(columns = [param,'control', 'GFP', 'N' ])
    temp = pd.DataFrame()
    params = {'alphaT': alpha0,
              'activity_ratio': 0.1,
              'Kd1': kd1,
              'Kd2': kd2}
    qcrispridf = convert2df()
    for i in param_range:
        params[param] = i
        qcrispridf.qcrispri(solve_qCRISPRi(**params))
        qcrispri_df = qcrispridf.returndf()
        
        temp[param] = [i]*len(t)
        temp['control'] = qcrispri_df.GFPc
        temp['GFP'] = qcrispri_df.GFP
        temp['N'] = qcrispri_df.N
        gfp_df = pd.concat([gfp_df, temp])
        
        s_density = switch_density(qcrispri_df, t)
        fold_rep = foldRepression(qcrispri_df)
        t_time = transition_time(qcrispri_df, t)
        new_row = {param:i, 'act_density': s_density, 'FoldRepression': fold_rep, 'Transition_time':t_time}
        param_df = param_df.append(new_row, ignore_index = True)
    return(gfp_df, param_df)


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
        qcrispridf.qcrispri_sponge(solve_qCRISPRi_sponge(**params))
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

def get_default_params(param):
    default_params = {'alphaT': alpha0,
                      'activity_ratio': 0.1,
                      'Kd1': kd1,
                      'Kd2': kd2,
                      'Stringency': 10}
    return(default_params[param])

#gets the index of a default parameter in a dataframe

def get_paramIndex(df, param):
    key = get_default_params(param)
    return(df[df[param] == key].index[0])
    