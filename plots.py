import matplotlib.pyplot as plt
from Calc_parameters import *
from solve_ode import *

#plots vertical stack of 3 plots. For sensitivity analysis
def vStack_3plots(df, param, xlabel):
    fig, (ax1, ax2, ax3) = plt.subplots(3,1, figsize = (4, 8))
    ax1.plot(df[param], df.act_density, "C1o-", mec = '1.0', color = 'crimson', linewidth = 2)
    ax1.plot(get_default_params(param), df.act_density.iloc[get_paramIndex(df, param)], marker = "o", color = 'crimson')
    ax1.set_ylabel('Switching Density (CFU/mL)')

    ax2.plot(df[param], df.FoldRepression, "C1o-", mec = '1.0', color = 'C2', linewidth = 2)
    ax2.plot(get_default_params(param), df.FoldRepression.iloc[get_paramIndex(df, param)], marker = "o", color = 'C2')
    ax2.set_ylabel('Fold Repression')

    ax3. plot(df[param], df.Transition_time, "C1o-", mec = '1.0', color = 'dodgerblue', linewidth = 2)
    ax3.plot(get_default_params(param), df.Transition_time.iloc[get_paramIndex(df, param)], marker = "o", color = 'dodgerblue')
    ax3.set_xlabel(xlabel)
    ax3.set_ylabel('Transition time (hr)')
    fig.tight_layout()
    return(fig)


#2 parameter varied sensitivity analysis data and plots y vs x for different leakiness. 
def plot2param_vary(params_df, Xparam, Yparam, xlabel = 'xlabel', ylabel = 'ylabel' ):
    fig, ax = plt.subplots()
    for i in [0, round(alpha0*0.5, 1), alpha0, round(alpha0*1.5, 1), round(alpha0*2, 1)  ]:
        subdf = params_df[params_df['alphaT'] == i]
        ax.plot(subdf[Xparam], subdf[Yparam], linewidth = 2, label = '%s' %i);
    for i in ['right', 'left', 'top', 'bottom']:
        ax.spines[i].set(linewidth = 1.25)
    ax.legend(prop={"size":12})
    plt.xlabel(xlabel, fontsize = 14)
    plt.ylabel(ylabel,  fontsize = 14)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    plt.legend(title = 'Leakiness', fancybox = True )
    plt.rc('legend', fontsize=8)
    return(plt)