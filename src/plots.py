import matplotlib.pyplot as plt
import pandas as pd
from .sensitivity_analysis import get_default_param_value, get_default_param_index
from .solve_ode import alpha0


def vertical_stack_plot_3(param_name: str, param_range: pd.Series, switching_density: pd.Series, fold_repression: pd.Series, transition_time: pd.Series, xlabel: str):
    """plots vertical stack of 3 plots. For sensitivity analysis

    Args:
        param_name: name of the parameter varied
        param_range: range of parameter values used for analysis
        switching_density: list of switching density values for different values of parameter
        fold_repression: list of fold repression values for different values of parameter
        transition_time: list of transition time values for different values of parameter
        xlabel: x label for the plot
    """
    fig, (ax1, ax2, ax3) = plt.subplots(3,1, figsize = (4, 8))
    
    ax1.plot(param_range, switching_density, "C1o-", mec = '1.0', color = 'crimson', linewidth = 2)
    ax1.plot(get_default_param_value(param_name), switching_density.iloc[get_default_param_index(param_range, param_name)], marker = "o", color = 'crimson')
    ax1.set_ylabel('Switching Density (CFU/mL)')

    ax2.plot(param_range, fold_repression, "C1o-", mec = '1.0', color = 'C2', linewidth = 2)
    ax2.plot(get_default_param_value(param_name), fold_repression.iloc[get_default_param_index(param_range, param_name)], marker = "o", color = 'C2')
    ax2.set_ylabel('Fold Repression')

    ax3. plot(param_range, transition_time, "C1o-", mec = '1.0', color = 'dodgerblue', linewidth = 2)
    ax3.plot(get_default_param_value(param_name), transition_time.iloc[get_default_param_index(param_range, param_name)], marker = "o", color = 'dodgerblue')
    ax3.set_xlabel(xlabel)
    ax3.set_ylabel('Transition time (hr)')
    fig.tight_layout()
    return(fig)


#2 parameter varied sensitivity analysis data and plots y vs x for different leakiness. 
def plot_2param_vary(metrics_df, Xparam, Yparam, xlabel = 'xlabel', ylabel = 'ylabel' ):
    fig, ax = plt.subplots()
    for i in [0, round(alpha0*0.5, 1), alpha0, round(alpha0*1.5, 1), round(alpha0*2, 1)  ]:
        subdf = metrics_df[metrics_df['alphaT'] == i]
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