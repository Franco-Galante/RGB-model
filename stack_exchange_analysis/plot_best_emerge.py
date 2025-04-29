import os
from src.common import compute_tte_distribution, add_ccdf_curve_to_ax
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(
    context="talk",
    style="ticks",
    palette="deep",
    font="sans-serif",
    color_codes=True,
    font_scale=1.2,
    rc={
        'figure.facecolor': 'white',
        'font.family': 'sans-serif',
        'xtick.major.pad': 1.5,
        'ytick.major.pad': 1.5,
        'legend.fontsize': 18,
        'lines.markersize': 4,
        'lines.linewidth': 3.0,
        'lines.linestyle': '-',
        'lines.marker': 'o'
    }
)

time_dict = {'m': 'Minutes', 'h': 'Hours', 'd': 'Days'} # global


def preprocess_tte(data_root_dir, field, views, delta_unit_selection, filter=True):

    tte, _, _ = compute_tte_distribution(data_root_dir, field, views, delta_unit_selection, filter=filter)

    # --------------- Time needed to Emerge for the Best Answer ----------------

    # distribution of the time needed for the best answer from its publication
    # (answer timestamp) to stably become the best answer (time to emerge)

    print('\n', '---- Plotting time to emerge CCDFs ----'.center(80))

    # time for the best to stably become the best, , the 'time_to_emerge' field 
    # has been fixed above we just have to remove eventual negative differences
    tte['delta_emergence'] = tte['time_to_emerge'] - tte['best_init_time']
    tte['delta_emergence_seconds'] = tte['delta_emergence'].dt.total_seconds()

    pre_len = len(tte)
    tte = tte[tte['delta_emergence_seconds'] >= 0]
    print(f'Dropped {pre_len - len(tte)}/{pre_len} rows with negative delta emergence'.center(80))

    # convert to all possible time units
    tte['delta_emergence_days'] = tte['delta_emergence'].dt.total_seconds() / 86400
    tte['delta_emergence_minutes'] = tte['delta_emergence'].dt.total_seconds() / 60
    tte['delta_emergence_hours'] = tte['delta_emergence'].dt.total_seconds() / 3600

    return tte


# plots the CCDF of the values passed in 'positive_time_delta' and saves the figure
def ccdf_and_save(fig_file, positive_time_delta, lab, color, xlog=True, ylog=False, time='m', ccdf=True, ax=None, figsize=None, linestyle='-'):

    plot = True
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize) # create a new figure
    else:
        plot = False # do not plot when the axis is passed as argument

    add_ccdf_curve_to_ax(positive_time_delta, lab, ccdf=ccdf, ax=ax, color=color, linestyle=linestyle)

    if xlog:
        ax.set_xscale('log')
    if ylog:
        ax.set_yscale('log')

    time_str = time_dict[time].lower()
    y_lab = 'CDF' if not ccdf else 'CCDF'
    ax.set_ylabel(y_lab)
    ax.set_xlabel(f'Time Elapsed ({time_str})')
    ax.set_ylim([1e-5, 1])

    if plot:
        fig.tight_layout()
        fig.savefig(os.path.join('figs', fig_file))



# ------------------------------------ MAIN ------------------------------------

if __name__ == '__main__':

    viewcount = 0
    delta_unit_selection = time_dict['m'].lower() # minutes
    EXTENSION = 'pdf' # 'png'
    
    if not os.path.isdir('figs'): # create the 'fig' folder if not present
        os.mkdir('figs')
    if not os.path.isdir(os.path.join('figs', 'paper')): # easier access to paper figs
        os.mkdir(os.path.join('figs', 'paper'))
    if not os.path.isdir(os.path.join('figs', 'others')):
        os.mkdir(os.path.join('figs', 'others'))

    # ----------------------------- Preprocess TTE -----------------------------

    print('\n', '* PLOTTING RESULTS TOGETHER *'.center(80, '*'), '\n', flush=True)

    tte_cs   = preprocess_tte('results', 'cs', viewcount, delta_unit_selection, filter=True)
    tte_math = preprocess_tte('results', 'math', viewcount, delta_unit_selection, filter=True)

    TIME_LAST = "d" # due to granularity needs to be days
    unit_selection = time_dict[TIME_LAST].lower()


    # -------------------------------- Plotting --------------------------------

    cs_color   = '#a406c4'
    math_color = '#0574c4'

    # ......................... y-logarithmic CCDF plot ........................
    #                 (plot with improved appearance for the paper)

    fig, ax = plt.subplots(figsize=(8.5, 4.5))

    ccdf_and_save('None', tte_math[f'delta_emergence_{unit_selection}'],
                   'MathStackExchange', math_color, xlog=False, ylog=True, time=TIME_LAST, ax=ax, linestyle='--')
    ccdf_and_save('None', tte_cs[f'delta_emergence_{unit_selection}'],
                   'StackOverflow', cs_color, xlog=False, ylog=True, time=TIME_LAST, ax=ax)

    ax.legend(loc="lower left") # add the legend
    ax.grid(True, which='major', linewidth=0.5, linestyle='--', color='darkgrey')
    ax.grid(False, which='minor')
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(1.5)

    # add custom ticks
    ax.tick_params(axis='both', which='both', direction='in', width=0.5, color='black', pad=3.0)
    
    fig.savefig(os.path.join('figs', 'paper', f'time_to_become_best_ccdf_ylog.{EXTENSION}'), pad_inches=0, bbox_inches='tight')
    fig.savefig(os.path.join('figs', 'paper', f'time_to_become_best_ccdf_ylog.png'), pad_inches=0, bbox_inches='tight')
    

    # ............................... loglog plot ..............................

    fig, ax = plt.subplots(figsize=(8.5, 4.5))

    ccdf_and_save('None', tte_math[f'delta_emergence_{unit_selection}'],
                        'MathStackExchange', math_color, xlog=True, ylog=True, time=TIME_LAST, ax=ax)
    ccdf_and_save('None', tte_cs[f'delta_emergence_{unit_selection}'],
                     'StackOverflow', cs_color, xlog=True, ylog=True, time=TIME_LAST, ax=ax)
    
    ax.legend() # add the legend
    ax.grid(True, which='major', linewidth=0.5, linestyle='--', color='darkgrey')
    # add custom ticks
    ax.tick_params(axis='both', which='both', direction='in', width=0.5, color='black', pad=3.0)

    fig.savefig(os.path.join('figs', 'others', f'time_to_become_best_ccdf_loglog.{EXTENSION}'), bbox_inches='tight')


    # ............................... linear plot ..............................

    fig, ax = plt.subplots(figsize=(8.5, 4.5))

    # linear (y-axes) CCDF plot
    ccdf_and_save('None', tte_cs[f'delta_emergence_{unit_selection}'],
                     'StackOverflow', cs_color, xlog=False, ylog=False, time=TIME_LAST, ax=ax)
    
    ccdf_and_save('None', tte_math[f'delta_emergence_{unit_selection}'],
                        'MathStackExchange', math_color, xlog=False, ylog=False, time=TIME_LAST, ax=ax)
    
    ax.legend() # add the legend
    ax.grid(True, which='major', linewidth=0.5, linestyle='--', color='darkgrey')
    # add custom ticks
    ax.tick_params(axis='both', which='both', direction='in', width=0.5, color='black', pad=3.0)

    fig.savefig(os.path.join('figs', 'others', f'time_to_become_best_ccdf.{EXTENSION}'), bbox_inches='tight')

    # NOTE: here I did not trim the distribution (raw distribution)
    print('\n', '* END *'.center(80, '*'), '\n')
