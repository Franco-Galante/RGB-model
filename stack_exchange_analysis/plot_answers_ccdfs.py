import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from src.common import *
from src.axes_utils import *
sns.reset_defaults()
plt.close()
sns.set_theme(
    style="ticks",
    palette="deep",
    font="sans-serif",
    color_codes=True,
    font_scale=1.1,
    rc={
        'figure.facecolor': 'white',
        'font.family': 'sans-serif',
        'xtick.major.pad': 1.5,
        'ytick.major.pad': 1.5,
        'lines.markersize': 4,
        'lines.linewidth': 2.0,
        'lines.linestyle': '-',
        'lines.marker': 'o'
    }
)

# --------------------------- Utilities ----------------------------------------

# reads from a Pandas dataframe column, a timestamp difference, returns a Series
def read_pd_delta_column(file_name, column_name):
    df = pd.read_csv(file_name)
    df[column_name] = pd.to_timedelta(df[column_name])
    return df[column_name]


def process_time_delta(time_delta, time='d', verbose=False):
    # check for negative values (errrors)
    positive_time_delta =  time_delta[time_delta >= pd.Timedelta(0)].tolist()
    negative_deltas = len(time_delta) - len(positive_time_delta)
    if negative_deltas > 0 and verbose:
        print(f'WARNING: {negative_deltas} negative time deltas found'.center(80), '\n', flush=True)

    # convert the delta time to minutes ('m'), hours ('h') or days ('d')
    if time == 'm':
        positive_time_delta = [v.total_seconds() / 60 for v in positive_time_delta]
    elif time == 'd':
        positive_time_delta = [v.total_seconds() / 86400 for v in positive_time_delta]
    elif time == 'h':
        positive_time_delta = [v.total_seconds() / 3600 for v in positive_time_delta]
    else:
        raise ValueError('Invalid time unit. Use "m", "h" or "d"'.center(80))
    return positive_time_delta


def trim_outliers(data, frac=0.01):
    sorted_data = np.sort(data)
    # trim the top-frac/2 and bottom-frac/2 of the data
    lower_bound = int(len(sorted_data) * frac / 2)
    upper_bound = int(len(sorted_data) * (1 - frac / 2))
    return sorted_data[lower_bound:upper_bound]


def set_axis_common_properties(ax, time_unit_label, y_lab='CCDF', xlab='Time'):
    """ Set common properties for axis,  to be reused for all figures in this script """

    ax.set_xlabel(f'{xlab} ({time_unit_label})')
    ax.set_ylabel(f'{y_lab}')

    # ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='major', axis='both', direction='in', width=0.5)
    ax.tick_params(which='minor', axis='both', direction='in', width=0.5, size=3)

    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(0.5)

    ax.grid(True, which="major", axis='both', linewidth=0.5, linestyle='--', color='darkgrey')


# ------------------------------------ MAIN ------------------------------------

if __name__ == '__main__':

    # ------------------------------- Load data --------------------------------
    res_root_dir = 'results'
    distrib_root_dir = 'distrib'

    argparser = argparse.ArgumentParser(description='Plotting CCDFs for time to first answer, best answer and emergence')
    argparser.add_argument('--dataset', type=str, default='cs', help='Dataset to use (cs or math)')
    argparser.add_argument('--viewcount', type=int, default=0, help='Filter by viewcount')
    args = argparser.parse_args()

    time_dict = {'m': 'Minutes', 'h': 'Hours', 'd': 'Days'}
    TIME = 'm'  
    unit_selection = time_dict[TIME].lower()
    REM_FRAC = 0.00

    COLOR_LIST = ['#f35b88', '#0094ff', '#00aea4']
    LINESTYLE_LIST = ['-', (0, (1,0.5)), (0, (3,1,1,1))]


    # 'Time to Emerge'
    tte, positive_time_diff_minutes, mask = compute_tte_distribution(res_root_dir, args.dataset, args.viewcount, unit_selection, filter=True)

    # 'Time to First Answer'
    file_root = 'time_to_first_answer' # also column name
    ttfa = read_pd_delta_column(
        os.path.join(distrib_root_dir, args.dataset, f'{file_root}.csv'), 
        file_root)
    if mask is not None:
        ttfa = ttfa[mask]
    positive_time_delta_ttfa = process_time_delta(ttfa, TIME)

    # Time to Best Answer
    file_root = 'time_to_best_answer' # also column name
    ttba = read_pd_delta_column(
        os.path.join(distrib_root_dir, args.dataset, f'{file_root}.csv'), 
        file_root)
    if mask is not None:
        ttba = ttba[mask]
    positive_time_delta_ttba = process_time_delta(ttba, TIME)


    # for each trace let us remove the top 0.005% and bottom 0.005% of the data (if REM_FRAC=0.01)
    trimmed_ttfa = trim_outliers(positive_time_delta_ttfa, frac=REM_FRAC)
    trimmed_ttba = trim_outliers(positive_time_delta_ttba, frac=REM_FRAC)
    trimmed_diff_minutes = trim_outliers(positive_time_diff_minutes, frac=REM_FRAC)


    # ---------------------------------- Plotting ----------------------------------

    # plot settings
    unit_selection = time_dict[TIME].lower()
    savefig=True
    figsize_one_column = (6, 3)

    fig, ax = plt.subplots(figsize=figsize_one_column)

    # plot CCDF curves 
    lw = plt.rcParams['lines.linewidth']
    lw = lw * .85
    add_ccdf_curve_to_ax(trimmed_ttfa, 'From question to first answer posted', ax=ax, linestyle=LINESTYLE_LIST[2], linewidth=lw, color=COLOR_LIST[0])
    add_ccdf_curve_to_ax(trimmed_ttba, 'From question to best answer posted', ax=ax, linestyle=LINESTYLE_LIST[1], linewidth=lw, color=COLOR_LIST[1])
    add_ccdf_curve_to_ax(trimmed_diff_minutes, 'From question to best answer established', ax=ax, linestyle=LINESTYLE_LIST[0], linewidth=lw, color=COLOR_LIST[2])

    # set a bunch of properties for the axis which can be reused in other figures with different scales on axis
    set_axis_common_properties(ax, 
                            time_unit_label=unit_selection,
                            xlab='Time Elapsed')

    # set properties specific to this figure
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(1e-4, 1.1)
    ax.set_xlim(1, 1e7)
    ax.legend(ncol=1, loc='lower left')
    ax2 = add_secondary_xaxis(ax, [60*24, 60*24*30, 60*24*30*12], ['1 day', '1 month', '1 year'], rotation=15, fontsize=10)
    ax2.tick_params(direction='in', width=0.5, size=10)
    no_minor_ticks(ax2, y=False)
    for spine in ax2.spines.values():
        spine.set_color('black')
        spine.set_linewidth(0.5)


    # plot legend in file figs/paper/legend.pdf
    #   - frameon = False will be passed to ax.legend() and avoids the frame in the legend box
    #   - remove_legend_from_original also removes the legend from the original plot
    if not os.path.exists(os.path.join('figs', 'paper')):
        os.makedirs(os.path.join('figs', 'paper'))
    legend_alone(ax, 
                legend_figsize=(figsize_one_column[0], 1), 
                saveto=os.path.join('figs', 'paper','legend.pdf'),
                remove_legend_from_original=True,
                ncol=3,
                loc='center',
                frameon=False)

    myplot(fig, os.path.join('figs', 'paper', f'summary_ccdf_loglog_{args.dataset}'), savefig=savefig)
