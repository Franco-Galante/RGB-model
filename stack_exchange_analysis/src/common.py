import matplotlib.pyplot as plt
import pandas as pd
import os


def compute_tte_distribution(data_root_dir, field, views, delta_unit_selection, filter=True, verbose=False):

    tte = pd.read_csv(os.path.join(data_root_dir, field, f'time_to_emerge_{views}.csv'))
    init_tte = len(tte)

    # ensure correct timestamp format
    tte['time_to_emerge'] = pd.to_datetime(tte['time_to_emerge'])
    tte['question_time']  = pd.to_datetime(tte['question_time'])
    tte['best_init_time'] = pd.to_datetime(tte['best_init_time'])

    # I need to take care of the None values (the best has never been overcome),
    # thus the best has always been the best. And the negative timestamp deltas
    # (up to negative 1 day), due to the difference in timestamp granularity
    tte['time_to_emerge'] = tte['time_to_emerge'].fillna(tte['best_init_time'])

    # add a filter that allows to exclude the records whose "time to emerge" has
    # been in the past 6 moths (they may not yet be the 'best' answer)
    mask = None
    if filter:
        # timestamp of the most recent question in the dataset
        most_recent_question = tte['question_time'].max()
        print(f'Most recent question timestamp: {most_recent_question}'.center(80))
        # timestamp of the most recent question in the dataset - 6 months
        x_months_ago = most_recent_question - pd.DateOffset(months=6)
        # filter the records that are older than 6 months
        mask = tte['time_to_emerge'] < x_months_ago
        tte = tte[mask]


    # compute timestamp differences and adjust those that are negative up to 1 day
    tte['delta'] = tte['time_to_emerge'] - tte['question_time']
    
    df_TOT_neg_delta = tte[tte['delta'] < pd.Timedelta(0)]
    TOT_neg_delta = len(tte)
    diff = len(tte) - len(df_TOT_neg_delta)
    if verbose:
        print(f'WARNING: there are {diff}/{len(tte)} TOTAL negative deltas'.center(80))

    tte['delta_seconds'] = tte['delta'].dt.total_seconds() # convert to seconds
    new_diff = len(tte) - sum([1 for v in tte['delta_seconds'] if v < 0]) # needs to be the same
    if verbose:
        print(f'WARNING: there are {new_diff}/{len(tte)} negative deltas (same)'.center(80))
        print('--> This are due to a granularity mismatch in the timestamps and can be solved!'.center(80), '\n')
    
    
    # keep positive values as they are, discard timestamps smaller than -1 day
    # modify those that are between -1 day and 0 days (granularity mismatch)
    minus_1day_seconds = -1 * 86400

    pre_tte = len(tte)
    tte = tte[(tte['delta_seconds'] > minus_1day_seconds)] # discard row with more than 1 day delta
    smaller1d = pre_tte - len(tte)
    if verbose:
        print(f'Dropped {smaller1d}/{pre_tte} rows with delta smaller than -1 day'.center(80))

    # for the remaining negatives, that by construction are greater than -1 day
    # I substitute the best init time (so we have the time between the question 
    # and the best answer) and we artificially get a better granularity
    mask_negatives = (tte['delta_seconds'] < 0)
    tte.loc[mask_negatives, 'time_to_emerge'] = tte.loc[mask_negatives, 'best_init_time']

    # we need to recompute the deltas with the new timestamps
    tte['delta'] = tte['time_to_emerge'] - tte['question_time']
    tte['delta_seconds'] = tte['delta'].dt.total_seconds()      # convert to seconds

    tte['delta_days'] = tte['delta'].dt.total_seconds() / 86400 # convert to days
    tte['delta_minutes'] = tte['delta'].dt.total_seconds() / 60 # convert to minutes
    tte['delta_hours'] = tte['delta'].dt.total_seconds() / 3600 # convert to hours

    # how many deltas have been fixed and how many are still there
    df_still_negatives = tte[tte['delta_seconds'] < 0]
    tot_remaining_neg = len(df_still_negatives)
    over = TOT_neg_delta - smaller1d
    if verbose:
        print(f'WARNING: {tot_remaining_neg}/{over} negative deltas still unsolved'.center(80), '\n')
    
    df_still_negatives[['question_time', 'best_init_time', 'first_answer_time']].head(5)

    assert delta_unit_selection == 'minutes', "Only minutes are supported"
    positive_time_diff_minutes = [v for v in tte[f'delta_{delta_unit_selection}'] if v >= 0]

    print(f'# records plotted: {len(positive_time_diff_minutes)}/{init_tte} (all records)'.center(80))
    
    return tte, positive_time_diff_minutes, mask


def myplot(fig, fig_file, savefig=False, **kwargs):
    """ 
    Expects global variable `savefig` to decide whether 
    to save the figure or show it, if saving to file
    it saves both a .png and a .pdf version and for the
    pdf version crops margins, sets dpi 300 for quality
    """

    fig.tight_layout()
    if savefig:
        print(f'Saving to {fig_file}'.center(80))
        fig.savefig(fig_file+'.png', pad_inches=0, bbox_inches='tight', dpi=300, **kwargs)
        fig.savefig(fig_file+'.pdf', pad_inches=0, bbox_inches='tight', dpi=300, **kwargs)
    else:
        print(f'Plotting {fig_file}'.center(80))
        plt.show()
    plt.close()


def add_ccdf_curve_to_ax(positive_time_delta, lab, ccdf=True, ax=None, **kwargs):

    # compute the CCDF (handling duplicates occurences)
    data_series = pd.Series(positive_time_delta)
    data_counts = data_series.value_counts().sort_index() # return values and counts, sort by value
    cumulative_counts = data_counts.cumsum()
    out = cumulative_counts / cumulative_counts.iloc[-1]  # renormalize by tot occurences

    # data may be zero, so add the first point
    vals = data_counts.index.tolist()
    if vals[0] == 0:
        start = vals[0] - vals[1] # computed on the basis of the values
        vals = [start] + vals     # zero occurences on the data

        out = [0] + out.tolist() # add the first point
            
    elif vals[0] < 0:
        raise ValueError('Negative values in the data'.center(80))
    
    if ccdf:
        out = [1-o for o in out] # ccdf
    
    ax.plot(vals, out, marker='', label=lab, **kwargs)
