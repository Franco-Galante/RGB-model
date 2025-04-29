import os
from matplotlib import pyplot as plt


def no_minor_ticks(ax, x=True, y=True):
    """ force remove of minor ticks """
    if x:
        ax.xaxis.set_minor_formatter(plt.NullFormatter())
        ax.xaxis.set_minor_locator(plt.NullLocator())
    if y:
        ax.yaxis.set_minor_formatter(plt.NullFormatter())
        ax.yaxis.set_minor_locator(plt.NullLocator())


def add_secondary_xaxis(ax, xticks, xlabels, rotation=0, fontsize=12):
    """ creates second x-axis on top """
    ax2 = ax.twiny()
    ax2.set_xscale(ax.get_xscale())
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xlabels, rotation=rotation, fontsize=fontsize)
    return ax2


def legend_alone(ax, legend_figsize, saveto='legend.pdf', remove_legend_from_original = False, **kwargs):
    """ 
        reads legend from current axis and saves it to a new figure 
    """
    
    # then create a new image with the legend only
    fig_leg = plt.figure(figsize=legend_figsize, dpi=300)
    ax_leg = fig_leg.add_subplot(111)
    # add the legend from the previous axes
    ax_leg.legend(*ax.get_legend_handles_labels(), **kwargs)
    # hide the axes frame and the x/y labels
    ax_leg.axis('off')
    print(f"Saving legend to file {saveto}".center(80))
    fig_leg.savefig(saveto, bbox_inches='tight', pad_inches=0, transparent=True)
    fig_leg.savefig(os.path.join(os.path.dirname(saveto), 'legend.png'), bbox_inches='tight', pad_inches=0)
    if remove_legend_from_original:
        print("Removing legend from original axes".center(80))
        ax.get_legend().remove()
    plt.close()
