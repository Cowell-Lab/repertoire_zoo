"""
Module Name: giraffe

**G**raphical **I**nteractive **R**epertoire **A**nalysis with **F**ast **F**eature **E**xploration.
This module provides plotting utilities for analyzing repertoire data.

Functions:
- plot_duplicate_frequency_average(df_combined, colors, figsize, title) : Plot
  grouped bar chart with error bars from combined dataframe.
  Returns a tuple containing the matplotlib figure and axes objects.
- plot_duplicate_frequency(df_combined, colors, figsize, title) : Plot grouped
  bar chart with error bars from combined dataframe. Returns a tuple containing
  the matplotlib figure and axes objects.
- plot_junction_aa_length : Plot grouped bar chart with error bars from combined
  dataframe. Returns a tuple containing the matplotlib figure and axes objects.

Usage:
For plotting repertoire and repertoire group data from a Pandas DataFrame. This
module is to be used in conjunction with data loaded from the `bird` module.
Import this module and call the desired functions with appropriate arguments.
"""

from matplotlib.axes import Axes
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pycirclize import Circos
import seaborn as sns

def plot_gene_usage_groups(
        df_combined: pd.DataFrame,
        relative:bool=True,
        abundance:bool=True,
        colors:list=None,
        figsize:tuple=(16,8),
        title:str=None
    ) -> tuple[Figure, Axes]:
    """
    Plot grouped bar chart with error bars from combined dataframe.

    Parameters
    ----------
    df_combined : pd.DataFrame
        - A Pandas DataFrame with columns `gene`, `duplicate_frequency_avg`,
        `duplicate_frequency_std`, `condition`.
    relative : bool
        - Flag to plot either relative frequencies (True) or absolute counts (False).
    abundance : bool
        - Flag to plot either abundance (True) or sequence (False) values.
    colors : list
        - A list of colors for each condition. (optional, Default: `None` uses a
        predefined color palette. This currently supports up to 5 groups).
    figsize : tuple
        - The figure size. (Default:`(16,8)`)
    title : str
        - The plot's title (optional, Default: `None`)

    Returns
    -------
    fig, ax : tuple
        - Matplotlib figure and axes objects
    """
    genes = df_combined['gene'].unique()
    conditions = df_combined['condition'].unique()

    if colors is None:
        ## Add more colors if needed
        default_palette = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854']
        colors = default_palette[:len(conditions)]

    n_groups = len(conditions)
    bar_width = 0.8 / n_groups
    x = np.arange(len(genes))

    fig, ax = plt.subplots(figsize=figsize)

    for i, cond in enumerate(conditions):
        cond_data = df_combined[df_combined['condition'] == cond]
        cond_data = cond_data.set_index('gene').reindex(genes, fill_value=0).reset_index()

        positions = x + i * bar_width

        if relative:
            if abundance:
                ylabel = 'Average Duplicate Frequency'
                mean = cond_data['duplicate_frequency_avg']
                std = cond_data['duplicate_frequency_std']
            else:
                ylabel = 'Average Sequence Frequency'
                mean = cond_data['sequence_frequency_avg']
                std = cond_data['sequence_frequency_std']
        else:
            if abundance:
                ylabel = 'Average Duplicate Count'
                mean = cond_data['duplicate_count_avg']
                std = cond_data['duplicate_count_std']
            else:
                ylabel = 'Average Sequence Count'
                mean = cond_data['sequence_count_avg']
                std = cond_data['sequence_count_std']
        lower_error = np.minimum(std, mean)
        upper_error = std
        asymmetric_error = [lower_error, upper_error]

        ax.bar(positions,
               mean,
               bar_width,
               label=cond,
               yerr=asymmetric_error,
               capsize=4,
               color=colors[i],
               error_kw={'elinewidth':1, 'ecolor':'black'})

    ax.set_xticks(x + bar_width * (n_groups - 1) / 2)
    ax.set_xticklabels(genes, rotation=90)
    ax.set_ylabel(ylabel)
    ax.set_xlabel('Gene Family')
    ax.grid(True)
    ax.legend(title='Condition', loc='upper left', bbox_to_anchor=(1, 1))

    if title is None:
        title = 'TCR Family Usage: ' + ' vs '.join(conditions)
    ax.set_title(title)

    plt.tight_layout()
    return fig, ax

def plot_duplicate_frequency(
        df:pd.DataFrame,
        average:bool=False,
        datapoints:bool=True,
        figsize:tuple=(16,8),
        title:str=None,
        hue=None,
        hue_order=None,
        palette:dict=None,
        threshold:float=None,
        ylim:tuple[float,float]=(None,None),
        ax:Axes=None,
        split:bool=False
    ) -> tuple[Figure, Axes]:
    """
    Plot grouped bar chart with error bars from combined dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        - A Pandas DataFrame with columns `gene`, `duplicate_frequency_avg`,
        `duplicate_frequency_std`, `condition`
    average : bool
        - Displays the average and the error bars. (Default: `True`.)
    datapoints : bool
        - Displays the data points on the bars. (Default: `True`.)
    palette : list
        - A list of colors for each condition. (optional, Default: `None`.)
        uses a predefined color palette. This currently supports up to 5 groups)
    figsize : tuple
        - The figure size. (Default:`(16,8)`)
    title : str
        - The plot's title (optional, Default: `None`)
    threshold : float
        - Define a duplicate frequency threshold to apear. Values are inclusive. (Default: `None`.)
    ylim : tuple[float, float]
        - Define the y-axim limit. If a tuple is provided, then read as [min, max]. 
        It is also possible to define the tuple using `None` to denote automatic limits. 
        (Default: `(None, None)`.)
    ax : Axes
        - A `matplotlib.axes.Axes` object to plot with. (Deafult: `None`.)

    Returns
    -------
    fig, ax : tuple
        - Matplotlib figure and axes objects
    """
    external_ax = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()
        external_ax = True


    # Determine hue order
    if hue_order is None and isinstance(palette, dict):
        hue_order = list(palette.keys())
    if hue is not None and palette is None:
        palette = 'Set2'
    if hue is None:
        hue='gene'

    # reduce df by threshold
    if threshold is not None and isinstance(threshold, float):
        df = df[df['duplicate_frequency']>=threshold]

    

    if not split:
        sns.barplot(
            data=df,
            x='gene',
            y='duplicate_frequency',
            hue=hue,
            legend=False,
            palette=palette,
            hue_order=hue_order,
            errorbar=None,
            ax=ax
        ) 
    else:
        sns.barplot(
            data=df,
            x='gene',
            y='duplicate_frequency',
            hue='condition',
            errorbar=None,
            ax=ax
        )
        ax.legend(title='Category', loc='upper left', fontsize='small')

    # sns.pointplot(
    #     data=df,
    #     linestyle='none',
    #     x='gene',
    #     y='duplicate_frequency',
    #     color='black',
    #     zorder=1
    # )

    ax.tick_params(axis='x', rotation=90)
    # ax.grid(True)
    # ax.legend(title='Tissue', loc='upper left', bbox_to_anchor=(1, 1))

    if title:
        ax.set_title(title, size=20)
    
    if ylim[0] is not None and isinstance(ylim[0], float):
        plt.ylim(bottom=ylim[0])
    if ylim[1] is not None and isinstance(ylim[1], float):
        plt.ylim(top=ylim[1])

    if external_ax:
        plt.tight_layout()

    return fig, ax

def plot_junction_aa_length(
        df: pd.DataFrame,
        x_col:str='CDR3_LENGTH',
        y_col:str='CDR3_COUNT',
        hue_col:str='condition',
        palette:str='Set1',
        figsize:tuple[int,int]=(18, 8),
        log_scale:bool=False,
        aa_range:tuple[int,int]=(6, 29),
        ax:Axes=None
    ) -> tuple[Figure, Axes]:
    """
    Plot grouped bar chart with error bars from combined dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        - A Pandas DataFrame with columns `gene`, `duplicate_frequency_avg`,
        `duplicate_frequency_std`, `condition`
    x_col : str
        - The x-variable column name. (Default: `CDR3_LENGTH`)
    y_col : str
        - The y-variable column name. (Default: `CDR3_COUNT`)
    hue_col : str
        - The coloumn name to determine the hue color. (Default: `condition`)
    palette : str
        - The seaborn color palette used. (Default: `Set1`)
    figsize : tuple
        - The figure size. (Default:`(18,8)`)
    log_scale : bool
        - Place the y axis into a log scale. (Default: `False`)
    aa_range : tuple[int,int]
        - The range of amino acid positions to display. (Default: `(6,29)`.)
    Returns
    -------
    fig, ax : tuple
        - Matplotlib figure and axes objects
    """
    external_ax = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()
        external_ax = True
    
    sns.barplot(data=df, x=x_col, y=y_col, hue=hue_col, palette=palette, ax=ax)
    if log_scale:
        ax.set_yscale('log')
    plt.xticks(rotation=90)
    plt.tight_layout()
    ax.set_xlim(aa_range[0], aa_range[1])

    if external_ax:
        plt.tight_layout()

    return fig, ax

def plot_vj_chord_diagram(
        df:pd.DataFrame,
        freq_col:str="duplicate_frequency_avg",
        ticks_interval:float=0.03,
        space:int=2,
        cmap:str='Set1',
        figsize:tuple=(12, 12),
        max_links:int=None,
        title:str=None,
        ax:Axes=None
    ) -> tuple[Figure, Axes]:
    """
    Plot the VJ Chord Diagram.

    Parameters
    ----------
    df : pd.DataFrame
        - The dataframe from `bird.load_and_prepare_data_vj_combo`.
    freq_col : str
        - The name of the frequency column. (Default: `duplicate_frequency_avg`.)
    ticks_interval : float
        - The interval between each tick on the plot. (Default: `0.03`.)
    space : int
        - The spacing on the diagram. (Default: `2`.)
    cmap : str
        - The colormap used. (Default: `Set1`)
    figsize : tuple[int, int]
        - The size of the figure. (Default: `(12, 12)`.)
    max_links : int
        - The maximum nunber of links. (Default: `None`.)
    title : str
        - The title of the plot. (Default: `None`.)
    
    Returns
    -------
    fig, ax : tuple[Figure, Axes]
        - Matplotlib Figure and Axes objects.
    """
    external_ax = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, subplot_kw={"projection": "polar"})
    else:
        fig = ax.get_figure()
        external_ax = True

    # Optionally reduce overcrowding
    if max_links and len(df) > max_links: ## Limit to top N V–J links
        df = df.nlargest(max_links, freq_col)

    # Build matrix
    nodes = sorted(set(df["v_level"]).union(set(df["j_level"])))
    matrix_df = pd.DataFrame(0, index=nodes, columns=nodes, dtype=float)
    
    for _, row in df.iterrows():
        v, j = row["v_level"], row["j_level"]
        matrix_df.loc[v, j] = row[freq_col]
    
    circos = Circos.chord_diagram(
        matrix_df,
        space=space,  # more spacing
        ticks_interval=ticks_interval,
        cmap=cmap,
        label_kws=dict(
            size=10,
            r=115,
            orientation='vertical'),
        link_kws=dict(
            direction=0,
            ec="white",
            lw=0.4),
        ticks_kws=dict(
            label_size=8,
            label_orientation='vertical')
    )

    if title:
        if external_ax:
            text_size = 15
        else:
            text_size = 20
        circos.text(
            title,
            size=text_size,
            deg=0,
            r=140
        )

    if not external_ax:
        plt.tight_layout()

    circos.plotfig(figsize=figsize, ax=ax)

    return fig, ax

def plot_diversity_curve(df, figsize=(16,8)):
    fig, ax = plt.subplots(figsize=(16,8))
    ax.plot(df['q'], df['d'])
    ax.set_xscale('log')
    ax.grid(axis='x')
    ax.fill_between(df['q'], df['d_lower'], df['d_upper'])
    return fig, ax

def plot_mutational_hedgehog(
        df:pd.DataFrame,
        pad:int=5,
        num_intervals:int=4,
        start:float=0.05,
        pallete:str='tab20',
        figsize:tuple[int,int]=(10, 10),
        title:str=None,
        y_max:int=None
    ) -> tuple[Figure, Axes]:
    """
    Plot grouped bar chart with error bars from mutational Data.

    Parameters
    ----------
    df : pd.DataFrame
        - DataFrame with columns: `value`, `group`, `place`, `se`, `se_start`, `se_end`.
    pad : int
        - Spacing between each group. (Default: `5`.)
    num_intervals : int
        - An integer number for showing the intervals for y limit (Default: `4`.)
    start : float
        - A float for starting number for the interval (Default: `0.05`.)
    pallete : str
        - Seaborn color pallete (Default: `tab20`.)
    figsize : tuple[int, int]
        - The figure size. (Default: `(10, 10)`.)
    title : str
        - The title of the plot. (Default: `None`.)
    y_max : int
        - The maximum y value displayed. Can be used to set plots to the same scale. (Default: `None`.)

    Returns
    -------
    fig, ax : tuple[Figure, Axes]
        - Matplotlib Figure and Axes objects.
    """
    values = df["value"].values
    groups = df["group"].values

    #Keep the group orientation as is.
    unique_groups = pd.unique(groups)
    groups_sizes = [len(group_df) for _, group_df in df.groupby("group", sort=False)]

    # Calculate max y value
    if y_max is None:
        y_max = np.max(values)
    y_min = -0.1

    # Calculate angles and width
    num_angles = len(values) + pad * len(unique_groups)
    angles = np.linspace(0, 2 * np.pi, num=num_angles, endpoint=False)
    width = (2 * np.pi) / len(angles)

    # Choose a seaborn palette
    my_palette = sns.color_palette(pallete, len(unique_groups))
    # Assign colors to groups based on the palette
    group_to_color = dict(zip(unique_groups, my_palette))
    colors = [group_to_color[group] for group in groups]

    # Generate points from start to y_max, with num_intervals steps
    ref_pts = np.linspace(start, y_max, num_intervals)
    ref_pts_rounded = np.round(ref_pts, 2)
    ref_pt = [[float(x)] for x in ref_pts_rounded]

    #start plotting
    fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(polar=True))
    offset = 0
    idxs = []
    for size in groups_sizes:
        idxs += list(range(offset + pad, offset + size + pad))
        offset += size + pad

    # Determines where to place the first bar.
    # By default, matplotlib starts at 0 (the first bar is horizontal)
    # but here we say we want to start at pi/2 (90 deg)
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1) #Change the plot direction from anticlockwise to clockwise
    ax.set_ylim([y_min, y_max])
    ax.bar(
        angles[idxs],
        values,
        width=width,
        yerr=df.se,
        color=colors,
        edgecolor="white",
        linewidth=1,
        error_kw={'elinewidth': 2}
    )

    # This iterates over the sizes of the groups adding reference
    offset = 0
    delta = 0.01
    for group, size in zip(unique_groups, groups_sizes):
        #for each group minimum and maximum id
        min_idx, max_idx = df[df.group == group].place.min(), df[df.group == group].place.max()

        # Add line below bars
        x1 = np.linspace(angles[offset + pad], angles[offset + size + pad - 1], num=50)
        ax.plot(x1, [-delta] * 50, color="black", linewidth = 2)   

        # Add text to indicate group
        ax.text(
            np.mean(x1),
            -3*delta,
            group,
            color="black",
            fontsize=10,
            fontweight="bold",
            ha="center",
            va="center"
        )

        # Min index label slightly left and below first bar
        ax.text(
            x1[0],
            -2*delta,
            str(min_idx),
            size=9,
            alpha=1,
            ha='center',
            va='top'
        )
        # Max index label slightly right and below last bar
        ax.text(
            x1[-1],
            -2*delta,
            str(max_idx),
            size=9,
            alpha=1,
            ha='center',
            va='top'
        )

        # Add reference lines at
        x2 = np.linspace(angles[offset], angles[offset + pad - 1], num=50)
        for rp in ref_pt:
            ax.plot(
                x2,
                rp*50,
                color="gray",
                lw=1
            )
            ax.text(
                0.09,
                rp[0]-delta,
                f'{rp[0]}',
                color="black",
                ha='center',
                va='top',
                fontsize=9
            )
        offset += size + pad

    ax.set_frame_on(False)
    ax.xaxis.grid(False)
    ax.yaxis.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    if title:
        ax.set_title(title, fontdict={'fontsize':20})
    plt.tight_layout()
    return fig, ax
