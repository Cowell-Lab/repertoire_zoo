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

def plot_duplicate_frequency_average(
        df_combined: pd.DataFrame,
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

        mean = cond_data['duplicate_frequency_avg']
        std = cond_data['duplicate_frequency_std']
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
    ax.set_ylabel('Average Duplicate Frequency')
    ax.set_xlabel('Gene Family')
    ax.grid(True)
    ax.legend(title='Condition', loc='upper left', bbox_to_anchor=(1, 1))

    if title is None:
        title = 'TCR Family Usage: ' + ' vs '.join(conditions)
    ax.set_title(title)

    plt.tight_layout()
    return fig, ax

def plot_duplicate_frequency(
        df_combined: pd.DataFrame,
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
        `duplicate_frequency_std`, `condition`
    colors : list
        - A list of colors for each condition. (optional, Default: `None`
        uses a predefined color palette. This currently supports up to 5 groups)
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
    conditions = df_combined['tissue'].unique()

    if colors is None:
        ## Add more colors if needed
        default_palette = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854']
        colors = default_palette[:len(conditions)]

    n_groups = len(conditions)
    bar_width = 0.8 / n_groups
    x = np.arange(len(genes))

    fig, ax = plt.subplots(figsize=figsize)

    for i, cond in enumerate(conditions):
        cond_data = df_combined[df_combined['tissue'] == cond]
        cond_data = cond_data.set_index('gene').reindex(genes, fill_value=0).reset_index()
        positions = x + i * bar_width
        freq_val = cond_data['duplicate_frequency']
        ax.bar(positions,
               freq_val,
               bar_width,
               label=cond,
               yerr=None,
               capsize=4,
               color=colors[i],
               error_kw={'elinewidth':1, 'ecolor':'black'})

    ax.set_xticks(x + bar_width * (n_groups - 1) / 2)
    ax.set_xticklabels(genes, rotation=90)
    ax.set_ylabel('Duplicate Frequency')
    ax.set_xlabel('Gene')
    ax.grid(True)
    ax.legend(title='Tissue', loc='upper left', bbox_to_anchor=(1, 1))

    if title is None:
        title = 'TCR Gene Usage: ' + ' vs '.join(conditions)
    ax.set_title(title)

    plt.tight_layout()
    return fig, ax

def plot_junction_aa_length(
        df: pd.DataFrame,
        x_col:str='CDR3_LENGTH',
        y_col:str='CDR3_COUNT',
        hue_col:str='condition',
        palette:str='Set1',
        figsize:tuple[int,int]=(18, 8),
        log_scale:bool=False
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
    Returns
    -------
    fig, ax : tuple
        - Matplotlib figure and axes objects
    """

    fig, ax = plt.subplots(figsize=figsize)
    sns.barplot(data=df, x=x_col, y=y_col, hue=hue_col, palette=palette, ax=ax)
    if log_scale:
        ax.set_yscale('log')
    plt.xticks(rotation=90)
    plt.tight_layout()
    return fig, ax

def plot_vj_chord_diagram(
        df:pd.DataFrame,
        freq_col:str="duplicate_frequency_avg",
        ticks_interval:float=0.03,
        space:int=2,
        cmap:str='Set1',
        figsize:tuple=(12, 12),
        max_links:int=None
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
    
    Returns
    -------
    fig, ax : tuple[Figure, Axes]
        - Matplotlib Figure and Axes objects.
    """
    # Optionally reduce overcrowding
    if max_links and len(df) > max_links: ## Limit to top N V–J links
        df = df.nlargest(max_links, freq_col)
    nodes = sorted(set(df["v_level"]).union(set(df["j_level"])))
    # Build matrix
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
        ticks_kws = dict(
            label_size=8,
            label_orientation='vertical')
    )

    fig = circos.plotfig(figsize=figsize)
    ax = fig.gca()  # get current axes (should be polar)
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
        title:str=None
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
    ylim = np.max(values)

    # Calculate angles and width
    num_angles = len(values) + pad * len(unique_groups)
    angles = np.linspace(0, 2 * np.pi, num=num_angles, endpoint=False)
    width = (2 * np.pi) / len(angles)

    # Choose a seaborn palette
    my_palette = sns.color_palette(pallete, len(unique_groups))
    # Assign colors to groups based on the palette
    group_to_color = dict(zip(unique_groups, my_palette))
    colors = [group_to_color[group] for group in groups]

    # Generate points from start to ylim, with num_intervals steps
    ref_pts = np.linspace(start, ylim, num_intervals)
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
    ax.set_ylim(-ylim, ylim)
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
    for group, size in zip(unique_groups, groups_sizes):
        #for each group minimum and maximum id
        min_idx, max_idx = df[df.group == group].place.min(), df[df.group == group].place.max()

        # Add line below bars
        x1 = np.linspace(angles[offset + pad], angles[offset + size + pad - 1], num=50)
        ax.plot(x1, [-0.07] * 50, color="black", linewidth = 2)

        # Add text to indicate group
        ax.text(
            np.mean(x1),
            -0.15,
            group,
            color="black",
            fontsize=10,
            fontweight="bold",
            ha="center",
            va="center"
        )

        delta = 0.035
        # Min index label slightly left and below first bar
        ax.text(
            x1[0],
            -delta,
            str(min_idx),
            size=9,
            alpha=1,
            ha='center',
            va='top'
        )
        # Max index label slightly right and below last bar
        ax.text(
            x1[-1],
            -delta,
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
                rp[0]-0.02,
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
        ax.set_title(title)
    plt.tight_layout()
    return fig, ax
