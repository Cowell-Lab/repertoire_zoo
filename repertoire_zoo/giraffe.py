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

# Gene Usage
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
    title : str | dict, optional
        Optional title for the axis. A dictionary can be passed to set other **kwargs for `ax.set_title()`.
        (optional, Default: `None`)

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
    
    if title:
        if isinstance(title, dict):
            ax.set_title(title.pop('t'), **title)
        else:
            ax.set_title(title)

    plt.tight_layout()
    return fig, ax

def plot_duplicate_frequency(
        df:pd.DataFrame,
        errorbar:str=None,
        figsize:tuple=(16,8),
        title:str|dict=None,
        hue=None,
        hue_order=None,
        palette:dict=None,
        threshold:float=None,
        ylim:tuple[float,float]=(None,None),
        ax:Axes=None,
        split:bool=False,
        categories:list=None,
        legend:bool=False
    ) -> tuple[Figure, Axes]:
    """
    Plot grouped bar chart with error bars from combined dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        - A Pandas DataFrame with columns `gene`, `N`, `duplicate_frequency_avg`,
        `duplicate_frequency_std`, `condition`
    errorbar : str
        - Display the error bar using defined method: `sd`, `se`. `sd` is for the standard deviation
         while `se` reports the standard error. (Default: None.)
    palette : list
        - A list of colors for each condition. (optional, Default: `None`.)
        uses a predefined color palette. This currently supports up to 5 groups)
    figsize : tuple
        - The figure size. (Default:`(16,8)`)
    title : str | dict, optional
        Optional title for the axis. A dictionary can be passed to set other **kwargs for `ax.set_title()`.
        (optional, Default: `None`)
    threshold : float
        - Define a duplicate frequency threshold to apear. Values are inclusive. (Default: `None`.)
    ylim : tuple[float, float]
        - Define the y-axim limit. If a tuple is provided, then read as [min, max]. 
        It is also possible to define the tuple using `None` to denote automatic limits. 
        (Default: `(None, None)`.)
    ax : Axes
        - A `matplotlib.axes.Axes` object to plot with. (Deafult: `None`.)
    categories : list
        - A list of strings of the names of categories to plot on x-axis

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
    if not split:
        hue_order=None

    # reduce df by threshold

    ## will keep the group if at least one value in 'gene' is not 0
    # if threshold is not None and isinstance(threshold, float):
    #     df.loc[:, 'duplicate_frequency_avg'] = df['duplicate_frequency_avg'].where(df['duplicate_frequency_avg'] >= threshold, 0)
    #     df = df.groupby('gene').filter(lambda g: not (g['duplicate_frequency_avg'] == 0).all())
    
    # will drop the group if all values of 'gene' in group are below threshold
    if threshold is not None and isinstance(threshold, float):
        df = df.groupby('gene').filter(
            lambda g: (g['duplicate_frequency'] >= threshold).any()
        )

    if categories is not None:
        df = df[df['gene'].isin(categories)].copy()
        # df.loc[:, 'gene'] = df['gene'].where(df['gene'].isin(categories))
        placeholder_rows = []
        for cat in categories:
            if cat not in df['gene'].unique():
                print(f"Warning: category {cat} has no data")
                row = {'gene':cat}
                row['duplicate_frequency_avg'] = 0
                row['condition'] = hue_order[0]
                placeholder_rows.append(row)
        df_placeholders = pd.DataFrame(placeholder_rows)
        df = pd.concat([df, df_placeholders], ignore_index=True)
                
    if not split:
        sns.barplot(
            data=df,
            x='gene',
            y='duplicate_frequency_avg',
            hue=hue,
            palette=palette,
            # hue_order=hue_order,
            order=categories if categories is not None else None,
            capsize=0.3,
            ax=ax,
            legend=legend,
            errorbar=None
        ) 
    else:
        sns.barplot(
            data=df,
            x='gene',
            y='duplicate_frequency_avg',
            hue='condition',
            ax=ax,
            legend=legend,
            hue_order=hue_order,
            palette=palette,
            order=categories if categories is not None else None
        )

    # errorbars
    def bar_width_from_data_coord_to_pt(ax, x, width):
        """
        Convert the width of a bar from data coordinates to points.

        Parameters
        ----------
        ax : Axes
        bar : 

        Returns
        -------
        width_pt : float
            - The width of the passed bar in points.
        """
        # data -> display coord sys
        x0, _ = ax.transData.transform((x, 0))
        x1, _ = ax.transData.transform((x + width, 0))
        width_px = x1 - x0

        # px -> in -> pt
        width_in = width_px/ax.figure.dpi
        width_pt = width_in*72

        return width_pt

    if errorbar is not None:
        total_width = 0.8  # typical bar width
        n_hues = len(hue_order) if hue_order else 1
        bar_width = total_width / n_hues

        x_order = [e.get_text() for e in ax.get_xticklabels()]
        for _, row in df.iterrows():
            x_index = x_order.index(row['gene'])
            
            if hue_order:
                hue_index = hue_order.index(row['condition'])
                x_center = x_index - total_width/2 + bar_width/2 + hue_index*bar_width
            else:
                x_center = x_index

            if errorbar=='sd':
                ax.errorbar(
                    x=x_center,
                    y=row['duplicate_frequency_avg'],
                    yerr=row['duplicate_frequency_std'],
                    fmt='none',
                    ecolor='black',
                    capsize=0.25*bar_width_from_data_coord_to_pt(ax, x_center - bar_width/2, bar_width)
                )

            if errorbar=='se':
                ax.errorbar(
                    x=x_center,
                    y=row['duplicate_frequency_avg'],
                    yerr=row['duplicate_frequency_std']/np.sqrt(row['N']),
                    fmt='none',
                    ecolor='black',
                    capsize=0.25*bar_width_from_data_coord_to_pt(ax, x_center - bar_width/2, bar_width)
                )
            
    ax.tick_params(axis='x', rotation=90)

    if title:
        if isinstance(title, dict):
            ax.set_title(title.pop('t'), **title)
        else:
            ax.set_title(title)
    
    if ylim[0] is not None or ylim[1] is not None:
        ax.set_ylim(bottom=ylim[0], top=ylim[1])

    return fig, ax

def plot_ratio(
        df:pd.DataFrame,
        region:str='mu_ratio_fwr',
        errorbar:str=None,
        figsize:tuple=(16,8),
        title:str|dict=None,
        hue=None,
        hue_order=None,
        palette:dict=None,
        threshold:float=None,
        ylim:tuple[float,float]=(None,None),
        ax:Axes=None,
        split:bool=False,
        categories:list=None,
        legend:bool=False
    ) -> tuple[Figure, Axes]:
    """
    Plot grouped bar chart with error bars from combined dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        - A Pandas DataFrame with columns `gene`, `N`, `duplicate_frequency_avg`,
        `duplicate_frequency_std`, `condition`
    errorbar : str
        - Display the error bar using defined method: `sd`, `se`. `sd` is for the standard deviation
         while `se` reports the standard error. (Default: None.)
    palette : list
        - A list of colors for each condition. (optional, Default: `None`.)
        uses a predefined color palette. This currently supports up to 5 groups)
    figsize : tuple
        - The figure size. (Default:`(16,8)`)
    title : str | dict, optional
        Optional title for the axis. A dictionary can be passed to set other **kwargs for `ax.set_title()`.
        (optional, Default: `None`)
    threshold : float
        - Define a duplicate frequency threshold to apear. Values are inclusive. (Default: `None`.)
    ylim : tuple[float, float]
        - Define the y-axim limit. If a tuple is provided, then read as [min, max]. 
        It is also possible to define the tuple using `None` to denote automatic limits. 
        (Default: `(None, None)`.)
    ax : Axes
        - A `matplotlib.axes.Axes` object to plot with. (Deafult: `None`.)
    categories : list
        - A list of strings of the names of categories to plot on x-axis

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
        hue='repertoire_group_name'
    if not split:
        hue_order=None

    # reduce df by threshold

    ## will keep the group if at least one value in 'gene' is not 0
    # if threshold is not None and isinstance(threshold, float):
    #     df.loc[:, 'duplicate_frequency_avg'] = df['duplicate_frequency_avg'].where(df['duplicate_frequency_avg'] >= threshold, 0)
    #     df = df.groupby('gene').filter(lambda g: not (g['duplicate_frequency_avg'] == 0).all())
    
    # will drop the group if all values of 'gene' in group are below threshold
    # if threshold is not None and isinstance(threshold, float):
    #     df = df.groupby('gene').filter(
    #         lambda g: (g['duplicate_frequency'] >= threshold).any()
    #     )

    # if categories is not None:
    #     df = df[df['gene'].isin(categories)].copy()
    #     # df.loc[:, 'gene'] = df['gene'].where(df['gene'].isin(categories))
    #     placeholder_rows = []
    #     for cat in categories:
    #         if cat not in df['gene'].unique():
    #             print(f"Warning: category {cat} has no data")
    #             row = {'gene':cat}
    #             row['duplicate_frequency_avg'] = 0
    #             row['condition'] = hue_order[0]
    #             placeholder_rows.append(row)
    #     df_placeholders = pd.DataFrame(placeholder_rows)
    #     df = pd.concat([df, df_placeholders], ignore_index=True)
                
    if not split:
        sns.barplot(
            data=df,
            x='repertoire_group_name',
            y=region+'_avg',
            hue=hue,
            palette=palette,
            # hue_order=hue_order,
            order=categories if categories is not None else None,
            capsize=0.3,
            ax=ax,
            legend=legend,
            errorbar=None
        ) 
    else:
        sns.barplot(
            data=df,
            x='repertoire_group_name',
            y=region+'_avg',
            hue='repertoire_group_name',
            ax=ax,
            legend=legend,
            hue_order=hue_order,
            order=categories if categories is not None else None
        )

    # errorbars
    def bar_width_from_data_coord_to_pt(ax, x, width):
        """
        Convert the width of a bar from data coordinates to points.

        Parameters
        ----------
        ax : Axes
        bar : 

        Returns
        -------
        width_pt : float
            - The width of the passed bar in points.
        """
        # data -> display coord sys
        x0, _ = ax.transData.transform((x, 0))
        x1, _ = ax.transData.transform((x + width, 0))
        width_px = x1 - x0

        # px -> in -> pt
        width_in = width_px/ax.figure.dpi
        width_pt = width_in*72

        return width_pt

    if errorbar is not None:
        total_width = 0.8  # typical bar width
        n_hues = len(hue_order) if hue_order else 1
        bar_width = total_width / n_hues

        x_order = [e.get_text() for e in ax.get_xticklabels()]
        for _, row in df.iterrows():
            x_index = x_order.index(row['repertoire_group_name'])
            
            if hue_order:
                hue_index = hue_order.index(row['repertoire_group_name'])
                x_center = x_index - total_width/2 + bar_width/2 + hue_index*bar_width
            else:
                x_center = x_index

            if errorbar=='sd':
                ax.errorbar(
                    x=x_center,
                    y=row[region+'_avg'],
                    yerr=row[region+'_std'],
                    fmt='none',
                    ecolor='black',
                    capsize=0.25*bar_width_from_data_coord_to_pt(ax, x_center - bar_width/2, bar_width)
                )

            if errorbar=='se':
                ax.errorbar(
                    x=x_center,
                    y=row[region+'_avg'],
                    yerr=row[region+'_std']/np.sqrt(row[region+'_N']),
                    fmt='none',
                    ecolor='black',
                    capsize=0.25*bar_width_from_data_coord_to_pt(ax, x_center - bar_width/2, bar_width)
                )
            
    ax.tick_params(axis='x', rotation=90)

    if title:
        if isinstance(title, dict):
            ax.set_title(title.pop('t'), **title)
        else:
            ax.set_title(title)
    
    if ylim[0] is not None or ylim[1] is not None:
        ax.set_ylim(bottom=ylim[0], top=ylim[1])

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
        if isinstance(title, str):
            circos.text(
                title,
                size=text_size,
                deg=0,
                r=140
            )
        else:
            circos.text(
                title['t'],
                size=title.get('fontsize', 16),
                deg=title.get('deg', 0),
                r=title.get('r', 140)
            )

    # if not external_ax:
    #     plt.tight_layout()

    circos.plotfig(figsize=figsize, ax=ax)

    return fig, ax

# Mutations
def plot_mutational_hedgehog(
        df:pd.DataFrame,
        pad:int=5,
        num_intervals:int=4,
        start:float=0.05,
        pallete:str='tab20',
        figsize:tuple[int, int]=(10, 10),
        title:str=None,
        y_max:int=None,
        ax:Axes=None
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
    title : str | dict, optional
        Optional title for the axis. A dictionary can be passed to set other **kwargs for `ax.set_title()`.
        (optional, Default: `None`)
    y_max : int
        - The maximum y value displayed. Can be used to set plots to the same scale. (Default: `None`.)
    ax : Axes
        - A Matplotlib.pyplot.Axes object to plot with. If `None` a new Figure and Axes 
        will be created to plot on to. (Default: `None`.)

    Returns
    -------
    fig, ax : tuple[Figure, Axes]
        - Matplotlib Figure and Axes objects.
    """
    # Set up or call in Figure
    external_ax = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(polar=True))
    else:
        fig = ax.get_figure()
        external_ax = True
    
    values = df["value"].values
    regions = df["region"].values

    # Calculate max y value
    if y_max is None:
        y_max = np.max(values)
    y_min = -0.1

    #Keep the group orientation as is.
    unique_regions = pd.unique(regions)
    regions_sizes = [len(group_df) for _, group_df in df.groupby("region", sort=False)]

    # Calculate angles and width
    num_angles = len(values) + pad * len(unique_regions)
    angles = np.linspace(0, 2 * np.pi, num=num_angles, endpoint=False)
    width = (2 * np.pi) / len(angles)

    # Choose a seaborn palette
    my_palette = sns.color_palette(pallete, len(unique_regions))
    # Assign colors to groups based on the palette
    group_to_color = dict(zip(unique_regions, my_palette))
    colors = [group_to_color[region] for region in regions]

    # Generate points from start to y_max, with num_intervals steps
    ref_pts = np.linspace(start, y_max, num_intervals)
    ref_pts_rounded = np.round(ref_pts, 2)
    ref_pt = [[float(x)] for x in ref_pts_rounded]

    #start plotting
    offset = 0
    idxs = []
    for size in regions_sizes:
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
    for region, size in zip(unique_regions, regions_sizes):
        #for each group minimum and maximum id
        min_idx, max_idx = df[df.region == region].place.min(), df[df.region == region].place.max()

        # Add line below bars
        x1 = np.linspace(angles[offset + pad], angles[offset + size + pad - 1], num=50)
        ax.plot(x1, [-delta] * 50, color="black", linewidth = 2)   

        # Add text to indicate group
        ax.text(
            np.mean(x1),
            -3*delta,
            region,
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
        # x2 = np.linspace(angles[offset], angles[offset + pad - 1], num=50)
        # for rp in ref_pt:
        #     ax.plot(
        #         x2,
        #         rp*50,
        #         color="gray",
        #         lw=1
        #     )
        #     ax.text(
        #         0.09,
        #         rp[0]-delta,
        #         f'{rp[0]}',
        #         color="black",
        #         ha='center',
        #         va='top',
        #         fontsize=9
        #     )
        offset += size + pad

    # ax.set_frame_on(False)
    # ax.xaxis.grid(False)
    # ax.yaxis.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    
    if title:
        if isinstance(title, dict):
            ax.set_title(title.pop('t'), **title)
        else:
            ax.set_title(title)
    
    # if external_ax:
    #     plt.tight_layout()
    
    return fig, ax

# Diversity
def plot_diversity_curve(
        df:pd.DataFrame,
        ax:Axes=None,
        figsize:tuple[int, int]=(16, 8),
        palette:str|dict='Set2',
        hue='condition',
        legend:bool|str='auto',
        legend_title:str=None,
        title:str=None,
        errorbar=None
    ) -> tuple[Figure, Axes]:
    """
    
    Parameters
    ----------
    df : pd.DataFrame
        - A pandas dataframe containing the diversity curve data
    ax : Axes
        - A Matplotlib.pyplot.Axes object to plot with. If `None` a new Figure and Axes 
        will be created to plot on to. (Default: `None`.)
    figsize : tuple[int, int]
        - The size of the figure. (Default: `(16,8)`.)
    palette : str|dict
        - The name of the palette to use, or a custom defined dictionary. (Default: `'Set2'`)
    hue : str
        - The columns to set the hue on. Can be set as 'repertoire_id' or 'condition'. (Default: `'condition'`)
    legend : bool|str
        - Display a legend. Use seaborn documentation. Can use “auto”, “brief”, “full”, or False. (Default: `'auto'`)
    legend_title : str
        - The title of the legend. (Default: `None`.)
    title : str | dict, optional
        Optional title for the axis. A dictionary can be passed to set other **kwargs for `ax.set_title()`.
        (optional, Default: `None`)

    Returns
    -------
    fig, ax : tuple[Figure, Axes]
        - A tuple containing a Figure and Axes object from Matplotlib.pyplot
    """
    # Set up or bring in Figure
    external_ax = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()
        external_ax = True

    # Draw lineplot
    sns.lineplot(
        data=df,
        x='q',
        y='d',
        markers=True,
        hue=hue,
        palette=palette,
        legend=legend,
        ax=ax,
        errorbar=errorbar
    )

    # Set log scale
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Fill between for each group
    # for condition, group_df in df.groupby('condition'):
    #     ax.fill_between(group_df['q'], group_df['d_lower'], group_df['d_upper'], alpha=0.3)

    # Legend and formatting
    if legend != False:
        ax.legend(title=legend_title, bbox_to_anchor=(1, 1))
    
    if title:
        if isinstance(title, dict):
            ax.set_title(title.pop('t'), **title)
        else:
            ax.set_title(title)

    # Only call tight_layout if fig was created
    if not external_ax:
        plt.tight_layout()

    return fig, ax

def plot_diversity_boxplot(
        df:pd.DataFrame,
        ax:Axes=None,
        figsize:tuple[int, int]=(16, 8),
        palette:str|dict='Set2',
        hue=None,
        legend:bool|str='auto',
        legend_title:str=None,
        title:str=None
    ) -> tuple[Figure, Axes]:
    """
    
    Parameters
    ----------
    df : pd.DataFrame
        - A pandas dataframe containing the diversity curve data
    ax : Axes
        - A Matplotlib.pyplot.Axes object to plot with. If `None` a new Figure and Axes 
        will be created to plot on to. (Default: `None`.)
    figsize : tuple[int, int]
        - The size of the figure. (Default: `(16,8)`.)
    palette : str|dict
        - The name of the palette to use, or a custom defined dictionary. (Default: `'Set2'`)
    hue : str
        - The columns to set the hue on. Can be set as 'repertoire_id', 'condition', or None. (Default: `None.`)
    legend : bool|str
        - Display a legend. Use seaborn documentation. Can use “auto”, “brief”, “full”, or False. (Default: `'auto'`)
    legend_title : str
        - The title of the legend. (Default: `None`.)
    title : str | dict, optional
        Optional title for the axis. A dictionary can be passed to set other **kwargs for `ax.set_title()`.
        (optional, Default: `None`)

    Returns
    -------
    fig, ax : tuple[Figure, Axes]
        - A tuple containing a Figure and Axes object from Matplotlib.pyplot
    """
    # Set up or bring in Figure
    external_ax = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()
        external_ax = True

    df_q1 = df[df['q']==1].reset_index(drop=True)
    df_q0 = df[df['q']==0].reset_index(drop=True)
    assert df_q0.shape[0]==df_q1.shape[0], 'Unequal quantity of values for q=0 and q=1'

    df_new = pd.DataFrame(columns=['repertoire_id', 'H_norm', 'condition'])
    for i in range(len(df_q0)):
        df_new.loc[i,'repertoire_id'] = df_q0.loc[i,'repertoire_id']
        # df_new.loc[i,'H_norm'] = df_q1[df_q1['repertoire_id']==df_q0.loc[i,'repertoire_id']].loc[i,'d']/df_q0.loc[i,'d']
        df_new.loc[i,'H_norm'] = np.log2(df_q1[df_q1['repertoire_id']==df_q0.loc[i,'repertoire_id']].loc[i,'d'])/np.log2(df_q0.loc[i,'d'])
        df_new.loc[i,'condition'] = df_q0.loc[i,'condition']
        
    # Draw lineplot
    sns.boxplot(
        data=df_new,
        x=hue,
        y='H_norm'
    )

    # Set log scale
    # ax.set_yscale('log')
    ax.tick_params('x', rotation=45)

    # Legend and formatting
    if legend != False:
        ax.legend(title=legend_title, bbox_to_anchor=(1, 1))
    
    if title:
        if isinstance(title, dict):
            ax.set_title(title.pop('t'), **title)
        else:
            ax.set_title(title)

    # Only call tight_layout if fig was created
    if not external_ax:
        plt.tight_layout()

    return fig, ax

# Clonal Abundance
def plot_clonal_abundance(
        df:pd.DataFrame,
        y_col:str='cumulative_abundance',
        ax:Axes=None,
        figsize:tuple[int,int]=(16, 8),
        title:str=None,
        legend_title:str=None,
        hue=None,
        errorbar=None
    ) -> tuple[Figure, Axes]:
    """
    Plot the clonal abundance

    Parameters
    ----------
    df : pd.DataFrame
        - The Pandas DataFrame with the relevant information.
    y_col : str
        - The name of the y-axis column to use. (Default: `'cumulative_abundance'`.)
    ax : Axes
        - A matplotlib.pyplot.Axes object to plot on to. If an Axes object is not provided
        a new Figure and Axes will be generated. (Default: `None`.)
    figsize : tuple[int, int]
        - The size of the figure. (Default: `(16, 8)`.)
    title : str | dict, optional
        Optional title for the axis. A dictionary can be passed to set other **kwargs for `ax.set_title()`.
        (optional, Default: `None`)
    legend_title : str
        - The title of the legend. (Default: `None`.)
    hue : str
        - The column in `df` to be used to split for groups. (Default: `None`.)
    errorbar : str
        - The type of errorbar to use as currently defined by seaborn lineplot. (Default: `None`.)

    Returns
    -------
    fig, ax : tuple[Figure, Axes]
        - A tuple containing a Figure and Axes object from Matplotlib.pyplot
    """
    # Set up or bring in Figure
    external_ax = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()
        external_ax = True

    sns.lineplot(data=df, x='rank', y=y_col, hue=hue, ax=ax, errorbar=errorbar)
    
    ax.set_xscale('log')

    # Legend and formatting
    if legend_title:
        ax.legend(title=legend_title, loc='upper center', bbox_to_anchor=(0.5,-0.1), ncol=1)
    
    if title:
        if isinstance(title, dict):
            ax.set_title(title.pop('t'), **title)
        else:
            ax.set_title(title)
    
    # Only call tight_layout if new fig was created
    if not external_ax:
        plt.tight_layout()
    
    return fig, ax

# CDR3
def plot_junction_aa_length(
        df: pd.DataFrame,
        x_col:str='CDR3_LENGTH',
        y_col:str='CDR3_COUNT',
        controls:list=[None],
        title:str|dict=None,
        hue_col:str='condition',
        palette:str='Set1',
        legend:str=None,
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
    controls : list
        - A list of control group names. This will display a baseline on the plot. If no controls are needed, 
        then the default value of `[None]` is used.
    title : str | dict, optional
        Optional title for the axis. A dictionary can be passed to set other **kwargs for `ax.set_title()`.
        (optional, Default: `None`)
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
    
    sns.barplot(
        data=df,
        x=x_col,
        y=y_col,
        hue=hue_col,
        palette=palette,
        ax=ax,
        legend=legend
    )

    if controls != [None]:
        for control in controls:
            sns.lineplot(
                df[df['condition']==control],
                x=x_col,
                y=y_col,
                label=control,
                color=palette[control],
                ls='-',
                lw=2,
                ax=ax
            )

    if log_scale:
        ax.set_yscale('log')
        
    # ax.tick_params(axis='x', labelrotation=45)
    ax.set_xlim(aa_range[0], aa_range[1])

    if title:
        if isinstance(title, dict):
            ax.set_title(title.pop('t'), **title)
        else:
            ax.set_title(title)

    if not external_ax:
        fig.tight_layout()

    return fig, ax

def plot_cdr3_aa_sharing_heatmap(
        df:pd.DataFrame,
        cmap:str|dict='viridis',
        figsize=(20, 15),
        title:str|dict=None,
        ax:Axes=None
    ) -> tuple[Figure, Axes]:
    """
    Plot heatmap from preprocessed matrix and labels.

    Parameters
    ----------
    df : pd.DataFrame
        - The Pandas DataFrame that contains the data to plot.
    cmap : str|dict
        - The colormap to use for plotting. Can either be a string that aligns
        with a seaborn colormap, or a custom defined dictionary. (Default: `'viridis'`.)
    figsize : tuple[int, int]
        - The size of the figure. (Default: `(20,15)`.)
    title : str | dict, optional
        Optional title for the axis. A dictionary can be passed to set other **kwargs for `ax.set_title()`.
        (optional, Default: `None`)
    ax : Axes
        - A matplotlib.pyplot.Axes object to plot on to. If an Axes object is not provided
        a new Figure and Axes will be generated. (Default: `None`.)

    Returns
    -------
    fig, ax : tuple[Figure, Axes]
        - A tuple containing a Figure and Axes object from Matplotlib.pyplot
    """
    # Set up or bring in Figure
    external_ax = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()
        external_ax = True
    
    if df.select_dtypes(include=['float']).empty:
        fmt = 'd'  # no floats found, assume integers
    else:
        fmt = '.2f'
    
    sns.heatmap(df, annot=True, fmt=fmt, cmap = cmap, ax=ax)
    
    ax.tick_params(axis = 'x', rotation=90)
    
    if title:
        if isinstance(title, dict):
            ax.set_title(title.pop('t'), **title)
        else:
            ax.set_title(title)

    # Only call tight_layout if new fig was created
    if not external_ax:
        plt.tight_layout()

    return fig, ax

def plot_cdr3_sequence_count(
        df:pd.DataFrame,
        x_col:str='junction_aa',
        abundance:bool=False,
        top_n:int=None,
        figsize:tuple[int,int]=(16,8),
        title:str|dict=None,
        ax:Axes=None
    ) ->tuple[Figure, Axes]:
    """
    Plot a barplot of top CDR3 sequences.
    
    Parameters
    ----------
    df : pd.DataFrame
        - DataFrame of top CDR3 sequences
    x_col : str
        - Column for x-axis (junction_aa for cdr3)
    abundance : bool
        - Whether to use duplicate_count (True) or sequence_count (False)
    top_n : int
        - Number of top sequences to plot. If no value is provided, then all available sequences will be plotted.
        (Default: `None`.)
    figsize : tuple[int, int]
        - Size of the Figure
    title : str | dict, optional
        Optional title for the axis. A dictionary can be passed to set other **kwargs for `ax.set_title()`.
        (optional, Default: `None`)
    ax : Axes
        - A matplotlib.pyplot.Axes object to plot on to. If an Axes object is not provided
        a new Figure and Axes will be generated. (Default: `None`.)

    Returns
    -------
    fig, ax : tuple[Figure, Axes]
        - A tuple containing a Figure and Axes object from Matplotlib.pyplot
    """
    # Set up or bring in Figure
    external_ax = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()
        external_ax = True
    
    y_col = 'duplicate_count' if abundance else 'sequence_count'

    if top_n is None:
        top_n = df.shape[0]
    
    df = df.sort_values(by=y_col, ascending=False).iloc[:top_n]
    
    sns.barplot(data=df, x=x_col, y=y_col, ax=ax)
    ax.tick_params(axis='x', rotation=90)

    if title:
        if isinstance(title, dict):
            ax.set_title(title.pop('t'), **title)
        else:
            ax.set_title(title)

    if not external_ax:
        plt.tight_layout()
    
    return fig, ax

def plot_unique_number_of_cdr3(
        df,
        palette='Set2',
        figsize=(16, 6)
    ):
    
    fig, (ax1, ax2) = plt.subplots(
        1,
        ncols=2,
        figsize=figsize,
        sharey=True,
        gridspec_kw={'width_ratios': [4, 1]}
    )
    
    sns.barplot(
        data=df,
        x='sample_id',
        y='n_unique_cdr3',
        hue='condition',
        palette=palette,
        ax=ax1,
        legend=False
    )
    
    sns.boxplot(
        data=df,
        x='condition',
        y='n_unique_cdr3',
        hue='condition',
        showmeans=True,
        palette=palette,
        ax=ax2,
        legend=True
    )

    sns.stripplot(
        data=df,
        x='condition',
        y='n_unique_cdr3',
        color='black',
        dodge=True,
        ax=ax2
    )

    ax1.tick_params(axis = 'x', rotation = 90)
    ax2.tick_params(axis = 'x', rotation = 90)
    ax2.legend(title = 'Condition', bbox_to_anchor = (1, 1))
    
    return fig, (ax1, ax2)


