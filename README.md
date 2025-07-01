# Repertoire Zoo

A set of python modules for repertoire analysis.

# Table of Contents

- [Installation](#installation)
- [Modules](#modules)
  - [bird: Bring in Repertoire Data](#bird--bring-in-repertoire-data)
  - [giraffe: Graphical Interactive Repertoire Analysis with Fast Feature Exploration)](#giraffe--graphical-interactive-repertoire-analysis-with-fast-feature-exploration)
- [Tests](#tests)
- [License](#license)

# Installation

# Modules

The modules are contained within the `repertoire_zoo` subdirectory folder.

### `bird` : **B**ring **I**n **R**epertoire **D**ata

- Contains functions that will provide Pandas DataFrames containing various type of information from the `repertoires.airr.json` file that can be exported from VDJ Server. Can be called by `import repertoire_zoo.bird as br`. Functions include:
  - `read_sample_info_airr(path)` : Creates and returns a Pandas DataFrame containing sample information from the repertoire file.
    - `path` : str. Path to `repertoires.airr.json` or alternative `*.airr.json` containin the repertoire information.

  - `read_repertoire_group_info_airr(path, silence)` : Creates a Pandas DataFrame containing sample information from the repertoire file.
    - `path` : str. Path to `repertoires.airr.json` or alternative `*.airr.json` containin the repertoire information.
    - `silence` : bool. When `False`, prints the total number of repertoire groups. (Default: `False`)

  - `load_and_prepare_data_groups(repcalc_dir, groups, call_type, level, mode, processing_stage)` : Load and preprocess data from files for multiple groups. Returns a Pandas DataFrame with columns: `gene`, `duplicate_frequency_avg`, `duplicate_frequency_std`, `condition`.
    - `repcalc_dir` : str. RepCalc base directory path.
    - `groups` : list. List of tuples (group_id, condition_name)
    - `call_type` : str. The type of call. Can be `v_call`, `d_call`, `j_call`. Default: `v_call`.
    - `level` : str. Type of gene level information. Can be `allele`, `gene`, `subgroup`. Default: `gene`.
    - `mode` : str. Can be `exists`, `proportion`, or `unique`. Default: `proportion`.
    - `processing_stage` : str. The middle part of the filename path

  - `load_and_prepare_data_subjects(repcalc_dir, repertoire_groups, group_id, sample_info_df, processing_stage, level, mode, call_type)` : Load and preprocess data from files one group with multiple repertoires. Returns a Pandas DataFrame with columns: `gene`, `duplicate_frequency_avg`, `duplicate_frequency_std`, `condition`.
    - `repcalc_dir` : str. The base directory path
    - `repertoire_groups` : dict. All repertoire groups
    - `group_id` : str. The group we want to look at
    - `sample_info_df` : pd.DataFrame. Pandas DataFrame containing metadata about samples
    - `processing_stage` : str. The middle part of the filename path. (Default: `.igblast.makedb.allele.clone.`.)
    - `level` : str. The type of gene level information. Can be `allele`, `gene`, `subgroup`. (Default: `subgroup`.)
    - `mode` : str. Can be `exists`, `proportion`, `unique`. (Default: `proportion`.)
    - `call_type` : str. The type of call. Can be `v_call`, `d_call`, `j_call`. (Default: `v_call`.)

  - `load_and_prepare_data_junction_aa(repcalc_dir, groups, processing_stage, call_type)` : Load and preprocess data for junction amino acid length. Returns Pandas DataFrame with columns: `CDR3_LENGTH`, `CDR3_COUNT`, `CDR3_RELATIVE`.
    - `repcalc_dir` : str. The base directory path.
    - `groups` : list. A list of tuples in the format of [(`group_id` , `name`)]
    - `processing_stage` : str. The middle part of the filename path. Default: `.igblast.makedb.allele.clone.`
    - `call_type` : str. The type of call. (Default: `junction_aa_length`.)

  - `load_and_prepare_data_vj_combo(repcalc_dir, group_id, combo_type, level, mode, productive, processing_stage)` : Load and prepare data to be used in the VJ Combo Circos Plot. Returns a Pandas DataFrame containings the coloumns: `v_level`, `j_level`, `duplicate_frequency_avg`, `duplicate_frequency_std`.
    - `repcalc_dir` : str. The base directory path.
    - `group_id` : str. The group ID.
    - `combo_type` : str. The VDJ combination. (Default `vj_combo`.)
    - `level` : str. The type of gene level information. Can be `allele`, `gene`, `subgroup`. (Default: `subgroup|subgroup`.)
    - `mode` : str. Can be `exists`, `proportion`, `unique`. Default: `proportion`.
    – `productive` : str. Filter for productive contigs. Can be `TRUE` or `FALSE`. (Default: `TRUE`.)
    - `processing_stage` : str. The middle part of the filename path. Default: `.igblast.makedb.allele.clone.group.`

  - `load_and_prepare_mutation_data(mutation_data_dir, group_name, trim_codons, processing_stage)` : Load and prepare mutation data for hedgehog plots. Returns a Pandas DataFrame containing the columns: `place`, `group`, `value`, `se`, `se_start`, `se_end`.
    - `mutation_data_dir` : str. Path to mutation data directory.
    - `group_name` : str. The name of the group.
    - `trim_codons` : int. The number of codons to trim. (Default: `16`.)
    - `processing_stage` : str. The current processing stage. Used for the file_path. (Default: `gene.mutations`.)

### `giraffe` : **G**raphical **I**nteractive **R**epertoire **A**nalysis with **F**ast **F**eature **E**xploration

- Contains functions that allow for imported repertoire and repertoire goup data to be plotted. Can be called by `import repertoire_zoo.giraffe as gr`. Functions include:

  - `plot_duplicate_frequency_average(df_combined, colors, figsize, title)` : Plot grouped bar chart with error bars from combined dataframe. Returns a tuple containing the matplotlib figure and axes objects.
    - `df_combined` : pd.DataFrame. A pandas dataframe with columns `gene`, `duplicate_frequency_avg`, `duplicate_frequency_std`, `condition`.
    - `colors` : list. A list of colors for each condition. (optional, Default: `None` uses a predefined color palette. This currently supports up to 5 groups)
    - `figsize` : tuple. The figure size. (Default:`(16,8)`)
    - `title` : str. The plot's title (optional, Default: `None`)

  - `plot_duplicate_frequency(df_combined, colors, figsize, title)` : Plot grouped bar chart with error bars from combined dataframe. Returns a tuple containing the matplotlib figure and axes objects.
    - `df_combined` : pd.DataFrame. A pandas dataframe with columns `gene`, `duplicate_frequency_avg`, `duplicate_frequency_std`, `condition`.
    - `colors` : list. A list of colors for each condition. (optional, Default: `None` uses a predefined color palette. This currently supports up to 5 groups)
    - `figsize` : tuple. The figure size. (Default:`(16,8)`)
    - `title` : str. The plot's title (optional, Default: `None`)

  - `plot_junction_aa_length` : Plot grouped bar chart with error bars from combined dataframe. Returns a tuple containing the matplotlib figure and axes objects.
    - `df` : pd.DataFrame. A Pandas DataFrame with columns `gene`, `duplicate_frequency_avg`, `duplicate_frequency_std`, `condition`
    - x_col : str. The x-variable column name. (Default: `CDR3_LENGTH`)
    - y_col : str. The y-variable column name. (Default: `CDR3_COUNT`)
    - hue_col : str. The coloumn name to determine the hue color. (Default: `condition`)
    - palette : str. The seaborn color palette used. (Default: `Set1`)
    - figsize : tuple. The figure size. (Default:`(18,8)`)
    - log_scale : bool. Place the y axis into a log scale. (Default: `False`)

  - `plot_vj_chord_diagram(df, freq_col, ticks_interval, space, cmap, figsize, max_links)` : Plots the VJ Chord Diagram. Returns a tuple containing Matplotlib Figure and Axes objects.
    - `df` : pd.DataFrame. The dataframe from `bird.load_and_prepare_data_vj_combo`.
    - `freq_col` : str. The name of the frequency column. (Default: `duplicate_frequency_avg`.)
    - `ticks_interval` : float. The interval between each tick on the plot. (Default: `0.03`.)
    - `space` : int. The spacing on the diagram. (Default: `2`.)
    - `cmap` : str. The colormap used. (Default: `Set1`)
    - `figsize` : tuple[int, int]. The size of the figure. (Default: `(12, 12)`.)
    - `max_links` : int. The maximum nunber of links. (Default: `None`.)

  - `plot_mutational_hedgehog(df, pad, num_intervals, start, pallete, figsize, title)` : Plot grouped bar chart with error bars from mutational Data. Returns a tuple containing Matplotlib Figure and Axes objects.
    - `df` : pd.DataFrame. DataFrame with columns: `value`, `group`, `place`, `se`, `se_start`, `se_end`.
    - `pad` : int. Spacing between each group. (Default: `5`.)
    - `num_intervals` : int. An integer number for showing the intervals for y limit (Default: `4`.)
    - `start` : float. A float for starting number for the interval (Default: `0.05`.)
    - `pallete` : str. Seaborn color pallete (Default: `tab20`.)
    - `figsize` : tuple[int, int]. The figure size. (Default: `(10, 10)`.)
    - `title` : str. The title of the plot. (Default: `None`.)

# Tests

Currently there are no tests at this time.

# To Do
- Fix "too many local variables" lint error in `giraffe`.
- Fix "too many argumnets/positional arguments" lint error in `giraffe`
- Fix "too many local variables" lint error in `bird`.