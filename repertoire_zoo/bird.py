"""
Module Name: bird

**B**ring **i**n **R**epertoire **D**ata.
This module provides plotting utilities for analyzing repertoire data.

Functions:
- 
- read_sample_info_airr(path) : Creates and returns a Pandas DataFrame 
  containing sample information from the repertoire file.
- read_repertoire_group_info_airr(path, silence) : Creates a Pandas 
  DataFrame containing sample information from the repertoire file.

Usage:
For reading in repertoire and repertoire group data from an `.airr.tsv` file.
Import this module and call the desired functions with appropriate arguments.
"""
import airr
import numpy as np
import pandas as pd

def read_repertoire_info_airr(
        path:str='repertoires.airr.json'
    ) -> pd.DataFrame:
    """
    Creates a Pandas DataFrame containing repertoire information from the AIRR JSON metadata file.

    Parameters
    ----------
    path : str
        - Path to `repertoires.airr.json` or alternative `*.airr.json` containing the
        repertoire information.
    
    Returns
    -------
    df : pd.DataFrame
        - Pandas DataFrame containing sample information
    """
    info_list = []
    data = airr.read_airr(path)
    for rep in data['Repertoire']:
        subject = rep.get("subject", {})
        sample = rep.get("sample", [{}])[0]

        repertoire_id = rep.get("repertoire_id", None)
        subject_id = subject.get('subject_id', None)
        sample_id = sample.get('sample_id', None)
        age = subject.get('age_min', None)
        ethnicity = subject.get('ethnicity', None)
        race = subject.get('race', None)
        diagnosis = subject.get('diagnosis', [{}])[0].get('study_group_description', None)
        disease_state = sample.get('disease_state_sample', None)
        tissue = sample.get('tissue', {}).get('label', None)
        library = sample.get('sequencing_run_id', None)
        info_list.append(
            {
                'repertoire_id':repertoire_id,
                'subject_id':subject_id,
                'sample_id':sample_id, 
                'age':age, 
                'ethnicity':ethnicity, 
                'race':race, 
                'tissue':tissue, 
                'library':library,
                'diagnosis':diagnosis,
                'disease_state':disease_state
            }
        )
    df = pd.DataFrame(info_list)
    df['ss_id'] = df['subject_id'] + '_' + df['sample_id']

    return df

def read_repertoire_group_info_airr(
        path:str='repertoires.airr.json',
        silence:bool=False
    ) -> pd.DataFrame:
    """
    Creates a Pandas DataFrame containing sample information from the repertoire file.

    Parameters
    ----------
    path : str
        - Path to `repertoires.airr.json` or alternative `*.airr.json` containing the
        repertoire groups information.
    silence : bool
        - Prints the total number of repertoire groups if `False`. (Default: `False`)
        
    Returns
    -------
    df : pd.DataFrame
        - Pandas DataFrame containing `repertoire_group_id` and group `name`
    """
    repertoire_group_data = airr.read_airr(path)['RepertoireGroup']
    repertoire_group_list = [
        {
            'repertoire_group_id':obj['repertoire_group_id'], 
            'name':obj['repertoire_group_name'],
            'description':obj.get('description'),
            'repertoires':obj['repertoires'],
            'repertoire_count':len(obj['repertoires'])
        }
        for obj in repertoire_group_data
    ]
    df = pd.DataFrame(repertoire_group_list)

    if not silence:
        print("Total Number of Repertoire Groups: ", len(df))

    return df

def load_gene_usage_data(
        repcalc_dir:str,
        repertoire_id:str,
        processing_stage:str,
        call_type:str='v_call',
        level:str='gene',
        mode:str='proportion',
        productive:bool=True
    ) -> pd.DataFrame:
    """
    Load gene usage data for a single repertoire.

    Parameters
    ----------
    repcalc_dir : str
        - RepCalc base directory path.
    repertoire_id : str
        - Repertoire ID use to determine file name
    processing_stage : str
        - The middle part of the filename path
    call_type : str
        - The type of call. Can be `v_call`, `d_call`, `j_call`. Default: `v_call`.
    level : str
        - Type of gene level information. Can be `allele`, `gene`, `subgroup`.
        Default: `gene`.
    mode : str
        - Can be `exists`, `proportion`, or `unique`. Default: `proportion`.
    productive: bool
        - Only productive rearrangements (True), both productive and non-productive (False)

    Returns
    -------
    df : pd.DataFrame
    """
    path = f"{repcalc_dir}{repertoire_id}.{processing_stage}.{call_type}.tsv"
    df = pd.read_csv(path, sep='\t')
    df = df[(df['level']==level) & (df['mode']==mode) & (df['productive']==productive)]
    return df

def load_gene_usage_group_data(
        repcalc_dir:str,
        groups:list,
        processing_stage:str,
        call_type:str='v_call',
        level:str='gene',
        mode:str='proportion',
        productive:bool=True,
        threshold:float=None
    ) -> pd.DataFrame:
    """
    Load gene usage statistical data for a single group.

    Parameters
    ----------
    repcalc_dir : str
        - RepCalc base directory path.
    groups : list
        - List of tuples (group_id, condition_name)
    call_type : str
        - The type of call. Can be `v_call`, `d_call`, `j_call`. Default: `v_call`.
    level : str
        - Type of gene level information. Can be `allele`, `gene`, `subgroup`.
        Default: `gene`.
    mode : str
        - Can be `exists`, `proportion`, or `unique`. Default: `proportion`.
    processing_stage : str
        - The middle part of the filename path
    threshold : float
        - Define a duplicate frequency threshold to apear. Values are inclusive. (Default: `None`.)

    Returns
    -------
    df_combined : pd.DataFrame
        - Combined DataFrame with columns: `gene`, `duplicate_frequency`, 
        `duplicate_frequency_avg`, `duplicate_frequency_std`, `condition`.
    """
    dfs = []
    for group_id, condition_name, *_ in groups:
        path = f"{repcalc_dir}{group_id}.{processing_stage}.group.{call_type}.tsv"
        df = pd.read_csv(path, sep='\t')
        df = df[(df['level']==level) & (df['mode']==mode) & (df['productive']==productive)]
        df = df.loc[:,['gene', 'duplicate_frequency_avg', 'duplicate_frequency_std', 'N']]
        df['condition'] = condition_name
        dfs.append(df)

    df_combined = pd.concat(dfs, ignore_index=True)
    df_combined = df_combined.sort_values(by=['gene', 'condition'])
    df_combined = df_combined.reset_index(drop=True)
    
    # reduce df by threshold

    ## will keep the group if at least one value in 'gene' is not 0
    # if threshold is not None and isinstance(threshold, float):
    #     df_combined.loc[:, 'duplicate_frequency_avg'] = df_combined['duplicate_frequency_avg'].where(df_combined['duplicate_frequency_avg'] >= threshold, 0)
    #     df_combined = df_combined.groupby('gene').filter(lambda g: not (g['duplicate_frequency_avg'] == 0).all())

    # will drop the group if all values of 'gene' in group are below threshold
    if threshold is not None and isinstance(threshold, float):
        df_combined = df_combined.groupby('gene').filter(
            lambda g: (g['duplicate_frequency_avg'] >= threshold).any()
        )

    return df_combined

def load_gene_usage_data_for_group(
        repcalc_dir:str,
        repertoire_group:dict,
        processing_stage:str,
        call_type:str='v_call',
        level:str='gene',
        mode:str='proportion',
        productive:bool=True
    ) -> pd.DataFrame:
    """
    Load gene usage data for each repertoire in a given group.

    Parameters
    ----------
    repcalc_dir : str
        - RepCalc base directory path.
    repertoire_group : dict
        - Repertoire group object with list of repertoires
    processing_stage : str
        - The middle part of the filename path
    call_type : str
        - The type of call. Can be `v_call`, `d_call`, `j_call`. Default: `v_call`.
    level : str
        - Type of gene level information. Can be `allele`, `gene`, `subgroup`.
        Default: `gene`.
    mode : str
        - Can be `exists`, `proportion`, or `unique`. Default: `proportion`.
    productive: bool
        - Only productive rearrangements (True), both productive and non-productive (False)

    Returns
    -------
    dfs : list of pd.DataFrame, one for each repertoire in group
    """
    dfs = []
    for rep in repertoire_group['repertoires']:
        df = load_gene_usage_data(repcalc_dir, rep['repertoire_id'], processing_stage, call_type, level, mode, productive)
        dfs.append(df)
    return dfs

def load_and_prepare_data_subjects(
        repcalc_dir:str,
        repertoire_groups:dict,
        group_id:str,
        sample_info_df:pd.DataFrame,
        processing_stage:str='.igblast.makedb.allele.clone.',
        level:str='subgroup',
        mode:str='proportion',
        call_type:str='v_call'
    ) -> pd.DataFrame:
    """
    Load and preprocess data from files one group with multiple repertoires.

    Parameters
    ----------
    repcalc_dir : str
        - The base directory path
    repertoire_groups : dict
        - All repertoire groups
    group_id : str
        - The group we want to look at
    sample_info_df : pd.DataFrame
        - Pandas DataFrame containing metadata about samples
    processing_stage : str
        - The middle part of the filename path. 
        Default: `.igblast.makedb.allele.clone.`
    level : str
        - The type of gene level information. Can be `allele`, `gene`, `subgroup`.
        Default: `subgroup`.
    mode : str
        - Can be `exists`, `proportion`, `unique`. Default: `proportion`.
    call_type : str
        - The type of call. Can be `v_call`, `d_call`, `j_call`. Default: `v_call`.

    Returns
    -------
    df_combined : pd.DataFrame
        - Combined DataFrame with columns: `gene`, `duplicate_frequency_avg`,
        `duplicate_frequency_std`, `condition`.
    """
    dfs = []
    for item in repertoire_groups[group_id]['repertoires']:
        repertoire_id = item['repertoire_id']
        tissue = sample_info_df.loc[sample_info_df.repertoire_id==repertoire_id, 'tissue'].values[0]
        path = f"{repcalc_dir}{repertoire_id}{processing_stage}{call_type}.tsv"
        df = pd.read_csv(path, sep='\t')
        df = df[(df['level'] == level) & (df['mode'] == mode) & (df['productive'])]
        df = df.loc[:, ['gene', 'duplicate_frequency']].copy()
        df['tissue'] = tissue
        dfs.append(df)

    df_combined = pd.concat(dfs, ignore_index=True)
    df_combined = df_combined.sort_values(by=['gene', 'tissue']).reset_index(drop=True)
    return df_combined

def load_and_prepare_data_junction_aa(
        repcalc_dir:str,
        groups:list=None,
        processing_stage:str='.igblast.makedb.allele.clone.group.',
        call_type:str='junction_aa_length'
    ) -> pd.DataFrame:
    """
    Load and preprocess data for junction amino acid length.
    
    Parameters
    ----------
    repcalc_dir : str
        - The base directory path.
    groups : list
        - A list of tuples in the format of [(`group_id` , `name`)]
    processing_stage : str
        - The middle part of the filename path. Default: `.igblast.makedb.allele.clone.`
    call_type : str
        - The type of call. (Default: `junction_aa_length`.)

    Returns
    -------
    df_combined : pd.DataFrame
        - Combined DataFrame with columns: `CDR3_LENGTH`, `CDR3_COUNT`, `CDR3_RELATIVE`.
    """
    dfs = []
    if groups is not None:
        for group_id, condition_name, *_ in groups:
            path = f"{repcalc_dir}{group_id}{processing_stage}{call_type}.tsv"
            df = pd.read_csv(path, sep='\t')
            df = df.loc[:, ['CDR3_LENGTH', 'CDR3_COUNT', 'CDR3_RELATIVE']].copy()
            df['condition'] = condition_name
            dfs.append(df)
        df_combined = pd.concat(dfs, ignore_index=True)


    return df_combined

# rename this to load_and_prepare_data_vdj_combo? seems like vj_combo is selected from other options
# need to update `level` parameter in doc_string below to have correct options listed.
def load_and_prepare_data_vj_combo(
        repcalc_dir:str,
        group_id:str=None,
        repertoire_groups:tuple[str, str]=None,
        combo_type:str='vj_combo',
        level:str='subgroup|subgroup',
        mode:str='proportion', 
        productive:str=True,
        processing_stage:str='.igblast.makedb.allele.clone.group.'
    ) -> pd.DataFrame:
    """
    Load and prepare data to be used in the VJ Combo Circos Plot

    Parameters
    ----------
    repcalc_dir : str
        - The base directory path.
    group_id : str
        - The repertoire group id.
    repertoire_groups : list[tuple[str, str]]
        - A list of tupels containg two strings: The repertoire group ID, and the repertoire group name.
    combo_type : str
        - The VDJ combination. (Default `vj_combo`.)
    level : str
        - The type of gene level information. Can be `allele|allele`, `gene|gene`, `subgroup|subgroup`.
        Default: `subgroup|subgroup`.
    mode : str
        - Can be `exists`, `proportion`, `unique`. Default: `proportion`.
    productive : str
        - Filter for productive contigs. Can be `True` or `False`. (Default: `True`.)
    processing_stage : str
        - The middle part of the filename path. Default: `.igblast.makedb.allele.clone.group.`
    
    Returns
    -------
    df : pd.DataFrame
        - Pandas DataFrame to be used for the VJ Combo Circos Plot.
        Contains the coloumns `v_level`, `j_level`, `duplicate_frequency_avg`,
        `duplicate_frequency_std`.
    """
    dfs=[]
    if group_id is not None:
        repertoire_groups = [(group_id, group_id)]
    for group_id, group_name, *_ in repertoire_groups:
        path = f"{repcalc_dir}{group_id}{processing_stage}{combo_type}.tsv"
        df = pd.read_csv(path, sep='\t')
        df['condition'] = [group_name]*df.shape[0]
        dfs.append(df)
    df = pd.concat(dfs)
        
    # check with Tanzira, added ==productive below since productive arg was not in use
    df = df[(df['level']==level) & (df['mode']==mode) & (df['productive']==productive)]
    cols = ['v_level', 'j_level', 'duplicate_frequency_avg', 'duplicate_frequency_std', 'condition']
    return df.loc[:, cols]

def load_diversity_data(
        data_dir:str,
        repertoire_id:str,
        processing_stage:str
    ):
    """
    Load in the diversity data for a repertoire.

    Parameters
    ----------
    data_dir : str
        The data directory path.
    repertoire_id :  str
        The repertoire ID.
    processing_stage : str
        The processing stage.
    
    Returns
    -------
    df : pd.DataFrame
    """
    filename = f"{data_dir}{repertoire_id}.{processing_stage}.diversity.tsv"
    df = pd.read_csv(filename, sep='\t')
    return df

def load_diversity_group_data(
        data_dir:str,
        repertoire_group:dict,
        processing_stage:str
    ):
    """
    Load in the diversity data for a repertoire group.

    Parameters
    ----------
    data_dir : str
        The data directory path.
    repertoire_group :  dict
        The repertoire group extracted from `repertoiers.airr.json`.
    processing_stage : str
        The processing stage.
    
    Returns
    -------
    df : pd.DataFrame
    """
    dfs = []
    for rep in repertoire_group['repertoires']:
        df = load_diversity_data(
            data_dir=data_dir,
            repertoire_id=rep['repertoire_id'],
            processing_stage=processing_stage)
        df['condition'] = repertoire_group['repertoire_group_name']
        dfs.append(df)

    return pd.concat(dfs)

def load_diversity_multiple_group_data(
        data_dir:str,
        repertoire_groups:list[tuple[str, str]],
        processing_stage:str,
        all_groups:dict
    ):
    """
    Load in the diversity data for multiple repertoire groups.

    Parameters
    ----------
    data_dir : str
        The data directory path.
    repertoire_groups :  list[tuple[str, str]]
        Containing the tuple of `repertoire_group_id` and `group_name`. Technically group name is not used here but this format was used for consistency.
    processing_stage : str
        The processing stage.
    all_groups : dict
        The repertoire groups from `repertoiers.airr.json`. Can be constructed using `{ obj['repertoire_group_id'] : obj for obj in data['RepertoireGroup'] }`
    
    Returns
    -------
    df : pd.DataFrame
    """
    dfs = []
    for (group_id, *_) in repertoire_groups:
        df = load_diversity_group_data(
            data_dir=data_dir,
            repertoire_group=all_groups[group_id],
            processing_stage=processing_stage
        )
        dfs.append(df)
    df = pd.concat(dfs)
    return df

# TO DO: Check if group_name is actually a name or an ID. 
# Using repcalc files, I noticed that there is group_id
# def load_and_prepare_mutation_data(
#         mutation_data_dir:str,
#         group_id:str=None,
#         repertoire_groups:tuple[str, str]=None,
#         trim_codons:int=16,
#         processing_stage:str='gene.mutations'
#     ) -> pd.DataFrame:
#     """
#     Load and prepare mutation data for hedgehog plots

#     Parameters
#     ----------
#     mutation_data_dir : str
#         - Path to mutation data directory.
#     group_id : str
#         - The name of the group.
#     trim_codons : int
#         - The number of codons to trim. (Default: `16`.)
#     processing_stage : str
#         -  The current processing stage. Used for the file_path. (Default: `gene.mutations`.)

#     Returns
#     -------
#     data : pd.DataFrame
#         - The data required for the mutational hedgehog plot. Contains columns: `place`,
#         `group`, `value`, `se`, `se_start`, `se_end`.
#     """
#     # --- Load the data ---
#     file_path = f"{mutation_data_dir}{processing_stage}.repertoire_group.frequency.mutational_report.csv"
#     #read the mutation data file
#     group_muts = pd.read_csv(file_path, index_col='repertoire_group_id')

#     # --- Helper functions ---
#     def add_suffix(cols, suffix):
#         return [f"{col}_{suffix}" for col in cols]

#     # --- Define regions ---
#     regions = {
#         "FWR1": range(1, 27),
#         "CDR1": range(27, 39),
#         "FWR2": range(39, 56),
#         "CDR2": range(56, 66),
#         "FWR3": range(66, 105)
#     }

#     # --- Build column lists ---
#     pos_cols_r_aa = [f"mu_freq_{i}_r_aa" for r in regions.values() for i in r]
#     pos_cols_r_aa_avg = add_suffix(pos_cols_r_aa, "avg")
#     pos_cols_r_aa_std = add_suffix(pos_cols_r_aa, "std")
#     pos_cols_r_aa_N = add_suffix(pos_cols_r_aa, "N")

#     # # --- Extract relevant data ---
#     if group_id is not None:
#         avg = group_muts.loc[group_id, pos_cols_r_aa_avg]
#         std = group_muts.loc[group_id, pos_cols_r_aa_std]
#         n = group_muts.loc[group_id, pos_cols_r_aa_N]
#     else:
#         group_ids = []
#         for group_id, group_name, *_ in repertoire_groups:
#             group_ids.append(group_id)
#         avg = group_muts.loc[group_ids, pos_cols_r_aa_avg]
#         std = group_muts.loc[group_ids, pos_cols_r_aa_std]
#         n = group_muts.loc[group_ids, pos_cols_r_aa_N]

#     np.seterr(divide='ignore', invalid='ignore')
#     # can not divide if the index does not match (below)
#     se = pd.Series((std.values / np.sqrt(n.values)), index = std.index)

#     # Replace NaNs
#     se = se.fillna(0)
    
#     # Trim codons
#     avg.iloc[:trim_codons] = 0
#     se.iloc[:trim_codons] = 0
    
#     # --- Build data for plotting ---
#     def build_region_data(region_name, codon_range, avg, se, empty_bar=4):
#         group, place, value, se_vals, se_start, se_end = [], [], [], [], [], []
#         for i in codon_range:
#             col = f"mu_freq_{i}_r_aa"
#             group.append(region_name)
#             place.append(i)
#             v = avg.get(f"{col}_avg", 0)
#             s = se.get(f"{col}_std", 0)
#             value.append(v)
#             se_vals.append(s)
#             se_start.append(max(0, v - s))
#             se_end.append(v + s)
#         return group, place, value, se_vals, se_start, se_end

#     # Compile final data
#     group, place, value, se_vals, se_start, se_end = [], [], [], [], [], []
#     for region_name, codon_range in regions.items():
#         g, p, v, s, ss, se_ = build_region_data(region_name, codon_range, avg, se)
#         group.extend(g)
#         place.extend(p)
#         value.extend(v)
#         se_vals.extend(s)
#         se_start.extend(ss)
#         se_end.extend(se_)

#     data = pd.DataFrame({
#         "place": place,
#         "group": group,
#         "value": value,
#         "se": se_vals,
#         "se_start": se_start,
#         "se_end": se_end
#     })
    
#     # View result
#     return data

def load_and_prepare_mutation_data(
        mutation_data_dir: str,
        group_id: str = None,
        repertoire_groups: list[tuple[str, str]] = None,
        trim_codons: int = 16,
        processing_stage: str = 'gene.mutations'
    ) -> pd.DataFrame:
    """
    Load and prepare mutation data for hedgehog plots

    Parameters
    ----------
    mutation_data_dir : str
        - Path to mutation data directory.
    group_id : str, optional
        - The name of a single group.
    repertoire_groups : list of (group_id, group_name), optional
        - A list of (id, name) tuples for multiple repertoire groups.
    trim_codons : int
        - The number of codons to trim. (Default: `16`.)
    processing_stage : str
        - The current processing stage. Used for the file_path. (Default: `gene.mutations`.)

    Returns
    -------
    data : pd.DataFrame
        - Columns: `place`, `region`, `value`, `se`, `se_start`, `se_end`, `group_name`
    """
    # --- Load the data ---
    file_path = f"{mutation_data_dir}{processing_stage}.repertoire_group.frequency.mutational_report.csv"
    group_muts = pd.read_csv(file_path, index_col='repertoire_group_id')

    # --- Helper functions ---
    def add_suffix(cols, suffix):
        return [f"{col}_{suffix}" for col in cols]

    # --- Define regions ---
    regions = {
        "FWR1": range(1, 27),
        "CDR1": range(27, 39),
        "FWR2": range(39, 56),
        "CDR2": range(56, 66),
        "FWR3": range(66, 105)
    }

    pos_cols_r_aa = [f"mu_freq_{i}_r_aa" for r in regions.values() for i in r]
    pos_cols_r_aa_avg = add_suffix(pos_cols_r_aa, "avg")
    pos_cols_r_aa_std = add_suffix(pos_cols_r_aa, "std")
    pos_cols_r_aa_N = add_suffix(pos_cols_r_aa, "N")

    # --- Function to build per-group data ---
    def build_group_data(group_id: str, group_name: str):
        avg = group_muts.loc[group_id, pos_cols_r_aa_avg]
        std = group_muts.loc[group_id, pos_cols_r_aa_std]
        n = group_muts.loc[group_id, pos_cols_r_aa_N]

        np.seterr(divide='ignore', invalid='ignore')
        se = pd.Series((std.values / np.sqrt(n.values)), index=std.index).fillna(0)

        # Trim codons
        avg.iloc[:trim_codons] = 0
        se.iloc[:trim_codons] = 0

        # Build region data
        rows = []
        for region_name, codon_range in regions.items():
            for i in codon_range:
                col = f"mu_freq_{i}_r_aa"
                v = avg.get(f"{col}_avg", 0)
                s = se.get(f"{col}_std", 0)
                rows.append({
                    "place": i,
                    "region": region_name, # new
                    "value": v,
                    "se": s,
                    "se_start": max(0, v - s),
                    "se_end": v + s,
                    "group_name": group_name
                })
        return pd.DataFrame(rows)

    # --- Single or multiple groups ---
    if group_id is not None:
        return build_group_data(group_id, group_id)
    elif repertoire_groups is not None:
        dfs = [build_group_data(group_id, group_name) for group_id, group_name, *_ in repertoire_groups]
        return pd.concat(dfs, ignore_index=True)
    else:
        raise ValueError("Provide either group_id or repertoire_groups")

def load_shared_cdr3_aa_data(
        repcalc_dir, 
        repertoire_groups, 
        sample_info_df, 
        groups, 
        percentage=False,
        similarity_metric='jaccard',
        sharing='cdr3_aa_sharing',
        processing_stage='repertoire.shared.matrix',
        plot_type='between'
    ):
    """
    Load and process cdr3 data for plotting.
    
    Returns:
    - matrix: processed DataFrame ready for heatmap
    """

    def similarity_calculation(df, similarity_metric):
        '''
        J(A,B)= ∣A∩B∣ / ∣A∪B∣ 
        simpson_coefficient = ∣A∩B∣ / min(∣A∣,∣B∣)
        a_in_b = ∣A∩B∣ / ∣A∣ # % in one group or the other
        b_in_a = ∣A∩B∣ / ∣B∣
        '''
        percentage_matrix = df.copy()
        percentage_matrix = percentage_matrix.astype(float)

        # Loop over the matrix and calculate the percentage of shared overlap
        for i in percentage_matrix.index:
            for j in percentage_matrix.columns:
                # Calculate percentage of shared overlap
                shared_cdr3s = df.at[i, j]
                if similarity_metric == 'jaccard':
                    if shared_cdr3s == 0:
                        percentage_matrix.at[i, j] = 0
                        continue
                    total_cdr3s_i = df.at[i, i] + df.at[j, j] -  shared_cdr3s # total values in i and j - common
                elif similarity_metric == 'simpson_coefficient':
                    total_cdr3s_i = min(df.at[i, i], df.at[j, j]) # Minimum diagonal value for sample i, j
                elif similarity_metric == 'a_in_b':
                    total_cdr3s_i = df.at[i, i] # Minimum diagonal value for sample i
                elif similarity_metric == 'b_in_a':
                    total_cdr3s_i =  df.at[j, j] # Minimum diagonal value for sample j
                    
                pct_overlap = (shared_cdr3s / total_cdr3s_i)*100
                percentage_matrix.at[i, j] = pct_overlap

        np.fill_diagonal(percentage_matrix.values, 0)

        return percentage_matrix

    path = f"{repcalc_dir}{processing_stage}.{sharing}.tsv"
    df = pd.read_csv(path, sep='\t', index_col=0)

    if percentage:
        df = similarity_calculation(df, similarity_metric)
        
    groups_repertoire_ids = {}
    groups_sample_ids = {}
    for group_id, condition, *_ in groups:
        repertoire_ids = []
        sample_ids = []
        for item in repertoire_groups[group_id]['repertoires']:
            repertoire_id = item['repertoire_id'].strip()
            sample_id = sample_info_df.loc[sample_info_df.repertoire_id == repertoire_id, 'ss_id'].values[0]
            repertoire_ids.append(repertoire_id)
            sample_ids.append(sample_id)
        groups_repertoire_ids[condition] = repertoire_ids
        groups_sample_ids[condition] = sample_ids

    if plot_type == 'within':
        # One group only
        assert len(groups_repertoire_ids) == 1, "Within group plot requires exactly 1 group"
        cond1 = cond2 = list(groups_repertoire_ids.keys())[0]
    elif plot_type == 'between':
        # Two groups
        assert len(groups_repertoire_ids) == 2, "Between group plot requires exactly 2 groups"
        cond1, cond2 = list(groups_repertoire_ids.keys())
    else:
        raise ValueError("plot_type must be 'within' or 'between'")
        
    ids1 = groups_repertoire_ids[cond1]
    ids2 = groups_repertoire_ids[cond2]
    matrix = df.loc[ids1, ids2]
    matrix.index = groups_sample_ids[cond1]
    matrix.columns = groups_sample_ids[cond2]
    matrix.index.name = cond1
    matrix.columns.name = cond2
    return matrix

def load_cdr3_sequences_each_repertoire(
        repcalc_dir,
        repertoire_id,
        sharing='cdr3_aa_sharing',
        processing_stage='repertoire'
    ):
    """
    Load and process CDR3 sequence data for a single repertoire.

    Parameters:
    - repcalc_dir (str): Directory where files are stored
    - repertoire_id (str): Repertoire UUID
    - sharing (str): Type of sharing metric (default 'cdr3_aa_sharing')
    - processing_stage (str): Processing stage in filename (default 'repertoire')

    Returns:
    - df (DataFrame): Processed top-N CDR3 sequences
    """
    path = f'{repcalc_dir}{repertoire_id}.{processing_stage}.{sharing}.tsv'
    df = pd.read_csv(path, sep='\t')
    return df

def load_unique_cdr3_sequence_numbers(repcalc_dir, repertoire_groups, groups, sample_info_df, sharing='cdr3_aa_sharing', processing_stage='repertoire'):
    groups_unique_cdr3_counts = []
    for group_id, condition, *_ in groups:
        for item in repertoire_groups[group_id]['repertoires']:
            repertoire_id = item['repertoire_id'].strip()
            df = load_cdr3_sequences_each_repertoire( repcalc_dir=repcalc_dir, repertoire_id=repertoire_id)
            sample_id = sample_info_df.loc[sample_info_df.repertoire_id == repertoire_id, 'sample_id'].values[0]
            n_unique_cdr3 = df['junction_aa'].nunique()
            groups_unique_cdr3_counts.append({
                'condition': condition,
                'repertoire_id': repertoire_id,
                'sample_id': sample_id,
                'n_unique_cdr3': n_unique_cdr3
            })
    # Convert list of dicts to DataFrame
    df = pd.DataFrame(groups_unique_cdr3_counts)
    
    return df
    
def load_and_prepare_group_summary_comparison(
        repcalc_dir,
        search_ids,
        processing_stage='group.summary_comparison',
        sharing='cdr3_aa_sharing'
    ):

    def parse_grouped_matrices(
            file_path
        ):

        groups = []
        with open(file_path, 'r') as f:
            lines = f.readlines()
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith("GROUPS"):
                # Extract group metadata
                group_info = line.split("\t")[1:]
                # Next line is the header
                i += 1
                header = lines[i].strip().split("\t")[1:]  # skip "SHARED"
                # Now read the matrix rows
                matrix_data = []
                i += 1
                while i < len(lines) and not lines[i].startswith("GROUPS"):
                    row = lines[i].strip().split("\t")
                    if len(row) < 2:
                        i += 1
                        continue
                    matrix_data.append([int(row[0])] + [int(x) for x in row[1:]])
                    i += 1
                # Build DataFrame
                df = pd.DataFrame(matrix_data, columns=["ID"] + header)
                df.set_index('ID', inplace = True)
                groups.append({
                    "group_ids": group_info,
                    "data": df
                })
            else:
                i += 1

        return groups

    path = f'{repcalc_dir}{processing_stage}.{sharing}.tsv'

    ## Parse the data
    parsed_groups = parse_grouped_matrices(path)
    df = None
    for group in parsed_groups:
        if search_ids == group["group_ids"]:
            df = group['data']
            return df
    if df == None:
        raise ValueError("Groups not found.")
    
    return df

def load_and_prepare_clonal_abundance_data(data_dir, repertoires, clone_tool='repcalc', processing_stage='igblast.allele.clone', internal=False):
    dfs = []
    for repertoire_id, condition_name, sample_id in repertoires:
        if clone_tool == 'repcalc':
            path = f"{data_dir}{repertoire_id}.{processing_stage}.count.tsv"
            col = 'copy_freq'
        if clone_tool == 'changeo':
            path = f"{data_dir}{repertoire_id}.{processing_stage}.abundance.tsv"
            col = 'p'
        df = pd.read_csv(path, sep='\t')
        df['rank'] = np.arange(1, len(df)+1)
        ## Rename the column name to have consistency in the dataframe
        df = df.rename(columns={col: 'abundance'})
        df['cumulative_abundance'] = df.abundance.transform('cumsum')
        df['condition'] = condition_name
        # df['clones'] = df.shape[0]
        if internal:
            df['sample_id'] = sample_id
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)    
    return df

def load_and_prepare_clonal_abundance_group_data(data_dir, groups, all_groups, all_repertoires, clone_tool='repcalc', processing_stage='igblast.allele.clone'):
    dfs = []
    for (group_id, *_) in groups:
        repertoires = []
        for rep in all_groups[group_id]['repertoires']:
            repertoires.append([rep['repertoire_id'], all_groups[group_id]['repertoire_group_name'], all_repertoires[rep['repertoire_id']]['sample'][0]['sample_id']])
        
        df = load_and_prepare_clonal_abundance_data(
            data_dir=data_dir,
            repertoires=repertoires,
            clone_tool=clone_tool,
            processing_stage=processing_stage,
            internal=True
        )
    
        dfs.append(df)
    df = pd.concat(dfs)
    return df

# def load_shared_cdr3_aa_data(
#         repcalc_dir, 
#         repertoire_groups, 
#         sample_info_df, 
#         groups, 
#         percentage=False,
#         similarity_metric='jaccard',
#         sharing='cdr3_aa_sharing',
#         processing_stage='repertoire.shared.matrix',
#         plot_type='between'
#     ):
#     """
#     Load and process cdr3 data for plotting.
    
#     Returns:
#     - matrix: processed DataFrame ready for heatmap
#     """

#     path = f"{repcalc_dir}{processing_stage}.{sharing}.tsv"
#     df = pd.read_csv(path, sep='\t', index_col=0)

#     if percentage:
#         df = similarity_calculation(df, similarity_metric)
        
#     groups_repertoire_ids = {}
#     groups_sample_ids = {}
#     for group_id, condition, *_ in groups:
#         repertoire_ids = []
#         sample_ids = []
#         for item in repertoire_groups[group_id]['repertoires']:
#             repertoire_id = item['repertoire_id'].strip()
#             sample_id = sample_info_df.loc[sample_info_df.repertoire_id == repertoire_id, 'ss_id'].values[0]
#             repertoire_ids.append(repertoire_id)
#             sample_ids.append(sample_id)
#         groups_repertoire_ids[condition] = repertoire_ids
#         groups_sample_ids[condition] = sample_ids

#     if plot_type == 'within':
#         # One group only
#         assert len(groups_repertoire_ids) == 1, "Within group plot requires exactly 1 group"
#         cond1 = cond2 = list(groups_repertoire_ids.keys())[0]
#     elif plot_type == 'between':
#         # Two groups
#         assert len(groups_repertoire_ids) == 2, "Between group plot requires exactly 2 groups"
#         cond1, cond2 = list(groups_repertoire_ids.keys())
#     else:
#         raise ValueError("plot_type must be 'within' or 'between'")
        
#     ids1 = groups_repertoire_ids[cond1]
#     ids2 = groups_repertoire_ids[cond2]
#     matrix = df.loc[ids1, ids2]
#     matrix.index = groups_sample_ids[cond1]
#     matrix.columns = groups_sample_ids[cond2]
#     matrix.index.name = cond1
#     matrix.columns.name = cond2
    
#     return matrix
