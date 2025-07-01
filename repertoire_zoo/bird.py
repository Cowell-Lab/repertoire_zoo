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

def read_sample_info_airr(
        path:str='repertoires.airr.json'
    ) -> pd.DataFrame:
    """
    Creates a Pandas DataFrame containing sample information from the repertoire file.

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
                'diagnosis':diagnosis
            }
        )
    df = pd.DataFrame(info_list)

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
            'description':obj['description'],
            'repertoires':obj['repertoires'],
            'repertoire_count':len(obj['repertoires'])
        }
        for obj in repertoire_group_data
    ]
    df = pd.DataFrame(repertoire_group_list)

    if not silence:
        print("Total Number of Repertoire Groups: ", len(df))

    return df

def load_and_prepare_data_groups(
        repcalc_dir:str,
        groups:list,
        call_type:str='v_call',
        level:str='gene',
        mode:str='proportion',
        processing_stage:str='.igblast.makedb.allele.clone.group.'
    ) -> pd.DataFrame:
    """
    Load and preprocess data from files for multiple groups.

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

    Returns
    -------
    df_combined : pd.DataFrame
        - Combined DataFrame with columns: `gene`, `duplicate_frequency_avg`, 
        `duplicate_frequency_std`, `condition`.
    """
    dfs = []
    for group_id, condition_name in groups:
        path = f"{repcalc_dir}{group_id}{processing_stage}{call_type}.tsv"
        df = pd.read_csv(path, sep='\t')
        df = df[(df['level']==level) & (df['mode']==mode) & df['productive']]
        df.loc[:,['gene', 'duplicate_frequency_avg', 'duplicate_frequency_std']]
        df['condition'] = condition_name
        dfs.append(df)

    df_combined = pd.concat(dfs, ignore_index=True)
    df_combined = df_combined.sort_values(by=['gene', 'condition'])
    df_combined = df_combined.reset_index(drop=True)

    return df_combined

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
        groups:list,
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
    for group_id, condition_name in groups:
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
        group_id:str,
        combo_type:str='vj_combo',
        level:str='subgroup|subgroup',
        mode:str='proportion', 
        productive:str='TRUE',
        processing_stage:str='.igblast.makedb.allele.clone.group.'
    ) -> pd.DataFrame:
    """
    Load and prepare data to be used in the VJ Combo Circos Plot

    Parameters
    ----------
    repcalc_dir : str
        - The base directory path.
    group_id : str
        - The group ID.
    combo_type : str
        - The VDJ combination. (Default `vj_combo`.)
    level : str
        - The type of gene level information. Can be `allele`, `gene`, `subgroup`.
        Default: `subgroup|subgroup`.
    mode : str
        - Can be `exists`, `proportion`, `unique`. Default: `proportion`.
    productive : str
        - Filter for productive contigs. Can be `TRUE` or `FALSE`. (Default: `TRUE`.)
    processing_stage : str
        - The middle part of the filename path. Default: `.igblast.makedb.allele.clone.group.`
    
    Returns
    -------
    df : pd.DataFrame
        - Pandas DataFrame to be used for the VJ Combo Circos Plot.
        Contains the coloumns `v_level`, `j_level`, `duplicate_frequency_avg`,
        `duplicate_frequency_std`.
    """
    path = f"{repcalc_dir}{group_id}{processing_stage}{combo_type}.tsv"
    df = pd.read_csv(path, sep='\t')
    # check with Tanzira, added ==productive below since productive arg was not in use
    df = df[(df['level']==level) & (df['mode']==mode) & (df['productive']==productive)]
    cols = ['v_level', 'j_level', 'duplicate_frequency_avg', 'duplicate_frequency_std']
    return df.loc[:, cols]

def load_and_prepare_mutation_data(
        mutation_data_dir:str,
        group_name:str,
        trim_codons:int=16,
        processing_stage:str='gene.mutations'
    ) -> pd.DataFrame:
    """
    Load and prepare mutation data for hedgehog plots

    Parameters
    ----------
    mutation_data_dir : str
        - Path to mutation data directory.
    group_name : str
        - The name of the group.
    trim_codons : int
        - The number of codons to trim. (Default: `16`.)
    processing_stage : str
        -  The current processing stage. Used for the file_path. (Default: `gene.mutations`.)

    Returns
    -------
    data : pd.DataFrame
        - The data required for the mutational hedgehog plot. Contains columns: `place`,
        `group`, `value`, `se`, `se_start`, `se_end`.
    """
    # --- Load the data ---
    file_path = f"{mutation_data_dir}{processing_stage}+\
        .repertoire_group.frequency.mutational_report.csv"
    #read the mutation data file
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

    # --- Build column lists ---
    pos_cols_r_aa = [f"mu_freq_{i}_r_aa" for r in regions.values() for i in r]
    pos_cols_r_aa_avg = add_suffix(pos_cols_r_aa, "avg")
    pos_cols_r_aa_std = add_suffix(pos_cols_r_aa, "std")
    pos_cols_r_aa_N = add_suffix(pos_cols_r_aa, "N")

    # # --- Extract relevant data ---
    avg = group_muts.loc[group_name, pos_cols_r_aa_avg]
    std = group_muts.loc[group_name, pos_cols_r_aa_std]
    n = group_muts.loc[group_name, pos_cols_r_aa_N]
    np.seterr(divide='ignore', invalid='ignore')
    # can not divide if the index does not match (below)
    se = pd.Series((std.values / np.sqrt(n.values)), index = std.index)

    # Replace NaNs
    se = se.fillna(0)
    # Trim codons
    avg.iloc[:trim_codons] = 0
    se.iloc[:trim_codons] = 0
    # --- Build data for plotting ---
    def build_region_data(region_name, codon_range, avg, se, empty_bar=4):
        group, place, value, se_vals, se_start, se_end = [], [], [], [], [], []
        for i in codon_range:
            col = f"mu_freq_{i}_r_aa"
            group.append(region_name)
            place.append(i)
            v = avg.get(f"{col}_avg", 0)
            s = se.get(f"{col}_std", 0)
            value.append(v)
            se_vals.append(s)
            se_start.append(max(0, v - s))
            se_end.append(v + s)
        return group, place, value, se_vals, se_start, se_end

    # Compile final data
    group, place, value, se_vals, se_start, se_end = [], [], [], [], [], []
    for region_name, codon_range in regions.items():
        g, p, v, s, ss, se_ = build_region_data(region_name, codon_range, avg, se)
        group.extend(g)
        place.extend(p)
        value.extend(v)
        se_vals.extend(s)
        se_start.extend(ss)
        se_end.extend(se_)

    data = pd.DataFrame({
        "place": place,
        "group": group,
        "value": value,
        "se": se_vals,
        "se_start": se_start,
        "se_end": se_end
    })
    # View result
    return data
