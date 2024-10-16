
import os
import pathlib
if "__file__" in globals():#only run in the translated python file, as __file__ is not defined with ipython
    INTABLE_CONFIG = os.path.join(pathlib.Path(__file__).parent.absolute(), "configs", "intable_config.yaml") #the yaml config is located one directory below the python library files
    CONFIG_PATH = os.path.join(pathlib.Path(__file__).parent.absolute(), "configs")

import logging
import directlfq.config as config
import directlfq.utils_fileread as utils_fileread

LOGGER = logging.getLogger(__name__)


def get_samples_used_from_samplemap_file(samplemap_file, cond1, cond2):
    samplemap_df = load_samplemap(samplemap_file)
    return get_samples_used_from_samplemap_df(samplemap_df, cond1, cond2)


def get_samples_used_from_samplemap_df(samplemap_df, cond1, cond2):
    samples_c1 = samplemap_df[[cond1 == x for x in samplemap_df["condition"]]]["sample"] #subset the df to the condition
    samples_c2 = samplemap_df[[cond2 == x for x in samplemap_df["condition"]]]["sample"]
    return list(samples_c1), list(samples_c2)

def get_all_samples_from_samplemap_df(samplemap_df):
    return list(samplemap_df["sample"])

# %% ../nbdev_nbs/04_utils.ipynb 6
import pandas as pd

def get_samplenames_from_input_df(data):
    """extracts the names of the samples of the AQ input dataframe"""
    names = list(data.columns)
    names.remove(config.PROTEIN_ID)
    names.remove(config.QUANT_ID)
    return names

# %% ../nbdev_nbs/04_utils.ipynb 7
import numpy as np
def filter_df_to_minrep(quant_df_wideformat, samples_c1, samples_c2, minrep):
    """filters dataframe in directlfq format such that each column has a minimum number of replicates
    """
    quant_df_wideformat = quant_df_wideformat.replace(0, np.nan)
    df_c1_minrep = quant_df_wideformat[samples_c1].dropna(thresh = minrep, axis = 0)
    df_c2_minrep = quant_df_wideformat[samples_c2].dropna(thresh = minrep, axis = 0)
    idxs_both = df_c1_minrep.index.intersection(df_c2_minrep.index)
    quant_df_reduced = quant_df_wideformat.iloc[idxs_both].reset_index()
    return quant_df_reduced

# %% ../nbdev_nbs/04_utils.ipynb 8
def get_condpairname(condpair):
    return f"{condpair[0]}_VS_{condpair[1]}"

# %% ../nbdev_nbs/04_utils.ipynb 9
def get_quality_score_column(acquisition_info_df):
    if "FG.ShapeQualityScore" in acquisition_info_df.columns:
        param = "FG.ShapeQualityScore"
    elif "Quantity.Quality" in acquisition_info_df.columns:
        param = "Quantity.Quality"
    return param

# %% ../nbdev_nbs/04_utils.ipynb 10
import os

def make_dir_w_existcheck(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

# %% ../nbdev_nbs/04_utils.ipynb 11
import os
def get_results_plot_dir_condpair(results_dir, condpair):
    results_dir_plots = f"{results_dir}/{condpair}_plots"
    make_dir_w_existcheck(results_dir_plots)
    return results_dir_plots

# %% ../nbdev_nbs/04_utils.ipynb 12
def get_middle_elem(sorted_list):
    nvals = len(sorted_list)
    if nvals==1:
        return sorted_list[0]
    middle_idx = nvals//2
    if nvals%2==1:
        return sorted_list[middle_idx]
    return 0.5* (sorted_list[middle_idx] + sorted_list[middle_idx-1])

# %% ../nbdev_nbs/04_utils.ipynb 13
import numpy as np
def get_nonna_array(array_w_nas):
    res = []
    isnan_arr = np.isnan(array_w_nas)

    for idx in range(len(array_w_nas)):
        sub_res = []
        sub_array = array_w_nas[idx]
        na_array = isnan_arr[idx]
        for idx2 in range(len(sub_array)):
            if not na_array[idx2]:
               sub_res.append(sub_array[idx2])
        res.append(np.array(sub_res))
    return np.array(res)


# %% ../nbdev_nbs/04_utils.ipynb 15
# %% ../nbdev_nbs/04_utils.ipynb 16
def invert_dictionary(my_map):
    inv_map = {}
    for k, v in my_map.items():
        inv_map[v] = inv_map.get(v, []) + [k]
    return inv_map


# %% ../nbdev_nbs/04_utils.ipynb 18
def count_fraction_outliers_from_expected_fc(result_df, threshold, expected_log2fc):
    num_outliers = sum([abs(x-expected_log2fc)> threshold for x in result_df["log2fc"]])
    fraction_outliers = num_outliers/len(result_df["log2fc"])
    LOGGER.info(f"{round(fraction_outliers, 2)} outliers")
    return fraction_outliers

# %% ../nbdev_nbs/04_utils.ipynb 19
import os
import shutil
def create_or_replace_folder(folder):
    if os.path.exists(folder):
        shutil.rmtree(folder)
    os.makedirs(folder)

# %% ../nbdev_nbs/04_utils.ipynb 20
def add_mq_protein_group_ids_if_applicable_and_obtain_annotated_file(mq_file, input_type_to_use ,mq_protein_group_file, columns_to_add):
    try:
        input_type = _get_input_type(mq_file, input_type_to_use)
        if ("maxquant_evidence" in input_type or "maxquant_peptides" in input_type) and ("aq_reformat" not in mq_file) and ("directlfq" not in input_type_to_use) and (input_type_to_use != "directlfq"):
            if mq_protein_group_file is None:
                LOGGER.info("You provided a MaxQuant peptide or evidence file as input. To have the identical ProteinGroups as in the MaxQuant analysis, please provide the ProteinGroups.txt file as well.")
                return mq_file
            else:
                mq_df = load_input_file_and_de_duplicate_if_evidence(mq_file, input_type, columns_to_add)
                id_column = determine_id_column_from_input_df(mq_df)
                id2protein_df = create_id_to_protein_df(mq_protein_group_file, id_column)
                annotated_mq_df = annotate_mq_df(mq_df, id2protein_df, id_column)
                annotated_mq_filename = f"{mq_file}.protgroup_annotated.tsv"
                save_annotated_mq_df(annotated_mq_df, annotated_mq_filename)
                return annotated_mq_filename
        else:
            return mq_file
    except:
        return mq_file


def _get_input_type(mq_file ,input_type_to_use):
    if input_type_to_use is not None:
        return input_type_to_use
    else:
        return get_input_type_and_config_dict(mq_file)[0]
    

def load_input_file_and_de_duplicate_if_evidence(input_file, input_type, columns_to_add):
    input_df = pd.read_csv(input_file, sep = "\t")
    if "maxquant_evidence" in input_type:
        subset_columns = ['id','Sequence','Modified sequence', 'Experiment','Charge', 'Raw file', 'Gene names', 'Intensity', 'Reverse', 'Potential contaminant'] + columns_to_add
        columns_to_group_by = ['Sequence','Modified sequence', 'Experiment','Charge', 'Raw file']
        input_df = input_df[subset_columns].set_index(columns_to_group_by)
        input_df_grouped = input_df.groupby(columns_to_group_by).Intensity.max()
        input_df_no_intensities = input_df.drop(columns=["Intensity"])

        input_df = input_df_no_intensities.merge(input_df_grouped, how= 'right', left_index=True, right_index=True).reset_index()
        input_df = input_df.drop_duplicates(subset=columns_to_group_by)

    return input_df

def create_id_to_protein_df(mq_protein_group_file, id_column):    
    id_mapping_df = pd.read_csv(mq_protein_group_file, sep = "\t", usecols=["Protein IDs", id_column])
    #apply lambda function to id column to split it into a list of ids
    id_mapping_df[id_column] = id_mapping_df[id_column].apply(lambda x: x.split(";"))
    #explode the id column
    id_mapping_df = id_mapping_df.explode(id_column) #https://stackoverflow.com/questions/12680754/split-explode-pandas-dataframe-string-entry-to-separate-rows
    return id_mapping_df


def determine_id_column_from_input_df(input_df):
    input_file_columns = input_df.columns
    num_cols_starting_w_intensity = sum([x.startswith("Intensity ") for x in input_file_columns])
    if num_cols_starting_w_intensity>0:
        return "Peptide IDs"
    else:
        return "Evidence IDs"


def annotate_mq_df(mq_df, id2protein_df, id_column):
    #set dtype of id to string
    mq_df["id"] = mq_df["id"].astype(str)
    id2protein_df = remove_ids_not_occurring_in_mq_df(id2protein_df, mq_df, id_column)
    return mq_df.merge(id2protein_df, how = "right",  left_on = "id", right_on = id_column, suffixes=('', '_y'))

def remove_ids_not_occurring_in_mq_df(id2protein_df, mq_df, id_column):
    mq_df_ids = set(mq_df["id"])
    id2protein_df = id2protein_df[id2protein_df[id_column].isin(mq_df_ids)]
    return id2protein_df

def save_annotated_mq_df(annotated_mq_df, annotated_mq_file):
    annotated_mq_df.to_csv(annotated_mq_file, sep = "\t", index = False)



# %% ../nbdev_nbs/04_utils.ipynb 21


def add_columns_to_lfq_results_table(lfq_results_df, input_file, columns_to_add):
    input_type, config_dict, _ = get_input_type_and_config_dict(input_file)

    input_file = clean_input_filename_if_necessary(input_file)

    protein_column_input_table = get_protein_column_input_table(config_dict)
    standard_columns_for_input_type = get_standard_columns_for_input_type(input_type)

    all_columns = columns_to_add + [protein_column_input_table] + standard_columns_for_input_type
    all_columns = filter_columns_to_existing_columns(all_columns, input_file)

    lfq_results_df = lfq_results_df[[x is not None for x in lfq_results_df[config.PROTEIN_ID]]]
    if (len(columns_to_add) == 0) and (len(standard_columns_for_input_type)==0) :
        return lfq_results_df 
    input_df = pd.read_csv(input_file, sep="\t", usecols=all_columns).drop_duplicates(subset=protein_column_input_table)

    length_before = len(lfq_results_df.index)
    lfq_results_df_appended = pd.merge(lfq_results_df, input_df, left_on=config.PROTEIN_ID, right_on=protein_column_input_table, how='left')
    length_after = len(lfq_results_df_appended.index)

    #lfq_results_df_appended = lfq_results_df_appended.set_index(config.PROTEIN_ID)
    

    assert length_before == length_after
    return lfq_results_df_appended

def clean_input_filename_if_necessary(input_file):
    if "aq_reformat.tsv" in input_file:
        input_file = get_original_file_from_aq_reformat(input_file)
    return input_file

def get_protein_column_input_table(config_dict):
    return config_dict["protein_cols"][0]

def get_standard_columns_for_input_type(input_type):
    
    if 'maxquant' in input_type:
        return ["Gene names"]
    elif 'diann' in input_type:
        return ["Protein.Names", "Genes"]
    elif 'spectronaut' in input_type:
        return ['PG.Genes']
    else:
        return []

def filter_columns_to_existing_columns(columns, input_file):
    existing_columns =  utils_fileread.read_columns_from_file(input_file)
    return [x for x in columns if x in existing_columns]



#function that shows the differing rows between two dataframes
def show_diff(df1, df2):
    return df1.merge(df2, indicator=True, how='outer').loc[lambda x : x['_merge']!='both']


def write_chunk_to_file(chunk, filepath ,write_header):
    """write chunk of pandas dataframe to a file"""
    chunk.to_csv(filepath, header=write_header, mode='a', sep = "\t", index = None)

def index_and_log_transform_input_df(data_df):
    data_df = data_df.set_index([config.PROTEIN_ID, config.QUANT_ID])
    return np.log2(data_df.replace(0, np.nan))

def remove_allnan_rows_input_df(data_df):
    return data_df.dropna(axis = 0, how = 'all')

def remove_potential_quant_id_duplicates(data_df : pd.DataFrame):
    """
    Remove duplicate entries from a DataFrame based on the QUANT_ID column.

    This function removes duplicate rows from the input DataFrame, keeping only the first
    occurrence of each unique QUANT_ID. It also logs a warning message if any duplicates
    are found and removed.

    Args:
        data_df (pd.DataFrame): dataframe in directLFQ format

    Returns:
        pd.DataFrame: dataframe in directLFQ format w duplicate QUANT_ID entries removed.
    """
    before_drop = len(data_df)
    data_df = data_df.drop_duplicates(subset=config.QUANT_ID, keep='first')
    after_drop = len(data_df)
    if before_drop != after_drop:
        entries_removed = before_drop - after_drop
        LOGGER.warning(f"Duplicate quant_ids detected. {entries_removed} rows removed from input df.")

    return data_df


def sort_input_df_by_protein_id(data_df):
    return data_df.sort_values(by = config.PROTEIN_ID,ignore_index=True)


    

# %% ../nbdev_nbs/04_utils.ipynb 29
import yaml
import itertools

def get_relevant_columns(protein_cols, ion_cols, sample_ID, quant_ID, filter_dict):
    filtcols = []
    for filtconf in filter_dict.values():
        filtcols.append(filtconf.get('param'))
    relevant_cols = protein_cols + ion_cols + [sample_ID] + [quant_ID] + filtcols
    relevant_cols = list(set(relevant_cols)) # to remove possible redudancies
    return relevant_cols


def get_relevant_columns_config_dict(config_typedict):
    filtcols = []
    dict_ioncols = []
    for filtconf in config_typedict.get('filters', {}).values():
        filtcols.append(filtconf.get('param'))

    if 'ion_hierarchy' in config_typedict.keys():
        for headr in config_typedict.get('ion_hierarchy').values():
            ioncols = list(itertools.chain.from_iterable(headr.get("mapping").values()))
            dict_ioncols.extend(ioncols)

    quant_ids = get_quant_ids_from_config_dict(config_typedict)
    sample_ids = get_sample_ids_from_config_dict(config_typedict)
    channel_ids = get_channel_ids_from_config_dict(config_typedict)
    relevant_cols = config_typedict.get("protein_cols") + config_typedict.get("ion_cols", []) + sample_ids + quant_ids + filtcols + dict_ioncols + channel_ids
    relevant_cols = list(set(relevant_cols)) # to remove possible redudancies
    return relevant_cols

def get_quant_ids_from_config_dict(config_typedict):
    quantID = config_typedict.get("quant_ID")
    if type(quantID) ==type("string"):
        return [config_typedict.get("quant_ID")]
    if quantID == None:
        return[]
    else:
        return list(config_typedict.get("quant_ID").values())

def get_sample_ids_from_config_dict(config_typedict):
    sampleID = config_typedict.get("sample_ID")
    if type(sampleID) ==type("string"):
        return [config_typedict.get("sample_ID")]
    if sampleID == None:
        return []
    else:
        return config_typedict.get("sample_ID")

def get_channel_ids_from_config_dict(config_typedict):
    return config_typedict.get("channel_ID", [])



def load_config(config_yaml):
    with open(config_yaml, 'r') as stream:
        config_all = yaml.safe_load(stream)
    return config_all

def get_type2relevant_cols(config_all):
    type2relcols = {}
    for type in config_all.keys():
        config_typedict = config_all.get(type)
        relevant_cols = get_relevant_columns_config_dict(config_typedict)
        type2relcols[type] = relevant_cols
    return type2relcols

# %% ../nbdev_nbs/04_utils.ipynb 31
def filter_input(filter_dict, input):
    if filter_dict == None:
        return input
    for filtname,filterconf in filter_dict.items():
        param = filterconf.get('param')
        comparator = filterconf.get('comparator')
        value = filterconf.get('value')

        if comparator not in [">",">=", "<", "<=", "==", "!="]:
            raise TypeError(f"cannot identify the filter comparator of {filtname} given in the longtable config yaml!")

        if comparator=="==":
            input = input[input[param] ==value]
            continue
        try:
            input = input.astype({f"{param}" : "float"})
        except:
            pass

        if comparator==">":
            input = input[input[param].astype(type(value)) >value]

        if comparator==">=":
            input = input[input[param].astype(type(value)) >=value]

        if comparator=="<":
            input = input[input[param].astype(type(value)) <value]

        if comparator=="<=":
            input = input[input[param].astype(type(value)) <=value]

        if comparator=="!=":
            input = input[input[param].astype(type(value)) !=value]

    return input

# %% ../nbdev_nbs/04_utils.ipynb 32
def merge_protein_and_ion_cols(input_df, config_dict):
    protein_cols =  config_dict.get("protein_cols")
    ion_cols = config_dict.get("ion_cols")
    input_df[config.PROTEIN_ID] = join_columns(input_df, protein_cols)
    input_df[config.QUANT_ID] = join_columns(input_df, ion_cols)

    input_df = input_df.rename(columns = {config_dict.get('quant_ID') : "quant_val"})
    return input_df

# %% ../nbdev_nbs/04_utils.ipynb 33
import copy
def merge_protein_cols_and_ion_dict(input_df, config_dict):
    """[summary]
    
    Args:
        input_df ([pandas dataframe]): longtable containing peptide intensity data
        confid_dict ([dict[String[]]]): nested dict containing the parse information. derived from yaml file

    Returns:
        pandas dataframe: longtable with newly assigned config.PROTEIN_ID and config.QUANT_ID columns
    """
    protein_cols = config_dict.get("protein_cols")
    ion_hierarchy = config_dict.get("ion_hierarchy")
    splitcol2sep = config_dict.get('split_cols')
    quant_id_dict = config_dict.get('quant_ID')

    ion_dfs = []
    #concatenate multiple protein columns into one
    input_df[config.PROTEIN_ID] = join_columns(input_df, protein_cols)

    input_df = input_df.drop(columns = [x for x in protein_cols if x!=config.PROTEIN_ID])
    for hierarchy_type in ion_hierarchy.keys():
        df_subset = input_df.copy()
        ion_hierarchy_local = ion_hierarchy.get(hierarchy_type).get("order")
        ion_headers_merged, ion_headers_grouped = get_ionname_columns(ion_hierarchy.get(hierarchy_type).get("mapping"), ion_hierarchy_local) #ion headers merged is just a helper to select all relevant rows, ionheaders grouped contains the sets of ionstrings to be merged into a list eg [[SEQ, MOD], [CH]]
        quant_columns = get_quantitative_columns(df_subset, hierarchy_type, config_dict, ion_headers_merged)
        headers = list(set(ion_headers_merged + quant_columns + [config.PROTEIN_ID]))
        if "sample_ID" in config_dict.keys():
            headers+=[config_dict.get("sample_ID")]
        df_subset = df_subset[headers].drop_duplicates()

        if splitcol2sep is not None:
            if quant_columns[0] in splitcol2sep.keys(): #in the case that quantitative values are stored grouped in one column (e.g. msiso1,msiso2,msiso3, etc.), reformat accordingly
                df_subset = split_extend_df(df_subset, splitcol2sep)
            ion_headers_grouped = adapt_headers_on_extended_df(ion_headers_grouped, splitcol2sep)

        #df_subset = df_subset.set_index(quant_columns)

        df_subset = add_merged_ionnames(df_subset, ion_hierarchy_local, ion_headers_grouped, quant_id_dict, hierarchy_type)
        ion_dfs.append(df_subset)
    input_df = pd.concat(ion_dfs, ignore_index=True)
    return input_df

def join_columns(df, columns, separator='_'):
    if len(columns) == 1:
        return df[columns[0]].fillna('nan').astype(str)
    else:
        return df[columns].fillna('nan').astype(str).agg(separator.join, axis=1)


def get_quantitative_columns(input_df, hierarchy_type, config_dict, ion_headers_merged):
    naming_columns = ion_headers_merged + [config.PROTEIN_ID]
    if config_dict.get("format") == 'longtable':
        quantcol = config_dict.get("quant_ID").get(hierarchy_type)
        return [quantcol]

    if config_dict.get("format") == 'widetable':
        quantcolumn_candidates = [x for x in input_df.columns if x not in naming_columns]
        if "quant_pre_or_suffix" in config_dict.keys():
            return [x for x in quantcolumn_candidates if x.startswith(config_dict.get("quant_pre_or_suffix")) or x.endswith(config_dict.get("quant_pre_or_suffix"))] # in the case that the quantitative columns have a prefix (like "Intensity " in MQ peptides.txt), only columns with the prefix are filtered
        else:
            return quantcolumn_candidates #in this case, we assume that all non-ionname/proteinname columns are quantitative columns


def get_ionname_columns(ion_dict, ion_hierarchy_local):
    ion_headers_merged = []
    ion_headers_grouped = []
    for lvl in ion_hierarchy_local:
        vals = ion_dict.get(lvl)
        ion_headers_merged.extend(vals)
        ion_headers_grouped.append(vals)
    return ion_headers_merged, ion_headers_grouped


def adapt_headers_on_extended_df(ion_headers_grouped, splitcol2sep):
    #in the case that one column has been split, we need to designate the "naming" column
    ion_headers_grouped_copy = copy.deepcopy(ion_headers_grouped)
    for vals in ion_headers_grouped_copy:
        if splitcol2sep is not None:
            for idx in range(len(vals)):
                if vals[idx] in splitcol2sep.keys():
                    vals[idx] = vals[idx] + "_idxs"
    return ion_headers_grouped_copy

def split_extend_df(input_df, splitcol2sep, value_threshold=10):
    """reformats data that is stored in a condensed way in a single column. For example isotope1_intensity;isotope2_intensity etc. in Spectronaut

    Args:
        input_df ([type]): [description]
        splitcol2sep ([type]): [description]
        value_threshold([type]): [description]

    Returns:
        Pandas Dataframe: Pandas dataframe with the condensed items expanded to long format
    """
    if splitcol2sep==None:
        return input_df

    for split_col, separator in splitcol2sep.items():
        idx_name = f"{split_col}_idxs"
        split_col_series = input_df[split_col].str.split(separator)
        input_df = input_df.drop(columns = [split_col])

        input_df[idx_name] = [list(range(len(x))) for x in split_col_series]
        exploded_input = input_df.explode(idx_name)
        exploded_split_col_series = split_col_series.explode()

        exploded_input[split_col] = exploded_split_col_series.replace('', 0) #the column with the intensities has to come after to column with the idxs

        exploded_input = exploded_input.astype({split_col: float})
        exploded_input = exploded_input[exploded_input[split_col]>value_threshold]
        #exploded_input = exploded_input.rename(columns = {'var1': split_col})
    return exploded_input



def add_merged_ionnames(df_subset, ion_hierarchy_local, ion_headers_grouped, quant_id_dict, hierarchy_type):
    """puts together the hierarchical ion names as a column in a given input dataframe"""
    all_ion_headers = list(itertools.chain.from_iterable(ion_headers_grouped))
    columns_to_index = [x for x in df_subset.columns if x not in all_ion_headers]
    df_subset = df_subset.set_index(columns_to_index)

    rows = df_subset[all_ion_headers].to_numpy()
    ions = []

    for row in rows: #iterate through dataframe
        count = 0
        ionstring = ""
        for lvl_idx in range(len(ion_hierarchy_local)):
            ionstring += f"{ion_hierarchy_local[lvl_idx]}"
            for sublvl in ion_headers_grouped[lvl_idx]:
                ionstring+= f"_{row[count]}_"
                count+=1
        ions.append(ionstring)
    df_subset[config.QUANT_ID] = ions
    df_subset = df_subset.reset_index()
    if quant_id_dict!= None:
        df_subset = df_subset.rename(columns = {quant_id_dict.get(hierarchy_type) : "quant_val"})
    return df_subset

# %% ../nbdev_nbs/04_utils.ipynb 34
import os.path
def reformat_and_write_longtable_according_to_config(input_file, outfile_name, config_dict_for_type, sep = "\t",decimal = ".", enforce_largefile_processing = False, chunksize =1000_000):
    """Reshape a long format proteomics results table (e.g. Spectronaut or DIA-NN) to a wide format table.
    :param file input_file: long format proteomic results table
    :param string input_type: the configuration key stored in the config file (e.g. "diann_precursor")
    """
    filesize = os.path.getsize(input_file)/(1024**3) #size in gigabyte
    file_is_large = (filesize>10 and str(input_file).endswith(".zip")) or filesize>50 or enforce_largefile_processing

    if file_is_large:
        tmpfile_large = f"{input_file}.tmp.longformat.columnfilt.tsv" #only needed when file is large
        #remove potential leftovers from previous processings
        if os.path.exists(tmpfile_large):
            os.remove(tmpfile_large)
        if os.path.exists(outfile_name):
            os.remove(outfile_name)
    
    relevant_cols = get_relevant_columns_config_dict(config_dict_for_type)
    input_df_it = utils_fileread.read_file_with_pandas(input_file=input_file, sep=sep, decimal=decimal, usecols=relevant_cols, chunksize=chunksize)
    input_df_list = []
    header = True
    for input_df_subset in input_df_it:
        input_df_subset = adapt_subtable(input_df_subset, config_dict_for_type)
        if file_is_large:
            write_chunk_to_file(input_df_subset,tmpfile_large, header)
        else:
            input_df_list.append(input_df_subset)
        header = False
        
    if file_is_large:
        process_with_dask(tmpfile_columnfilt=tmpfile_large , outfile_name = outfile_name, config_dict_for_type=config_dict_for_type)
    else:
        input_df = pd.concat(input_df_list)
        input_reshaped = reshape_input_df(input_df, config_dict_for_type)
        input_reshaped.to_csv(outfile_name, sep = "\t", index = None)
    

def adapt_subtable(input_df_subset, config_dict):
    input_df_subset = filter_input(config_dict.get("filters", {}), input_df_subset)
    if "ion_hierarchy" in config_dict.keys():
        return merge_protein_cols_and_ion_dict(input_df_subset, config_dict)
    else:
        return merge_protein_and_ion_cols(input_df_subset, config_dict)

# %% ../nbdev_nbs/04_utils.ipynb 35
import dask.dataframe as dd
import pandas as pd
import glob
import os
import shutil 

def process_with_dask(*, tmpfile_columnfilt, outfile_name, config_dict_for_type):
    df = dd.read_csv(tmpfile_columnfilt, sep = "\t")
    allcols = df[config_dict_for_type.get("sample_ID")].drop_duplicates().compute() # the columns of the output table are the sample IDs
    allcols = extend_sample_allcolumns_for_plexdia_case(allcols_samples=allcols, config_dict_for_type=config_dict_for_type)
    allcols = [config.PROTEIN_ID, config.QUANT_ID] + sorted(allcols)
    df = df.set_index(config.PROTEIN_ID)
    sorted_filedir = f"{tmpfile_columnfilt}_sorted"
    df.to_csv(sorted_filedir, sep = "\t")
    #now the files are sorted and can be pivoted chunkwise (multiindex pivoting at the moment not possible in dask)
    files_dask = glob.glob(f"{sorted_filedir}/*part")
    header = True
    for file in files_dask:
        input_df = pd.read_csv(file, sep = "\t")
        if len(input_df.index) <2:
            continue
        input_reshaped = reshape_input_df(input_df, config_dict_for_type)
        input_reshaped = sort_and_add_columns(input_reshaped, allcols)
        write_chunk_to_file(input_reshaped, outfile_name, header)
        header = False
    os.remove(tmpfile_columnfilt)
    shutil.rmtree(sorted_filedir)

def reshape_input_df(input_df, config_dict):
    input_df = input_df.astype({'quant_val': 'float'})
    input_df = adapt_input_df_columns_in_case_of_plexDIA(input_df=input_df, config_dict_for_type=config_dict)
    input_reshaped = pd.pivot_table(input_df, index = [config.PROTEIN_ID, config.QUANT_ID], columns = config_dict.get("sample_ID"), values = 'quant_val', fill_value=0)

    input_reshaped = input_reshaped.reset_index()
    return input_reshaped


def sort_and_add_columns(input_reshaped, allcols):
    missing_cols = set(allcols) - set(input_reshaped.columns)
    input_reshaped[list(missing_cols)] = 0
    input_reshaped = input_reshaped[allcols]
    return input_reshaped


def extend_sample_allcolumns_for_plexdia_case(allcols_samples, config_dict_for_type):
    if is_plexDIA_table(config_dict_for_type):
        new_allcols = []
        channels = ['mTRAQ-n-0', 'mTRAQ-n-4', 'mTRAQ-n-8']
        for channel in channels:
            for sample in allcols_samples:
                new_allcols.append(merge_channel_and_sample_string(sample, channel))
        return new_allcols
    else:
        return allcols_samples

# %% ../nbdev_nbs/04_utils.ipynb 36
#PLEXDIA case

def adapt_input_df_columns_in_case_of_plexDIA(input_df,config_dict_for_type):
    if is_plexDIA_table(config_dict_for_type):
        input_df = extend_sampleID_column_for_plexDIA_case(input_df, config_dict_for_type)
        input_df = set_mtraq_reduced_ion_column_into_dataframe(input_df)
        return input_df
    else:
        return input_df


def extend_sampleID_column_for_plexDIA_case(input_df,config_dict_for_type):
    channels_per_peptide = parse_channel_from_peptide_column(input_df)
    return merge_sample_id_and_channels(input_df, channels_per_peptide, config_dict_for_type)


def set_mtraq_reduced_ion_column_into_dataframe(input_df):
    new_ions = remove_mtraq_modifications_from_ion_ids(input_df[config.QUANT_ID])
    input_df[config.QUANT_ID] = new_ions
    return input_df

def remove_mtraq_modifications_from_ion_ids(ions):
    new_ions = []
    all_mtraq_tags = ["(mTRAQ-K-0)", "(mTRAQ-K-4)", "(mTRAQ-K-8)", "(mTRAQ-n-0)", "(mTRAQ-n-4)", "(mTRAQ-n-8)"]
    for ion in ions:
        for tag in all_mtraq_tags:
            ion = ion.replace(tag, "")
        new_ions.append(ion)
    return new_ions


def is_plexDIA_table(config_dict_for_type):
    return config_dict_for_type.get('channel_ID') == ['Channel.0', 'Channel.4', 'Channel.8']


import re
def parse_channel_from_peptide_column(input_df):
    channels = []
    for pep in input_df['Modified.Sequence']:
        pattern = "(.*)(\(mTRAQ-n-.\))(.*)"
        matched = re.match(pattern, pep)
        num_appearances = pep.count("mTRAQ-n-")
        if matched and num_appearances==1:
            channels.append(matched.group(2))
        else:
            channels.append("NA")
    return channels

def merge_sample_id_and_channels(input_df, channels, config_dict_for_type):
    sample_id = config_dict_for_type.get("sample_ID")
    sample_ids = list(input_df[sample_id])
    input_df[sample_id] = [merge_channel_and_sample_string(sample_ids[idx], channels[idx]) for idx in range(len(sample_ids))]
    return input_df
            
def merge_channel_and_sample_string(sample, channel):
    return f"{sample}_{channel}"

# %% ../nbdev_nbs/04_utils.ipynb 38
def reformat_and_write_wideformat_table(peptides_tsv, outfile_name, config_dict):
    input_df = pd.read_csv(peptides_tsv,sep="\t", encoding ='latin1')
    filter_dict = config_dict.get("filters")
    protein_cols = config_dict.get("protein_cols")
    ion_cols = config_dict.get("ion_cols")
    input_df = filter_input(filter_dict, input_df)
    #input_df = merge_protein_and_ion_cols(input_df, config_dict)
    input_df = merge_protein_cols_and_ion_dict(input_df, config_dict)
    if 'quant_pre_or_suffix' in config_dict.keys():
        quant_pre_or_suffix = config_dict.get('quant_pre_or_suffix')
        headers = [config.PROTEIN_ID, config.QUANT_ID] + list(filter(lambda x: x.startswith(quant_pre_or_suffix) or x.endswith(quant_pre_or_suffix), input_df.columns))
        input_df = input_df[headers]
        input_df = input_df.rename(columns = lambda x : x.replace(quant_pre_or_suffix, ""))

    #input_df = input_df.reset_index()
    
    input_df.to_csv(outfile_name, sep = '\t', index = None)


import os
def check_for_processed_runs_in_results_folder(results_folder):
    contained_condpairs = []
    folder_files = os.listdir(results_folder)
    result_files = list(filter(lambda x: "results.tsv" in x ,folder_files))
    for result_file in result_files:
        res_name = result_file.replace(".results.tsv", "")
        if ((f"{res_name}.normed.tsv" in folder_files) & (f"{res_name}.results.ions.tsv" in folder_files)):
            contained_condpairs.append(res_name)
    return contained_condpairs





import pandas as pd
import os
import pathlib

def import_data(input_file, input_type_to_use = None, samples_subset = None, filter_dict = None):
    """
    Function to import peptide level data. Depending on available columns in the provided file,
    the function identifies the type of input used (e.g. Spectronaut, MaxQuant, DIA-NN), reformats if necessary
    and returns a generic wide-format dataframe
    :param file input_file: quantified peptide/ion -level data
    :param file results_folder: the folder where the directlfq outputs are stored
    """

    samples_subset = add_ion_protein_headers_if_applicable(samples_subset)
    file_is_already_formatted = ("aq_reformat" in input_file) | (input_type_to_use == "directlfq")
    if file_is_already_formatted:
        file_to_read = input_file
    else:
        file_to_read = reformat_and_save_input_file(input_file=input_file, input_type_to_use=input_type_to_use, filter_dict=filter_dict)
    
    input_reshaped = pd.read_csv(file_to_read, sep = "\t", encoding = 'latin1', usecols=samples_subset)
    input_reshaped = adapt_table_for_alphabaseformat_backward_compatibility(file_is_already_formatted, input_reshaped)
    input_reshaped = input_reshaped.drop_duplicates(subset=config.QUANT_ID)

    return input_reshaped

def add_ion_protein_headers_if_applicable(samples_subset):
    if samples_subset is not None:
        return samples_subset + [config.QUANT_ID, config.PROTEIN_ID]
    else:
        return None


def reformat_and_save_input_file(input_file, input_type_to_use = None, filter_dict = None):
    
    input_type, config_dict_for_type, sep = get_input_type_and_config_dict(input_file, input_type_to_use)

    if filter_dict is not None:
        config_dict_for_type['filters']=  dict(config_dict_for_type.get('filters', {}),**filter_dict)
    LOGGER.info(f"using input type {input_type}")
    format = config_dict_for_type.get('format')
    outfile_name = f"{input_file}.{input_type}.aq_reformat.tsv"

    if format == "longtable":
        reformat_and_write_longtable_according_to_config(input_file, outfile_name,config_dict_for_type, sep = sep)
    elif format == "widetable":
        reformat_and_write_wideformat_table(input_file, outfile_name, config_dict_for_type)
    else:
        raise Exception('Format not recognized!')
    return outfile_name


def adapt_table_for_alphabaseformat_backward_compatibility(file_is_already_formatted, input_reshaped):
    if file_is_already_formatted:
        if 'quant_id' in input_reshaped.columns:
            return input_reshaped.rename(columns = {'quant_id' : config.QUANT_ID})

    return input_reshaped





# %% ../nbdev_nbs/04_utils.ipynb 45
import pandas as pd
import os.path
import pathlib

def get_input_type_and_config_dict(input_file, input_type_to_use = None):
    #parse the type of input (e.g. Spectronaut Fragion+MS1Iso) out of the input file


    config_dict = load_config(INTABLE_CONFIG)
    type2relevant_columns = get_type2relevant_cols(config_dict)

    if "aq_reformat.tsv" in input_file:
        input_file = get_original_file_from_aq_reformat(input_file)

    filename = str(input_file)
    if '.csv' in filename:
        sep=','
    if '.tsv' in filename:
        sep='\t'
    if '.txt' in filename:
        sep='\t'
    else:
        sep="\t"


    uploaded_data_columns = utils_fileread.read_columns_from_file(input_file, sep=sep)

    for input_type in type2relevant_columns.keys():
        if (input_type_to_use is not None) and (input_type!=input_type_to_use):
            continue
        relevant_columns = type2relevant_columns.get(input_type)
        relevant_columns = [x for x in relevant_columns if x] #filter None values
        if set(relevant_columns).issubset(uploaded_data_columns):
            config_dict_type =  config_dict.get(input_type)
            return input_type, config_dict_type, sep
    raise TypeError("format not specified in intable_config.yaml!")

import re
def get_original_file_from_aq_reformat(input_file):
    matched = re.match("(.*)(\..*\.)(aq_reformat\.tsv)",input_file)
    return matched.group(1)

# %% ../nbdev_nbs/04_utils.ipynb 47
def import_config_dict():
    config_dict = load_config(INTABLE_CONFIG)
    return config_dict

# %% ../nbdev_nbs/04_utils.ipynb 48
import pandas as pd

def load_samplemap(samplemap_file):
    file_ext = os.path.splitext(samplemap_file)[-1]
    if file_ext=='.csv':
        sep=','
    if (file_ext=='.tsv') | (file_ext=='.txt'):
        sep='\t'

    if 'sep' not in locals():
        LOGGER.info(f"neither of the file extensions (.tsv, .csv, .txt) detected for file {samplemap_file}! Trying with tab separation. In the case that it fails, please add the appropriate extension to your file name.")
        sep = "\t"

    return pd.read_csv(samplemap_file, sep = sep, encoding ='latin1', dtype='str')

# %% ../nbdev_nbs/04_utils.ipynb 49
def prepare_loaded_tables(data_df, samplemap_df):
    """
    Integrates information from the peptide/ion data and the samplemap, selects the relevant columns and log2 transforms intensities.
    """
    samplemap_df = samplemap_df[samplemap_df["condition"]!=""] #remove rows that have no condition entry
    filtvec_not_in_data = [(x in data_df.columns) for x in samplemap_df["sample"]] #remove samples that are not in the dataframe
    samplemap_df = samplemap_df[filtvec_not_in_data]
    headers = [config.PROTEIN_ID] + samplemap_df["sample"].to_list()
    data_df = data_df.set_index(config.QUANT_ID)
    for sample in samplemap_df["sample"]:
        data_df[sample] = np.log2(data_df[sample].replace(0, np.nan))
    return data_df[headers], samplemap_df

# %% ../nbdev_nbs/04_utils.ipynb 50
class LongTableReformater():
    """Generic class to reformat tabular files in chunks. For the specific cases you can inherit the class and specify reformat and iterate function
    """
    def __init__(self, input_file):
        self._input_file = input_file
        self._reformatting_function = None
        self._iterator_function = self.__initialize_df_iterator__
        self._concat_list = []

    def reformat_and_load_acquisition_data_frame(self):

        input_df_it = self._iterator_function()
        
        input_df_list = []
        for input_df_subset in input_df_it:
            input_df_subset = self._reformatting_function(input_df_subset)
            input_df_list.append(input_df_subset)
        input_df = pd.concat(input_df_list)
        
        return input_df

    def reformat_and_save_acquisition_data_frame(self, output_file):
        
        input_df_it = self._iterator_function()
        write_header = True
        
        for input_df_subset in input_df_it:
            input_df_subset = self._reformatting_function(input_df_subset)
            self.__write_reformatted_df_to_file__(input_df_subset, output_file, write_header)
            write_header = False

    def __initialize_df_iterator__(self):
        return pd.read_csv(self._input_file, sep = "\t", encoding ='latin1', chunksize=1000000)
    
    @staticmethod
    def __write_reformatted_df_to_file__(reformatted_df, filepath ,write_header):
        reformatted_df.to_csv(filepath, header=write_header, mode='a', sep = "\t", index = None)

