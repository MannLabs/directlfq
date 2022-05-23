import numpy as np
import pandas as pd
import itertools
import os
os.chdir(".")


def load_perseus_df():
    df_file = os.path.join("..","..","..","test_data", "ptmsite_mapping", "shortened_plugin_outputv4_loc_prob_75.txt.zip")
    perseus_res = pd.read_csv(df_file , sep = "\t")
    perseus_res = perseus_res.iloc[1:]
    perseus_gene = [x.split("_")[0] for x in perseus_res["PTM_collapse_key"]]
    perseus_start_keys = [x.split("_")[1][1:] for x in perseus_res["PTM_collapse_key"]]
    perseus_multiplicity = list(perseus_res["PTM_0_num"].astype("str"))
    perseus_res["pseudo_collapse_key"] = [f'{perseus_gene[x]}_{perseus_start_keys[x]}_{perseus_multiplicity[x]}' for x in range(len(perseus_gene))]
    return perseus_res


def load_aq_df():
    df_file = os.path.join(".", "results", "ptm_ids.tsv")
    aq_res = pd.read_csv(df_file, sep = "\t")
    aq_res["ptm_id"] = aq_res["REFPROT"] + aq_res["site"]
    aq_start_keys = [x.split(" ")[0].replace("[", "").replace("]", "") for x in aq_res["site"]]
    aq_multiplicity = [str(len(x.split(" "))) for x in aq_res["site"]]
    aq_genes = [x for x in aq_res["gene"]]
    aq_res["pseudo_collapse_key"] = [f'{aq_genes[x]}_{aq_start_keys[x]}_{aq_multiplicity[x]}' for x in range(len(aq_genes))]
    return aq_res


def get_exp_precursor_id(samplename, precursor):
    return f"{samplename}_{precursor}"


def get_key2precursors_perseus(perseus_df):
    ckey2ions = {}
    for _,row in perseus_df.iterrows():
        ckey = row["pseudo_collapse_key"]
        ions = row["EG.PrecursorId"].split(";")
        ionlist = ckey2ions.get(ckey, set())
        ionlist.update(ions)
        ckey2ions[ckey] = ionlist

    return ckey2ions

def get_key2expprecursors_perseus(perseus_res, key2prec_pers, all_expprecs):
    key2expprecs = {}
    perseus_res = perseus_res.set_index("pseudo_collapse_key").sort_index()
    experiment_columns = np.array(list(filter(lambda x: "MaTa_SA_phosphatase_screen_phospho" in x, perseus_res.columns)))
    for key_perseus in perseus_res.index:
        perseus_key_df = perseus_res.loc[key_perseus]#automatic conversion to series by pandas, if only one row exists
        if type(perseus_key_df) != type(pd.Series()): #continue if more than one row exists
            continue
        perseus_key_df_vals = perseus_key_df[experiment_columns].to_numpy(dtype = "float")
        nonan_experiments_filt = ~np.isnan(perseus_key_df_vals)
        nonan_experiments_names = experiment_columns[nonan_experiments_filt]
        precursors_perseus = list(key2prec_pers.get(key_perseus))
        exp_precs_pers = [get_exp_precursor_id(exp, prec) for exp, prec in itertools.product(nonan_experiments_names, precursors_perseus)]
        exp_precs_pers = set(exp_precs_pers).intersection(all_expprecs)
        key2expprecs[key_perseus] = list(exp_precs_pers)
    print("precursor assign finished")

    return key2expprecs


def get_key2expprecursors_aq(aq_res, all_expprecs):
    key2expprecs = {}
    aq_res = aq_res.set_index("pseudo_collapse_key").sort_index()
    for key_aq in aq_res.index.unique():
        aq_selected = aq_res.loc[[key_aq]]
        exp_precs_aq = {get_exp_precursor_id(exp.replace(".htrms", ""), prec) for exp, prec in aq_selected[["R.Label", "FG.Id"]].itertuples(index=False)}
        key2expprecs[key_aq] = list(exp_precs_aq)
    return key2expprecs

def get_collapse_df(collapse_df):
    keylist = [x for x in set(collapse_df["pseudo_collapse_key"])]
    genelist  = [x.split("_")[0] for x in keylist]
    genesite_list  = [x.split("_")[0] + "_" + x.split("_")[1] for x in keylist]
    keydf = pd.DataFrame({"key" : keylist, "genesite" : genesite_list, "gene" : genelist})
    return keydf
    
def get_exp_precursors_for_collapse_key(aq_res, key_aq, perseus_res, key_perseus, key2prec_pers, all_expprecs):
    filtvec_perseus = (perseus_res["pseudo_collapse_key"] == key_perseus)
    perseus_key_df = perseus_res[filtvec_perseus]
    experiment_columns = np.array(list(filter(lambda x: "MaTa_SA_phosphatase_screen_phospho" in x, perseus_res.columns)))
    perseus_key_df_vals = perseus_key_df[experiment_columns].to_numpy(dtype = "float")[0]
    nonan_experiments_filt = ~np.isnan(perseus_key_df_vals)
    nonan_experiments_names = experiment_columns[nonan_experiments_filt]
    precursors_perseus = list(key2prec_pers.get(key_perseus))
    exp_precs_pers = [get_exp_precursor_id(exp, prec) for exp, prec in itertools.product(nonan_experiments_names, precursors_perseus)]

    exp_precs_pers = set(exp_precs_pers).intersection(all_expprecs)

    aq_selected = aq_res[aq_res["pseudo_collapse_key"] == key_aq]
    exp_precs_aq = {get_exp_precursor_id(exp.replace(".htrms", ""), prec) for exp, prec in aq_selected[["R.Label", "FG.Id"]].itertuples(index=False)}

    return exp_precs_pers, exp_precs_aq



def get_site_probabilites_exp_precursors(shortened_aq_input, exp_precs_pers, exp_precs_aq, site_prob_column = 'EG.PTMLocalizationProbabilities'):
    filtvec_pers = [ get_exp_precursor_id(exp, prec) in exp_precs_pers for exp, prec in shortened_aq_input[["R.FileName", "FG.Id"]].itertuples(index=False)]
    filtvec_aq = [ get_exp_precursor_id(exp, prec) in exp_precs_aq for exp, prec in shortened_aq_input[["R.FileName", "FG.Id"]].itertuples(index=False)]

    result_perseus = shortened_aq_input[filtvec_pers][site_prob_column]
    result_aq = shortened_aq_input[filtvec_aq][site_prob_column]
    return result_perseus, result_aq


def get_shortened_aq_input():
    df_file = os.path.join("..","..","..","test_data", "ptmsite_mapping", "shortened_aq.tsv.zip")
    print(os.path.abspath(df_file))
    shortened_aq_input = pd.read_csv(df_file, sep = "\t")
    return shortened_aq_input

def get_samplemap():
    samplemap_file = os.path.join("..","..","..","test_data", "ptmsite_mapping", "samples.map")
    samplemap_df = pd.read_csv(samplemap_file, sep = "\t")
    return samplemap_df

def get_all_expprecs(shortened_aq_input):
    all_expprecs = {get_exp_precursor_id(exp, prec) for exp, prec in shortened_aq_input[["R.FileName", "FG.Id"]].itertuples(index=False)}
    return all_expprecs
