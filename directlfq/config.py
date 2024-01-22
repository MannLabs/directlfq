import logging


def setup_logging():
    logging.basicConfig(level=logging.INFO, 
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

##########################
LOG_PROCESSED_PROTEINS = True

def set_log_processed_proteins(log_processed_proteins = True):
    global LOG_PROCESSED_PROTEINS
    LOG_PROCESSED_PROTEINS = log_processed_proteins


##########################
PROTEIN_ID = 'protein'
QUANT_ID = 'ion'

def set_global_protein_and_ion_id(protein_id = 'protein', quant_id = 'ion'):
    global PROTEIN_ID
    global QUANT_ID
    PROTEIN_ID = protein_id
    QUANT_ID = quant_id

##########################
COMPILE_NORMALIZED_ION_TABLE = True

def set_compile_normalized_ion_table(compile_normalized_ion_table = True):
    global COMPILE_NORMALIZED_ION_TABLE
    COMPILE_NORMALIZED_ION_TABLE = compile_normalized_ion_table

