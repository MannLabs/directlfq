import logging
import pandas as pd


def setup_logging():
    logging.basicConfig(level=logging.INFO, 
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
setup_logging()
##########################
LOG_PROCESSED_PROTEINS = True
LOG_PROCESSED_PROTEINS_PERIOD = 100

def set_log_processed_proteins(log_processed_proteins: bool = True, period: int = 100) -> None:
    global LOG_PROCESSED_PROTEINS, LOG_PROCESSED_PROTEINS_PERIOD
    LOG_PROCESSED_PROTEINS = log_processed_proteins
    LOG_PROCESSED_PROTEINS_PERIOD = period


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

##########################
COPY_NUMPY_ARRAYS_DERIVED_FROM_PANDAS = False

def check_wether_to_copy_numpy_arrays_derived_from_pandas():
    global COPY_NUMPY_ARRAYS_DERIVED_FROM_PANDAS
    try:
        _manipulate_numpy_array_without_copy()
        COPY_NUMPY_ARRAYS_DERIVED_FROM_PANDAS = False
    except:
        logging.info('Some numpy arrays derived from pandas will be copied.')
        COPY_NUMPY_ARRAYS_DERIVED_FROM_PANDAS = True

def _manipulate_numpy_array_without_copy():
    
    protein_profile_df = pd.DataFrame({
    'ProteinA': [10, 20, 30, 40],
    'ProteinB': [15, 25, 35, 45],
    'ProteinC': [20, 30, 40, 50]
    }, index=['Sample1', 'Sample2', 'Sample3', 'Sample4'])
    
    protein_profile_df = protein_profile_df.iloc[1:3]
    protein_profile_numpy = protein_profile_df.to_numpy(copy=False)

    protein_profile_numpy[0] = protein_profile_numpy[0] +2
