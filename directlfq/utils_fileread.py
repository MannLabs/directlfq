import os
import pathlib
import logging
import pyarrow.parquet
import pandas as pd

from constants import REQUIRED_D210_COLUMNS, D210_QUAN_COLUMNS

if "__file__" in globals():
    INTABLE_CONFIG = os.path.join(pathlib.Path(__file__).parent.absolute(), "configs", "intable_config.yaml") #the yaml config is located one directory below the python library files
    CONFIG_PATH = os.path.join(pathlib.Path(__file__).parent.absolute(), "configs")

LOGGER = logging.getLogger(__name__)

def check_parquet_database(
    input_file: str,
) -> bool:
    """Check parquet database
    
    If input_file contains '.parquet', or if input_file is a directory containing .parquet files, it is treated as a parquet database.

    Args:
        input_file (str): Path to the input file or directory.

    Returns:
        bool: True if the input is a parquet database, False otherwise.
    
    """
    if input_file.endswith(".parquet"):
        return True
    elif os.path.isdir(input_file):
        for file in os.listdir(input_file):
            if file.endswith(".parquet"):
                return True
    return False

def read_file_with_pandas(input_file, decimal='.', usecols=None, chunksize=None, sep = None):
    filename = str(input_file)
    if check_parquet_database(filename):
        return read_parquet_file(input_file, usecols=usecols, chunksize=chunksize)
    else:
        if sep is None:
            if '.csv' in filename:
                sep=','
            elif '.tsv' in filename:
                sep='\t'
            else:
                sep='\t'
            LOGGER.info(f"neither of the file extensions (.tsv, .csv) detected for file {input_file}! Trying with tab separation. In the case that it fails, please provide the correct file extension")
        return pd.read_csv(input_file,sep=sep, decimal=decimal, usecols=usecols, encoding='latin1', chunksize=chunksize)


def read_parquet_file(input_file, usecols=None, chunksize=None):
    if chunksize is not None:
        return read_parque_file_chunkwise(input_file, usecols=usecols, chunksize=chunksize)
    else:
        return pd.read_parquet(input_file, columns=usecols)

def read_parque_file_chunkwise(input_file, usecols=None, chunksize=None):
    parquet_file = pyarrow.parquet.ParquetFile(input_file)
    for batch in parquet_file.iter_batches(columns=usecols, batch_size=chunksize):
        yield batch.to_pandas()


def read_columns_from_file(file, sep="\t"):
    if file.endswith(".parquet"):
        parquet_file = pyarrow.parquet.ParquetFile(file)
        return parquet_file.schema.names
    elif os.path.isdir(file):
        # If it's a directory, read schema from first parquet file
        for filename in os.listdir(file):
            if filename.endswith(".parquet"):
                first_parquet = os.path.join(file, filename)
                parquet_file = pyarrow.parquet.ParquetFile(first_parquet)
                return parquet_file.schema.names
        # If no parquet files found, fall back to CSV reading (will fail appropriately)
        return pd.read_csv(file, sep=sep, nrows=1).columns.tolist()
    else:
        return pd.read_csv(file, sep=sep, nrows=1).columns.tolist()