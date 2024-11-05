import os
import pandas as pd
import logging
from pathlib import Path
from preprocess.plinder_preparator import plinder_based_native_complexes_files_iterator
from preprocess.preparator import openbabel_based_preprocessing_iterator
from docking.smina import smina_based_docking_iterator

if __name__ == "__main__":
    # Initialize logger
    logging_file = Path(os.getenv('LOG_DIR')) / 'plinder.log'
    logging.basicConfig(filename=logging_file, filemode="w", level=logging.INFO,
                        format="%(asctime)s - %(levelname)s - %(message)s")

    # Get Dataframe with complexes
    complexes_csv_filepath = Path(os.getenv('DATA_DIR')) / 'complexes.csv'
    if complexes_csv_filepath.exists():
        complexes_df = pd.read_csv(complexes_csv_filepath,index_col=0)
        logging.info("Complexes Dataframe Loaded")
    else:
        logging.error(f"Complexes Dataframe Filepath {complexes_csv_filepath} doesn't exists!")

    # Iterator to get raw files
    complexes_df = plinder_based_native_complexes_files_iterator(complexes = complexes_df, create_df = True)

    # Iterator to prepare Proteins and Ligands for Docking,
    complexes_df = openbabel_based_preprocessing_iterator(complexes = complexes_df, create_df = True)

    # Molecular Docking
    smina_based_docking_iterator(complexes = complexes_df, create_df = True)