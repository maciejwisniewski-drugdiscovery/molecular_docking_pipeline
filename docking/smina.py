from pathlib import Path
import pandas as pd
import logging
import subprocess
import os
os.environ['DATA_DIR'] = '/Users/maciejwisniewski/Projects/ingenix/molecular_docking/data'
os.environ['DOCKING_RESULTS_DIR'] = '/Users/maciejwisniewski/Projects/ingenix/molecular_docking/docking_results'
os.environ['LOG_DIR'] = '/Users/maciejwisniewski/Projects/ingenix/molecular_docking/log'
os.environ['LIGANDS_DIR'] = '/Users/maciejwisniewski/Projects/ingenix/molecular_docking/data/ligands'
os.environ['NATIVE_LIGANDS_DIR'] = '/Users/maciejwisniewski/Projects/ingenix/molecular_docking/data/native_ligands'
os.environ['PROTEINS_DIR'] = '/Users/maciejwisniewski/Projects/ingenix/molecular_docking/data/proteins'


def check_files_in_directory(directory_path):
    # Define the directory path
    dir_path = Path(directory_path)

    # Define the required files
    required_files = ["poses.pdbqt", "log.txt", "atom_terms.txt"]

    # Check if each file exists and is non-empty
    for file_name in required_files:
        file_path = dir_path / file_name
        if not file_path.is_file():
            print(f"File {file_name} is missing.")
            return False
        elif file_path.stat().st_size == 0:
            print(f"File {file_name} is empty.")
            return False

    print("All required files are present and non-empty.")
    return True
def smina_docking(protein_filepath: Path, native_ligand_filepath: Path, ligand_filepath: Path, output_dirpath: Path):

    pdbqt_output_filepath = output_dirpath / 'poses.pdbqt'
    atom_terms_output_filepath = output_dirpath / 'atom_terms.txt'
    log_output_filepath = output_dirpath / 'log.txt'

    smina_command = ['smina',
                     '-r', str(protein_filepath),
                     '-l', str(ligand_filepath),
                     '--autobox_ligand', str(native_ligand_filepath),
                     '--autobox_add', '4',
                     '--exhaustiveness', '32',
                     '--num_modes', '10',
                     '-o', str(pdbqt_output_filepath),
                     '--atom_terms', str(atom_terms_output_filepath),
                     '--log', str(log_output_filepath),
                     '--atom_term_data',
                     '--cpu', '4',
                     ]
    docking = subprocess.run(smina_command, shell=False, capture_output=True, text=True)


def smina_based_docking_iterator_on_diagonal(complexes: pd.DataFrame):
    for _, complex in complexes.iterrows():
        protein_pdbqt_filepath = Path(complex['protein_processed_pdbqt_filepath'])
        ligand_pdbqt_filepath = Path(complex['ligand_processed_pdbqt_filepath'])
        native_ligand_pdbqt_filepath = Path(complex['native_ligand_processed_pdbqt_filepath'])

        protein_id = protein_pdbqt_filepath.stem
        ligand_id = ligand_pdbqt_filepath.stem
        native_ligand_id = native_ligand_pdbqt_filepath.stem

        docking_results_id = protein_id + '___' + ligand_id
        output_dirpath = Path(os.getenv('DOCKING_RESULTS_DIR')) / docking_results_id
        output_dirpath.mkdir(parents=True, exist_ok=True)
        logging.info(f'Perform Docking for {docking_results_id}')
        smina_docking(protein_filepath=protein_pdbqt_filepath,
                      ligand_filepath=ligand_pdbqt_filepath,
                      native_ligand_filepath=native_ligand_pdbqt_filepath,
                      output_dirpath=output_dirpath)

        if check_files_in_directory(output_dirpath):
            logging.info(f'Docking successful for {docking_results_id}')
        else:
            logging.info(f'Docking failed for {docking_results_id}')
def smina_based_docking_iterator_off_diagonal(complexes: pd.DataFrame):
    for _, protein_info in complexes.iterrows():
        for _, ligand_info in complexes.iterrows():
            protein_pdbqt_filepath = Path(protein_info['protein_processed_pdbqt_filepath'])
            ligand_pdbqt_filepath = Path(ligand_info['ligand_processed_pdbqt_filepath'])
            native_ligand_pdbqt_filepath = Path(protein_info['native_ligand_processed_pdbqt_filepath'])

            protein_id = protein_pdbqt_filepath.stem
            ligand_id = ligand_pdbqt_filepath.stem
            native_ligand_id = native_ligand_pdbqt_filepath.stem

            docking_results_id = protein_id + '___' + ligand_id
            output_dirpath = Path(os.getenv('DOCKING_RESULTS_DIR')) / docking_results_id
            output_dirpath.mkdir(parents=True, exist_ok=True)
            logging.info(f'Perform Docking for {docking_results_id}')
            smina_docking(protein_filepath=protein_pdbqt_filepath,
                          ligand_filepath=ligand_pdbqt_filepath,
                          native_ligand_filepath=native_ligand_pdbqt_filepath,
                          output_dirpath=output_dirpath)

            if check_files_in_directory(output_dirpath):
                logging.info(f'Docking successful for {docking_results_id}')
            else:
                logging.info(f'Docking failed for {docking_results_id}')

def smina_based_docking_iterator(complexes: pd.DataFrame, is_off_diagonal: bool = True):
    if is_off_diagonal == False:
        smina_based_docking_iterator_on_diagonal(complexes)
    elif is_off_diagonal == True:
        smina_based_docking_iterator_off_diagonal(complexes)
