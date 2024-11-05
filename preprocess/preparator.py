import pandas as pd
import rdkit
from openbabel import pybel
from pathlib import Path
from collections import defaultdict
import logging
import os
os.environ['DATA_DIR'] = '/Users/maciejwisniewski/Projects/ingenix/molecular_docking/data'
os.environ['LOG_DIR'] = '/Users/maciejwisniewski/Projects/ingenix/molecular_docking/log'
os.environ['LIGANDS_DIR'] = '/Users/maciejwisniewski/Projects/ingenix/molecular_docking/data/ligands'
os.environ['NATIVE_LIGANDS_DIR'] = '/Users/maciejwisniewski/Projects/ingenix/molecular_docking/data/native_ligands'
os.environ['PROTEINS_DIR'] = '/Users/maciejwisniewski/Projects/ingenix/molecular_docking/data/proteins'

def openbabel_based_preprocess_ligand(raw: Path, processed_sdf: Path, processed_pdbqt: Path, pH: float = 7.4):
    # Load Ligand to Pybel
    pybel_ligand = next(pybel.readfile('sdf',str(raw)))
    pybel_ligand.OBMol.CorrectForPH(pH)
    pybel_ligand.addh()
    for atom in pybel_ligand.atoms:
        atom.OBAtom.GetPartialCharge()
    a=1
    pybel_ligand.write("sdf", str(processed_sdf), overwrite=True)
    pybel_ligand.write("pdbqt", str(processed_pdbqt), overwrite=True)
def openbabel_based_preprocess_protein(raw: Path,processed_pdb: Path,processed_pdbqt: Path, pH: float = 7.4):
    # Load Protein to Pybel
    pybel_protein = next(pybel.readfile('pdb',str(raw)))
    pybel_protein.OBMol.CorrectForPH(pH)
    pybel_protein.addh()
    pybel_protein.write('pdb',str(processed_pdb), overwrite=True)
    for i, atom in enumerate(pybel_protein.atoms):
        print(i)
        atom.OBAtom.GetPartialCharge()
    pybel_protein.write("pdbqt", str(processed_pdbqt), overwrite=True)

def openbabel_based_preprocessing_iterator(complexes: pd.DataFrame, create_df: bool = True):
    logging.info(f'Initializing OpenBabel Preprocessing Iterator!')
    new_dict = defaultdict(list)
    for _, complex in complexes[:3].iterrows():
        pdb_id = complex['pdb_id'].lower()
        ligand_ccd_code = complex['ligand_id']

        protein_raw_pdb_filepath = Path(complex['protein_raw_pdb_filepath'])
        ligand_raw_sdf_filepath = Path(complex['ligand_raw_sdf_filepath'])
        native_ligand_raw_sdf_filepath = Path(complex['native_ligand_raw_sdf_filepath'])

        protein_id = protein_raw_pdb_filepath.stem
        ligand_id = ligand_raw_sdf_filepath.stem
        native_ligand_id = native_ligand_raw_sdf_filepath.stem

        logging.info(f'OpenBabel Processing Protein and Native Ligand from {native_ligand_id} and Ligand {ligand_id}')

        protein_processed_pdb_filepath = Path(os.getenv('PROTEINS_DIR')) / 'processed_pdb' / str(protein_id+'.pdb')
        protein_processed_pdbqt_filepath = Path(os.getenv('PROTEINS_DIR')) / 'processed_pdbqt' / str(protein_id+'.pdbqt')

        ligand_processed_sdf_filepath = Path(os.getenv('LIGANDS_DIR')) / 'processed_sdf' / str(ligand_id+'.sdf')
        ligand_processed_pdbqt_filepath = Path(os.getenv('LIGANDS_DIR')) / 'processed_pdbqt' / str(ligand_id + '.pdbqt')

        native_ligand_processed_sdf_filepath = Path(os.getenv('NATIVE_LIGANDS_DIR')) / 'processed_sdf' / str(native_ligand_id+'.sdf')
        native_ligand_processed_pdbqt_filepath = Path(os.getenv('NATIVE_LIGANDS_DIR')) / 'processed_pdbqt' / str(native_ligand_id+'.pdbqt')

        # Protein Preprocessing
        if (not protein_processed_pdb_filepath.exists()) & (not protein_processed_pdbqt_filepath.exists()):
            openbabel_based_preprocess_protein(raw=protein_raw_pdb_filepath,
                                               processed_pdb=protein_processed_pdb_filepath,
                                               processed_pdbqt=protein_processed_pdbqt_filepath)
        if (not protein_processed_pdb_filepath.exists()) & (not protein_processed_pdbqt_filepath.exists()):
            logging.error(f'Ligand {protein_id} preprocessing failed')

        # Ligand Preprocessing
        if (not ligand_processed_sdf_filepath.exists()) & (not ligand_processed_pdbqt_filepath.exists()):
            openbabel_based_preprocess_ligand(raw=ligand_raw_sdf_filepath,
                                              processed_sdf=ligand_processed_sdf_filepath,
                                              processed_pdbqt=ligand_processed_pdbqt_filepath)

        if (not ligand_processed_sdf_filepath.exists()) & (not ligand_processed_pdbqt_filepath.exists()):
            logging.error(f'Ligand {ligand_id} preprocessing failed')

        # Native Ligand Preprocessing
        if (not native_ligand_processed_pdbqt_filepath.exists()) and (not native_ligand_processed_sdf_filepath.exists()):
            openbabel_based_preprocess_ligand(raw=native_ligand_raw_sdf_filepath,
                                              processed_sdf=native_ligand_processed_sdf_filepath,
                                              processed_pdbqt=native_ligand_processed_pdbqt_filepath)
        if (not native_ligand_processed_pdbqt_filepath.exists()) and (not native_ligand_processed_sdf_filepath.exists()):
            logging.error(f'Native Ligand {native_ligand_id} preprocessing failed')

        new_dict['pdb_id'].append(pdb_id)
        new_dict['ligand_id'].append(ligand_ccd_code)
        new_dict['ligand_raw_sdf_filepath'].append(ligand_raw_sdf_filepath)
        new_dict['native_ligand_raw_sdf_filepath'].append(native_ligand_raw_sdf_filepath)
        new_dict['protein_raw_pdb_filepath'].append(protein_raw_pdb_filepath)
        new_dict['ligand_processed_sdf_filepath'].append(ligand_processed_sdf_filepath)
        new_dict['ligand_processed_pdbqt_filepath'].append(ligand_processed_pdbqt_filepath)
        new_dict['native_ligand_processed_sdf_filepath'].append(native_ligand_processed_sdf_filepath)
        new_dict['native_ligand_processed_pdbqt_filepath'].append(native_ligand_processed_pdbqt_filepath)
        new_dict['protein_processed_pdb_filepath'].append(protein_processed_pdb_filepath)
        new_dict['protein_processed_pdbqt_filepath'].append(protein_processed_pdbqt_filepath)

    if create_df:
        df = pd.DataFrame.from_dict(new_dict)
        new_df_filepath = Path(os.getenv('DATA_DIR')) / 'processed_complexes.csv'
        df.to_csv(new_df_filepath)
        logging.info(f'Dataframe with Processed Filepaths Created!')
        return df
    else:
        return None
