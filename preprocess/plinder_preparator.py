import os

os.environ['DATA_DIR'] = '/Users/maciejwisniewski/Projects/ingenix/molecular_docking/data'
os.environ['LOG_DIR'] = '/Users/maciejwisniewski/Projects/ingenix/molecular_docking/log'
os.environ['LIGANDS_DIR'] = '/Users/maciejwisniewski/Projects/ingenix/molecular_docking/data/ligands'
os.environ['NATIVE_LIGANDS_DIR'] = '/Users/maciejwisniewski/Projects/ingenix/molecular_docking/data/native_ligands'
os.environ['PROTEINS_DIR'] = '/Users/maciejwisniewski/Projects/ingenix/molecular_docking/data/proteins'

import logging
import shutil
from pathlib import Path
from typing import Optional

import pandas as pd
from tqdm import tqdm

from plinder.core.scores import query_index
from plinder.core import PlinderSystem
from rdkit import Chem
from rdkit.Chem import Mol
from collections import defaultdict

def get_full_plindex() -> pd.DataFrame:
    columns_to_query = ["system_id", "entry_pdb_id", "system_biounit_id", "entry_release_date",
                        "ligand_smiles", "ligand_rdkit_canonical_smiles", "ligand_ccd_code",
                        "system_num_interactions", "ligand_plip_type", "ligand_asym_id",
                        "ligand_instance", "ligand_instance_chain", "system_ligand_chains_asym_id",
                        "system_protein_chains_asym_id", "system_num_ligand_chains", "system_num_protein_chains",
                        "ligand_num_resolved_heavy_atoms", "ligand_num_unresolved_heavy_atoms"]
    plindex = query_index(columns=columns_to_query,
                          filters=[("ligand_is_proper", "==", True)],  # Filter to only include proper ligands
                          splits=["*"])
    return plindex


def save_rdkit_mol_to_sdf(mol: Mol, filepath: Path):
    writer = Chem.SDWriter(filepath)
    writer.write(mol)
    writer.close()


def copy_protein_file(old_filepath: Path, new_filepath: Path):
    shutil.copy(old_filepath, new_filepath)


def plinder_based_native_complex_files_generator(system_id: str, ligand_id: str, return_filepaths: bool = True):
    # Load PlinderSystem
    plinder_system = PlinderSystem(system_id=system_id)

    plinder_protein_pdb_filepath = Path(plinder_system.receptor_pdb)
    holo_structure = plinder_system.holo_structure
    native_rdkit_mol = holo_structure.resolved_ligand_mols[ligand_id]
    template_rdkit_mol = holo_structure.input_ligand_templates[ligand_id]

    # Save for Docking
    protein_id_to_save = '__'.join(system_id.split('__')[:-1])
    native_ligand_id_to_save = system_id
    template_ligand_id_to_save = None
    for ligand in plinder_system.system['ligands']:
        a=1
        if (str(ligand['instance']) == ligand_id.split('.')[0]) & (ligand['asym_id'] == ligand_id.split('.')[1]):
            a=1
            template_ligand_id_to_save = ligand['ccd_code']
    if template_ligand_id_to_save is None:
        a=1
        template_ligand_id_to_save = system_id.split('__')[-1]

    raw_protein_pdb_filepath = Path(os.getenv('PROTEINS_DIR')) / 'raw' / str(protein_id_to_save+'.pdb')
    raw_native_mol_sdf_filepath = Path(os.getenv('NATIVE_LIGANDS_DIR')) / 'raw' / str(native_ligand_id_to_save + '.sdf')
    raw_template_mol_sdf_filepath = Path(os.getenv('LIGANDS_DIR')) / 'raw' / str(template_ligand_id_to_save + '.sdf')

    if not raw_template_mol_sdf_filepath.exists():
        save_rdkit_mol_to_sdf(mol=template_rdkit_mol, filepath=raw_template_mol_sdf_filepath)
    if not raw_native_mol_sdf_filepath.exists():
        save_rdkit_mol_to_sdf(mol=native_rdkit_mol, filepath=raw_native_mol_sdf_filepath)
    if not raw_protein_pdb_filepath.exists():
        copy_protein_file(old_filepath=plinder_protein_pdb_filepath, new_filepath=raw_protein_pdb_filepath)

    if return_filepaths:
        return raw_template_mol_sdf_filepath, raw_native_mol_sdf_filepath, raw_protein_pdb_filepath
    else:
        return None, None, None

def plinder_based_native_complexes_files_iterator(complexes: pd.DataFrame, create_df: bool = True):
    new_dict = defaultdict(list)
    plindex = get_full_plindex()
    # Iterate over all mentioned complexes
    for _, complex in complexes[:3].iterrows():
        pdb_id = complex['pdb_id'].lower()
        ligand_ccd_code = complex['ligand_id']
        pdb_plindex = plindex[plindex['entry_pdb_id']==pdb_id]
        complex_plindex = pdb_plindex[pdb_plindex['ligand_ccd_code']==ligand_ccd_code].reset_index()
        if len(complex_plindex) > 0:
            system_id = complex_plindex['system_id'].values[0]
            plinder_ligand_id = complex_plindex['system_ligand_chains_asym_id'].values[0][0]
            ligand_raw_sdf_file, native_ligand_raw_sdf_file, protein_raw_pdb_file = plinder_based_native_complex_files_generator(system_id=system_id,ligand_id=plinder_ligand_id)
            new_dict['pdb_id'].append(pdb_id)
            new_dict['ligand_id'].append(ligand_ccd_code)
            new_dict['ligand_raw_sdf_filepath'].append(ligand_raw_sdf_file)
            new_dict['native_ligand_raw_sdf_filepath'].append(native_ligand_raw_sdf_file)
            new_dict['protein_raw_pdb_filepath'].append(protein_raw_pdb_file)
            logging.info(f'Raw Files for {pdb_id} and {ligand_ccd_code} created from Plinder {system_id} system.')

    if create_df:
        df = pd.DataFrame.from_dict(new_dict)
        new_df_filepath = Path(os.getenv('DATA_DIR')) / 'raw_complexes.csv'
        df.to_csv(new_df_filepath)
        logging.info(f'Dataframe with Raw Filepaths Created!')
        return df
    else:
        return None



