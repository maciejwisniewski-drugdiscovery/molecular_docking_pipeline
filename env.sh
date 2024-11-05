export DATA_DIR=data

export LIGANDS_DIR=${DATA_DIR}/ligands
export PROTEINS_DIR=${DATA_DIR}/proteins
export NATIVE_LIGANDS_DIR=${DATA_DIR}/native_ligands
export DOCKING_RESULTS_DIR=${DATA_DIR}/docking_results


export DATASET_CACHE_DIR=dataset_cache
export LOG_DIR=log


mkdir -p $DATA_DIR $DATASET_CACHE_DIR $LOG_DIR
mkdir -p $LIGANDS_DIR $PROTEINS_DIR $NATIVE_LIGANDS_DIR $DOCKING_RESULTS_DIR

mkdir -p $LIGANDS_DIR/raw
mkdir -p $PROTEINS_DIR/raw
mkdir -p $NATIVE_LIGANDS_DIR/raw

mkdir -p $LIGANDS_DIR/processed_pdbqt
mkdir -p $PROTEINS_DIR/processed_pdbqt
mkdir -p $NATIVE_LIGANDS_DIR/processed_pdbqt

mkdir -p $LIGANDS_DIR/processed_sdf
mkdir -p $NATIVE_LIGANDS_DIR/processed_sdf
mkdir -p $PROTEINS_DIR/processed_pdb


export $SMINA=smina
export PINDER_BASE_DIR='/ephemeral/users/maciej.wisniewski/data'
export PLINDER_MOUNT='/ephemeral/users/maciej.wisniewski/data'
export PLINDER_RELEASE='2024-06'
export PLINDER_ITERATION='v2'
export PLINDER_OFFLINE='true'