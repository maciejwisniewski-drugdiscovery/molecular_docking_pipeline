conda create --name dataset_builders -c conda-forge python=3.10.5 rdkit pandas pydantic erdantic -y
conda activate dataset_builders
pip install seaborn click biopython rich datasets pubchempy scikit-learn "huggingface_hub[cli]" --no-cache-dir
conda install conda-forge::openbabel
pip install torch
pip install torch-cluster
pip install git+https://github.com/yusuf1759/prodigy-cryst.git
pip install pinder[all]
pip install git+https://github.com/plinder-org/plinder.git@lig_loader_minor_fix