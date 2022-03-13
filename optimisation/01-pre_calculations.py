"""Pre-calculation to speed up the parameter fitting of Jazzy."""
# optimisation/01-pre_calculations.py
import _pickle as cpickle
import bz2
import os

import config

from jazzy.core import get_charges_from_kallisto_molecule
from jazzy.core import get_covalent_atom_idxs
from jazzy.core import kallisto_molecule_from_rdkit_molecule
from jazzy.core import rdkit_molecule_from_smiles

# Create dictionary with the data from the paper and the SMILES
data_path = os.path.abspath(os.path.join(os.getcwd(), "..", config.DATA_PATH))
original_dataset_filepath = os.path.join(
    data_path, config.GERBER_FREE_ENERGY_DIRNAME, config.ORIGINAL_DATA_FILENAME
)
f = open(original_dataset_filepath, "r")
next(f)
input_data = dict()
for line in f:
    split_results = line.split("\t")
    mol_data = dict()
    mol_data["name"] = split_results[1]
    mol_data["smiles"] = split_results[5].split("\n")[0]
    mol_data["exp_deltag"] = float(split_results[2])
    mol_data["mab_deltag"] = float(split_results[3])
    input_data[int(split_results[0])] = mol_data

# Iterate through dictionary and add pre-calculated data
for idx in input_data:
    mol = input_data[idx]
    rdkit_mol = rdkit_molecule_from_smiles(
        mol["smiles"], minimisation_method=config.MINIMISATION_METHOD
    )
    kallisto_mol = kallisto_molecule_from_rdkit_molecule(rdkit_mol)
    atoms_and_nbrs = get_covalent_atom_idxs(rdkit_mol)
    kallisto_charges = get_charges_from_kallisto_molecule(kallisto_mol, charge=0)
    input_data[idx]["name"] = mol["name"]
    input_data[idx]["smiles"] = mol["smiles"]
    input_data[idx]["rdkit_mol"] = rdkit_mol
    input_data[idx]["kallisto_mol"] = kallisto_mol
    input_data[idx]["atoms_and_nbrs"] = atoms_and_nbrs
    input_data[idx]["charges"] = kallisto_charges

# Drop any failures from the dictionary
for idx in input_data:
    mol = input_data[idx]
    if mol["charges"] is None or len(mol["charges"]) == 0:
        input_data.pop(idx)

# Write data for optimisation out
optuna_path = os.path.join(data_path, config.OPTUNA_DIRNAME)
os.makedirs(optuna_path) if not os.path.exists(optuna_path) else None
pickle_filepath = os.path.join(optuna_path, config.PRECALCULATED_DATA_FILENAME)
with bz2.BZ2File(pickle_filepath, "w") as f:
    cpickle.dump(input_data, f)
