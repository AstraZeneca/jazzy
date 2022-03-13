"""Configuration for the parameter fitting of Jazzy."""
# optimisation/config.py
# directory configuration
DATA_PATH = "data"
GERBER_FREE_ENERGY_DIRNAME = "free_energy_gerber"
OPTUNA_DIRNAME = "optuna_fitting"

# file configuration
ORIGINAL_DATA_FILENAME = "03-free_energies_gerber_curated.txt"
PRECALCULATED_DATA_FILENAME = "01-precalculated_mol_data.pbz2"
OPTUNA_DELTAG_STUDY_FILENAME = "02-optuna_deltag_fitting_study.pbz2"

# logic configuration
MINIMISATION_METHOD = "MMFF94"
VERBOSE = True
OPTUNA_EARLY_STOPPING = 300
LOSS_FUNCTION = "mean_absolute_error"  # "mean_squared_error" # "root_mean_square_error"
