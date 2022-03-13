"""Parameter fitting of Jazzy using Optuna and sklearn."""
# optimisation/02-deltag_fitting.py
import config
from helpers import load_data_configuration
from helpers import run_optimisation
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error

from jazzy.core import calculate_delta_apolar
from jazzy.core import calculate_delta_interaction
from jazzy.core import calculate_delta_polar
from jazzy.core import calculate_polar_strength_map

"""
Notes:
The parameter optimisation was carried out using optuna==2.3.0
(https://pypi.org/project/optuna/) and scikit-learn==0.24.2
(https://pypi.org/project/scikit-learn/).

Best trial:
[I 2021-11-30 08:54:47,298] Trial 859 finished with value: 3.817880136986301
and parameters: {'gd': -0.9076664181509553, 'ga': -16.13143317348325,
'g0': 1.8842791231644282, 'gs': 0.046698073586512054, 'gr': -3.6429426011852524,
'gpi1': -1.6022604035599135, 'gpi2': -1.1740853831959548, 'gi': 4.999560447432514,
'F': 0.5143359987215719, 'expd': 0.5041243951898653, 'expa': 0.3436880638589098}.
Best is trial 859 with value: 3.817880136986301.
"""


def deltag_objective(trial):
    """Minimisation function for Jazzy.

    Calculates of the free energy of hydration on a set of molecules
    and computes the loss between predicted and annotated values.

    Args:
    trial: optuna.trial.Trial object

    Returns:
    Loss between predicted and annotated values (float)

    """
    # select parameter ranges                     # parameters from the paper
    gd = trial.suggest_uniform("gd", -30.0, 0.0)  # gd=-139.0
    ga = trial.suggest_uniform("ga", -20.0, 0.0)  # ga=-32.0
    g0 = trial.suggest_uniform("g0", 0.0, 3.0)  # g0=5.53
    gs = trial.suggest_uniform("gs", 0.0, 1.0)  # gs=0.031
    gr = trial.suggest_uniform("gr", -6.0, 0.0)  # gr=-4.39
    gpi1 = trial.suggest_uniform("gpi1", -10.0, 0.0)  # gpi1=-1.82
    gpi2 = trial.suggest_uniform("gpi2", -10.0, 0.0)  # gpi2=-1.29
    gi = trial.suggest_uniform("gi", 0.0, 5.0)  # gi=37.2
    f = trial.suggest_uniform("F", 0.1, 0.7)  # F=0.325
    expd = trial.suggest_uniform("expd", 0.5, 0.7)  # expd=0.68
    expa = trial.suggest_uniform("expa", 0.3, 0.5)  # expa=0.51

    # generate y_true and y_pred
    y_exp = []
    y_calc = []
    for idx in input_data:
        mol = input_data[idx]
        rdkit_mol = mol["rdkit_mol"]
        kallisto_mol = mol["kallisto_mol"]
        atoms_and_nbrs = mol["atoms_and_nbrs"]
        charges = mol["charges"]
        y_exp.append(mol["exp_deltag"])

        # calculate polar strength map
        atomic_map = calculate_polar_strength_map(
            rdkit_mol, kallisto_mol, atoms_and_nbrs, charges
        )

        # calculate individual terms and finally append their sum
        dgp = round(
            calculate_delta_polar(
                atomic_map, atoms_and_nbrs, gd=gd, ga=ga, expd=expd, expa=expa
            ),
            4,
        )

        dga = round(
            calculate_delta_apolar(
                rdkit_mol, atomic_map, g0=g0, gs=gs, gr=gr, gpi1=gpi1, gpi2=gpi2
            ),
            4,
        )

        dgi = round(
            calculate_delta_interaction(
                rdkit_mol, atomic_map, atoms_and_nbrs, gi=gi, expa=expa, f=f
            ),
            4,
        )

        calc_deltag = dgp + dga + dgi
        y_calc.append(calc_deltag)

    # metric to minimise
    if config.LOSS_FUNCTION == "mean_absolute_error":
        loss = mean_absolute_error(y_exp, y_calc)
    elif config.LOSS_FUNCTION == "mean_squared_error":
        loss = mean_squared_error(y_exp, y_calc)
    elif config.LOSS_FUNCTION == "root_mean_square_error":
        loss = mean_squared_error(y_exp, y_calc, squared=False)
    else:
        raise ValueError("Please define a valid loss function")
    return loss


# load input/ouput configuration
input_data, study_filepath = load_data_configuration(
    config.OPTUNA_DELTAG_STUDY_FILENAME
)

# run optimisation
study = run_optimisation(deltag_objective, study_filepath, verbose=config.VERBOSE)
