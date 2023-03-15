"""Helpers for the parameter fitting of Jazzy."""
# optimisation/helpers.py
import _pickle as cpickle
import bz2
import json
import os

import config
import optuna
from optuna.samplers import TPESampler


def load_data_configuration(study_filename, params_filename):
    """Abstraction for input and output loading and configuration setting.

    The method only accepts the name of the output study and best parameters
    files, whilst the rest of the parameters must be defined in a config file.

    """
    # load the input
    data_path = os.path.abspath(os.path.join(os.getcwd(), "..", config.DATA_PATH))
    optuna_path = os.path.join(data_path, config.OPTUNA_DIRNAME)
    data_filepath = os.path.join(optuna_path, config.PRECALCULATED_DATA_FILENAME)
    data = bz2.BZ2File(data_filepath, "rb")
    input_data = cpickle.load(data)

    # configure output
    study_filepath = os.path.join(optuna_path, study_filename)
    params_filepath = os.path.join(optuna_path, params_filename)
    return input_data, study_filepath, params_filepath


def run_optimisation(objective, study_filepath, verbose=False):
    """Abstraction for Optuna fitting.

    Includes early stopping logic and fixed seed for reproducibility.
    The fitting can be either verbose or not but always dumps a
    pickle file with the full logs of the process.

    """
    # early stopping logic (https://github.com/optuna/optuna/issues/1001)
    class EarlyStoppingExceeded(optuna.exceptions.OptunaError):
        early_stop = config.OPTUNA_EARLY_STOPPING
        early_stop_count = 0
        best_score = None

    def early_stopping_opt(study, trial):
        if EarlyStoppingExceeded.best_score is None:
            EarlyStoppingExceeded.best_score = study.best_value

        if study.best_value < EarlyStoppingExceeded.best_score:
            EarlyStoppingExceeded.best_score = study.best_value
            EarlyStoppingExceeded.early_stop_count = 0

        else:
            if (
                EarlyStoppingExceeded.early_stop_count
                > EarlyStoppingExceeded.early_stop
            ):
                EarlyStoppingExceeded.early_stop_count = 0
                raise EarlyStoppingExceeded()
            else:
                EarlyStoppingExceeded.early_stop_count = (
                    EarlyStoppingExceeded.early_stop_count + 1
                )
        return

    # run the optimisation with sampler for reproducibility
    if not verbose:
        optuna.logging.set_verbosity(optuna.logging.WARNING)
    sampler = TPESampler(seed=5)
    study = optuna.create_study(sampler=sampler)

    try:
        study.optimize(objective, timeout=None, callbacks=[early_stopping_opt])
    except EarlyStoppingExceeded:
        print(
            f"Early stopping exceeded: No new best scores \
            after {config.OPTUNA_EARLY_STOPPING} iterations"
        )
    print(f"The best parameters were {study.best_params}")

    # write the results out
    with bz2.BZ2File(study_filepath, "w") as f:
        cpickle.dump(study, f)


def dump_parameters_to_json(study_filepath, params_filepath):
    """Helper for dumping the study best parameters.

    Accepts the path to a study file, unpickles it,
    and dumps the best parameters into a serialisable.

    """
    # read the parameters
    with bz2.BZ2File(study_filepath, "rb") as f:
        study = cpickle.load(f)
    best_params = study.best_trial.params

    # dump the parameters
    with open(params_filepath, "w") as f:
        json.dump(best_params, f, indent=4)
