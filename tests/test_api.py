"""Test cases for the API module."""
import numpy as np
import pytest

from jazzy.api import atomic_map_from_smiles
from jazzy.api import atomic_strength_vis_from_smiles
from jazzy.api import atomic_tuples_from_smiles
from jazzy.api import deltag_from_smiles
from jazzy.api import molecular_vector_from_smiles
from jazzy.config import Config
from jazzy.utils import JazzyError

# global jazzy config (parameter)
config = Config()


def test_api_molecular_vector_from_smiles():
    """Correcly calculates molecular vector from API."""
    smiles = "CC"
    # first only hydrogen bond strength
    minimisation_method = "MMFF94"
    only_strengths = True
    vector = molecular_vector_from_smiles(smiles, minimisation_method, only_strengths)
    assert len(vector.values()) == 3
    assert np.isclose(vector["sdc"], 0.5304)
    assert np.isclose(vector["sdx"], 0.0)
    # first only hydrogen bond strength
    minimisation_method = "MMFF94"
    only_strengths = False
    vector = molecular_vector_from_smiles(smiles, minimisation_method, only_strengths)
    assert len(vector.values()) == 6
    assert np.isclose(vector["sdc"], 0.5304)
    assert np.isclose(vector["sdx"], 0.0)
    assert np.isclose(vector["dga"], 4.877683888869463)
    assert np.isclose(vector["dgp"], -0.27805373716258514)
    assert np.isclose(vector["dgtot"], 4.599630151706878)

    smiles = "C(Cl)#C"
    minimisation_method = "MMFF94"
    only_strengths = False
    vector = molecular_vector_from_smiles(smiles, minimisation_method, only_strengths)
    assert len(vector.values()) == 6
    assert np.isclose(vector["sdc"], 0.7748)
    assert np.isclose(vector["sdx"], 0.0)
    assert np.isclose(vector["dga"], -0.5214872778921444)
    assert np.isclose(vector["dgp"], -4.050171713395187)
    assert np.isclose(vector["dgtot"], -4.5716589912873316)


def test_api_molecular_vector_from_smiles_fails_for_invalid_smiles():
    """Correcly fail for incorrect smiles."""
    smiles = "xxx"
    # first only hydrogen bond strength
    minimisation_method = "MMFF94"
    only_strengths = True
    with pytest.raises(JazzyError) as error:
        molecular_vector_from_smiles(smiles, minimisation_method, only_strengths)
    assert error.value.args[0] == "The SMILES 'xxx' appears to be invalid."


def test_deltag_from_smiles():
    """Correcly calculates free energy scalar."""
    smiles = "CC"
    # first only hydrogen bond strength
    minimisation_method = "MMFF94"
    scalar = deltag_from_smiles(smiles, minimisation_method)
    assert np.isclose(scalar, 4.599630151706878)


def test_api_deltag_from_smiles_fails_for_invalid_smiles():
    """Correcly fails for incorrect smiles."""
    smiles = "xxx"
    # first only hydrogen bond strength
    minimisation_method = "MMFF94"
    with pytest.raises(Exception) as error:
        deltag_from_smiles(smiles, minimisation_method)
    assert error.value.args[0] == "The SMILES 'xxx' appears to be invalid."


def test_atomic_tuples_from_smiles():
    """Correcly calculates atomic tuples from atomic map."""
    smiles = "C1CC2=C3C(=CC=C2)C(=CN3C1)"
    # first only hydrogen bond strength
    minimisation_method = "MMFF94"
    tuple_map = atomic_tuples_from_smiles(smiles, minimisation_method)
    assert tuple_map[0][0] == ("z", 6)
    assert tuple_map[0][1] == ("q", 0)
    assert tuple_map[0][2] == ("eeq", -0.1556)
    assert tuple_map[0][3] == ("alp", 7.4365)
    assert tuple_map[0][4] == ("hyb", "sp3")
    assert tuple_map[0][5] == ("num_lp", 0)
    assert tuple_map[0][6] == ("sdc", 0.0)
    assert tuple_map[0][7] == ("sdx", 0.0)
    assert tuple_map[0][8] == ("sa", 0.0)
    assert tuple_map[3][0] == ("z", 6)
    assert tuple_map[3][1] == ("q", 0)
    assert tuple_map[3][2] == ("eeq", 0.1174)
    assert tuple_map[3][3] == ("alp", 8.4498)
    assert tuple_map[3][4] == ("hyb", "sp2")
    assert tuple_map[3][5] == ("num_lp", 0)
    assert tuple_map[3][6] == ("sdc", 0.0)
    assert tuple_map[3][7] == ("sdx", 0.0)
    assert tuple_map[3][8] == ("sa", 0.0)


def test_atomic_tuples_from_smiles_fails_for_invalid_smiles():
    """It fails for en invalid SMILES."""
    smiles = "xxx"
    with pytest.raises(Exception) as error:
        atomic_tuples_from_smiles(smiles)
    assert error.value.args[0] == "The SMILES 'xxx' appears to be invalid."


def test_atomic_map_from_smiles():
    """It correctly condenses an atomic map."""
    smiles = "C1CC2=C3C(=CC=C2)C(=CN3C1)"
    # first only hydrogen bond strength
    minimisation_method = "MMFF94"
    condensed_map = atomic_map_from_smiles(smiles, minimisation_method)
    assert np.isclose(condensed_map[0]["alp"], 7.4365)
    assert np.isclose(condensed_map[0]["eeq"], -0.1556)
    assert condensed_map[0]["hyb"] == "sp3"


def test_atomic_map_from_smiles_fails_for_invalid_smiles():
    """It fails for an invalid SMILES."""
    smiles = "xxx"
    with pytest.raises(Exception) as error:
        atomic_map_from_smiles(smiles)
    assert error.value.args[0] == "The SMILES 'xxx' appears to be invalid."


def test_atomic_strength_vis_from_smiles():
    """It creates an SVG image from SMILES."""
    smiles = "C1CC2=C3C(=CC=C2)C(=CN3C1)"
    minimisation_method = "MMFF94"
    img_txt = atomic_strength_vis_from_smiles(
        smiles=smiles,
        minimisation_method=minimisation_method,
        encode=True,
        fig_size=(100, 100),
        flatten_molecule=False,
        highlight_atoms=False,
        ignore_sa=False,
        ignore_sdc=False,
        ignore_sdx=False,
        sdc_treshold=0.0,
        sdx_treshold=0.0,
        sa_treshold=0.0,
    )
    # extract small amount from encoded image
    got = img_txt[:10]
    # first 10 characters in base64 encoded image
    want = b"PD94bWwgdm"
    assert got == want
