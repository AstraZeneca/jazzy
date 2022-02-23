"""Test cases for the rdkit methods."""
import pytest

from jazzy.core import get_all_neighbors
from jazzy.core import get_covalent_atom_idxs
from jazzy.core import get_lone_pairs
from jazzy.core import rdkit_molecule_from_smiles


def rdkit_molecule_atomic_numbers(rdkit_molecule):
    """Return atomic number list for RDKit molecule."""
    return list(map(lambda x: x.GetAtomicNum(), list(rdkit_molecule.GetAtoms())))


def test_minimisation_works_for_valid_method() -> None:
    """It exists with no error when a valid method is entered."""
    want = [6, 6, 1]
    # define valid minimisation methods
    valid_methods = [None, "MMFF94", "MMFF94s", "UFF"]

    smiles = "CC"
    # None case
    m = rdkit_molecule_from_smiles(
        smiles=smiles,
        minimisation_method=valid_methods[0],
    )
    atomic_nums = rdkit_molecule_atomic_numbers(m)
    assert atomic_nums[0] == want[0]
    assert atomic_nums[1] == want[1]
    assert atomic_nums[2] == want[2]

    # MMFF94 case
    m = rdkit_molecule_from_smiles(
        smiles=smiles,
        minimisation_method=valid_methods[1],
    )
    atomic_nums = rdkit_molecule_atomic_numbers(m)
    assert atomic_nums[0] == want[0]
    assert atomic_nums[1] == want[1]
    assert atomic_nums[2] == want[2]

    # MMFF94s case
    m = rdkit_molecule_from_smiles(
        smiles=smiles,
        minimisation_method=valid_methods[2],
    )
    atomic_nums = rdkit_molecule_atomic_numbers(m)
    assert atomic_nums[0] == want[0]
    assert atomic_nums[1] == want[1]
    assert atomic_nums[2] == want[2]

    # UFF case
    m = rdkit_molecule_from_smiles(
        smiles=smiles,
        minimisation_method=valid_methods[3],
    )
    atomic_nums = rdkit_molecule_atomic_numbers(m)
    assert atomic_nums[0] == want[0]
    assert atomic_nums[1] == want[1]
    assert atomic_nums[2] == want[2]


def test_minimisation_fails_for_nonvalid_method() -> None:
    """It exists with a ValueError when a nonvalid method is entered."""
    with pytest.raises(Exception):
        smiles = "CC"
        minimisation_method = "MMFF95"
        m = rdkit_molecule_from_smiles(
            smiles=smiles,
            minimisation_method=minimisation_method,
        )
        assert m is None


def test_embedding_fails_for_particular_smiles() -> None:
    """It returns a None when a particular molecule cannot be embedded."""
    smiles = "COc1cc2cc(OC(=C2C=C)N)c1OC"
    m = rdkit_molecule_from_smiles(smiles=smiles)
    assert m is None


def test_calculate_correct_number_of_lone_pairs():
    """Calculate the correct number of lone pairs for different molecules."""
    want = {
        "O": [{"idx": 0, "lps": 2}],
        "[OH3+]": [{"idx": 0, "lps": 1}],
        "[NH3+]": [{"idx": 0, "lps": 0}],
        "O=S=O": [{"idx": 0, "lps": 2}],
        "F[Xe](F)(F)F": [{"idx": 0, "lps": 3}, {"idx": 1, "lps": 2}],
    }
    # unpack dict to smiles
    for _, smiles in enumerate(want):
        cpd = want[smiles]
        rdkit_molecule = rdkit_molecule_from_smiles(
            smiles=smiles, minimisation_method=None
        )
        # unpack atomic details
        for atom in cpd:
            atom_idx = atom["idx"]
            atom_lps = atom["lps"]
            calc_lps = get_lone_pairs(rdkit_molecule.GetAtomWithIdx(atom_idx))
            assert atom_lps == calc_lps


def test_calculate_correct_covalent_idx_list():
    """It calculates the correct covalent bounding list for a molecule."""
    want = [[1, 11, 12, 13], [0, 2, 14, 15]]
    smiles = "C1CC2=C3C(=CC=C2)C(=CN3C1)"
    rdkit_molecule = rdkit_molecule_from_smiles(smiles)
    atoms_and_idxs = get_covalent_atom_idxs(rdkit_molecule=rdkit_molecule)
    assert sorted(atoms_and_idxs[0]) == sorted(want[0])
    assert sorted(atoms_and_idxs[1]) == sorted(want[1])


def test_calculate_correct_all_neighbours_for_atom():
    """Calculates alpha, beta, and gamma neighbours for specific atom."""
    want = [[4, 9, 19], [3, 5, 10, 20], [2, 3, 6, 10, 11, 16]]
    smiles = "C1CC2=C3C(=CC=C2)C(=CN3C1)"
    rdkit_molecule = rdkit_molecule_from_smiles(smiles=smiles)
    atoms_and_idxs = get_covalent_atom_idxs(rdkit_molecule=rdkit_molecule)
    alpha, beta, gamma = get_all_neighbors(rdkit_molecule, atoms_and_idxs)
    assert sorted(alpha[8]) == sorted(want[0])
    assert sorted(beta[8]) == sorted(want[1])
    assert sorted(gamma[8]) == sorted(want[2])
