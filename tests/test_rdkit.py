"""Test cases for the rdkit methods."""
import pytest
from rdkit.Chem.rdchem import Mol

from jazzy.core import get_all_neighbours
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


def test_invalid_minimisation_method() -> None:
    """It exits with a ValueError when an invalid method is entered."""
    with pytest.raises(ValueError):
        smiles = "CC"
        minimisation_method = "MMFF95"
        rdkit_molecule_from_smiles(
            smiles=smiles,
            minimisation_method=minimisation_method,
        )


def test_invalid_embedding_type() -> None:
    """It exits with a ValueError when an invalid method is entered."""
    with pytest.raises(ValueError):
        smiles = "CC"
        rdkit_molecule_from_smiles(
            smiles=smiles,
            embedding_type="4d",
        )


def test_warning_2d_embedding_with_minimisation() -> None:
    """It ignores the clash between 2D embedding and energy minimisation."""
    smiles = "Brc1ccccc1OCCCOc1cccc2cccnc12"
    m = rdkit_molecule_from_smiles(
        smiles=smiles, embedding_type="2D", minimisation_method="MMFF94"
    )
    assert isinstance(m, Mol)


def test_embedding_fails_with_fewer_iterations() -> None:
    """It returns a None when a particular molecule cannot be embedded."""
    smiles = "Bc1ccc(cc1)S(=O)(=O)c2ccc3SC(CC=C)[C@@H](C(=O)N4[C@@H]5C[C@@H]4[C@H]6C[C@@H]56)c3c2"  # noqa
    m = rdkit_molecule_from_smiles(smiles=smiles, embedding_max_iterations=1)
    assert m is None


def test_embedding_succeeds_with_more_iterations() -> None:
    """It returns a None when a particular molecule cannot be embedded."""
    smiles = "Bc1ccc(cc1)S(=O)(=O)c2ccc3SC(CC=C)[C@@H](C(=O)N4[C@@H]5C[C@@H]4[C@H]6C[C@@H]56)c3c2"  # noqa
    m = rdkit_molecule_from_smiles(smiles=smiles, embedding_max_iterations=100)
    assert m is not None


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
    alpha, beta, gamma = get_all_neighbours(rdkit_molecule, atoms_and_idxs)
    assert sorted(alpha[8]) == sorted(want[0])
    assert sorted(beta[8]) == sorted(want[1])
    assert sorted(gamma[8]) == sorted(want[2])
