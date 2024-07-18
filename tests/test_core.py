"""Test cases for the core methods."""
import numpy as np

from jazzy.core import get_atom_and_nbrs_charges
from jazzy.core import get_atom_neighbours
from jazzy.core import get_charges_from_atom_list
from jazzy.core import get_charges_from_kallisto_molecule
from jazzy.core import get_covalent_atom_idxs
from jazzy.core import get_donor_atom_strength
from jazzy.core import kallisto_molecule_from_rdkit_molecule
from jazzy.core import rdkit_molecule_from_smiles


def test_calculate_correct_neighbours_for_atom():
    """Calculates alpha, beta, and gamma neighbours for specific atom."""
    want = [[4, 9, 19], [3, 5, 10, 20], [2, 3, 6, 10, 11, 16]]
    smiles = "C1CC2=C3C(=CC=C2)C(=CN3C1)"
    rdkit_molecule = rdkit_molecule_from_smiles(smiles=smiles)
    atoms_and_idxs = get_covalent_atom_idxs(rdkit_molecule=rdkit_molecule)
    alpha, beta, gamma = get_atom_neighbours(8, atoms_and_idxs)
    assert sorted(alpha) == sorted(want[0])
    assert sorted(beta) == sorted(want[1])
    assert sorted(gamma) == sorted(want[2])


def test_get_charges_from_atom_list():
    """It returns the charges for different atoms in a molecule."""
    want = [0.5, 0.02]
    charges = [0, 1, 0.5, 3, 0.02, 5, 6]
    atom = [2, 4]
    got = get_charges_from_atom_list(atom, charges)
    assert got[0] == want[0]
    assert got[1] == want[1]


def test_donor_atom_strength():
    """It calculates the correct donor strength."""
    want = [-8.749496986634238]
    smiles = "C1CC2=C3C(=CC=C2)C(=CN3C1)"
    rdkit_molecule = rdkit_molecule_from_smiles(smiles, minimisation_method="MMFF94")
    atoms_and_nbrs = get_covalent_atom_idxs(rdkit_molecule)
    kallisto_molecule = kallisto_molecule_from_rdkit_molecule(
        rdkit_molecule=rdkit_molecule
    )
    charge = 0
    eeq = get_charges_from_kallisto_molecule(kallisto_molecule, charge)
    strength = get_donor_atom_strength(1, atoms_and_nbrs, eeq)
    assert np.isclose(strength, want[0])


def test_get_atom_and_nbrs_charges():
    """It calculates the correct charge and delta."""
    want = [-0.14491844176443144, 0.007563609354876652]
    smiles = "C1CC2=C3C(=CC=C2)C(=CN3C1)"
    rdkit_molecule = rdkit_molecule_from_smiles(smiles, minimisation_method="MMFF94")
    atoms_and_idxs = get_covalent_atom_idxs(rdkit_molecule)
    kallisto_molecule = kallisto_molecule_from_rdkit_molecule(
        rdkit_molecule=rdkit_molecule
    )
    charge = 0
    eeq = get_charges_from_kallisto_molecule(kallisto_molecule, charge)
    q_at, q_delta = get_atom_and_nbrs_charges(1, atoms_and_idxs, eeq)
    assert np.isclose(q_at, want[0], 5)
    assert np.isclose(q_delta, want[1], 5)
