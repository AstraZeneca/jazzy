"""Test cases for the visualisations methods."""
import base64
import pytest

import numpy as np
from rdkit import Chem

from jazzy.core import calculate_polar_strength_map
from jazzy.core import get_charges_from_kallisto_molecule
from jazzy.core import get_covalent_atom_idxs
from jazzy.core import kallisto_molecule_from_rdkit_molecule
from jazzy.core import rdkit_molecule_from_smiles
from jazzy.visualisation import _create_color_scale
from jazzy.visualisation import depict_strengths
from jazzy.visualisation import _exclude_hydrogens
from jazzy.visualisation import _get_highlighted_atoms_and_strength_colors
from jazzy.visualisation import _increase_explicit_hydrogens
from jazzy.visualisation import _increase_explicit_hydrogen_for_bond_atom
from jazzy.visualisation import _remove_excluded_hydrogens
from jazzy.visualisation import _remove_strong_acceptor_hydrogens
from jazzy.visualisation import _set_acceptor_props
from jazzy.visualisation import _set_donor_props
from jazzy.visualisation import _zero_positive_value_check


def test_depict_strengths():
    """It correctly depicts strengths."""
    rdkit_molecule = rdkit_molecule_from_smiles("OO")
    kallisto_molecule = kallisto_molecule_from_rdkit_molecule(rdkit_molecule)
    atoms_and_nbrs = get_covalent_atom_idxs(rdkit_molecule)
    charges = get_charges_from_kallisto_molecule(kallisto_molecule, 0)
    atomic_map = calculate_polar_strength_map(
        rdkit_molecule, kallisto_molecule, atoms_and_nbrs, charges
    )
    img_txt = depict_strengths(
        rdkit_molecule,
        atomic_map,
        flatten_molecule=True,
    )

    base64_hash = base64.b64encode(img_txt.encode("utf-8"))
    assert str(base64_hash[:20]) == "b'PD94bWwgdmVyc2lvbj0n'"


def test_get_highlighted_atoms_and_strength_colors():
    """It correctly gets highlighted atoms and strength colors."""
    rdkit_molecule = rdkit_molecule_from_smiles("OO")
    kallisto_molecule = kallisto_molecule_from_rdkit_molecule(rdkit_molecule)
    atoms_and_nbrs = get_covalent_atom_idxs(rdkit_molecule)
    charges = get_charges_from_kallisto_molecule(kallisto_molecule, 0)
    atomic_map = calculate_polar_strength_map(
        rdkit_molecule, kallisto_molecule, atoms_and_nbrs, charges
    )
    mw = Chem.RWMol(rdkit_molecule)

    sa_threshold = 0.5
    sdc_threshold = 0.5
    sdx_threshold = 0.5
    ignore_sa = True
    ignore_sdc = False
    ignore_sdx = False

    exclude_hydrogens = _exclude_hydrogens(
        mw,
        atomic_map,
        sa_threshold,
        sdc_threshold,
        sdx_threshold,
        ignore_sa,
        ignore_sdc,
        ignore_sdx,
    )

    mol = mw.GetMol()
    colors = _get_highlighted_atoms_and_strength_colors(mol, True)
    assert colors[0][0] == 2
    assert colors[0][1] == 3

    sa_threshold = 0.5
    sdc_threshold = 0.5
    sdx_threshold = 0.5
    ignore_sa = False
    ignore_sdc = False
    ignore_sdx = False

    exclude_hydrogens = _exclude_hydrogens(
        mw,
        atomic_map,
        sa_threshold,
        sdc_threshold,
        sdx_threshold,
        ignore_sa,
        ignore_sdc,
        ignore_sdx,
    )

    mol = mw.GetMol()
    colors = _get_highlighted_atoms_and_strength_colors(mol, True)
    assert colors[0][0] == 0
    assert colors[0][1] == 1
    assert colors[0][2] == 2
    assert colors[0][3] == 3


def test_exclude_hydrogens():
    """It successfully creates an `excluded_hydrogen` list from RDKit molecule."""
    rdkit_molecule = rdkit_molecule_from_smiles("OO")
    kallisto_molecule = kallisto_molecule_from_rdkit_molecule(rdkit_molecule)
    atoms_and_nbrs = get_covalent_atom_idxs(rdkit_molecule)
    charges = get_charges_from_kallisto_molecule(kallisto_molecule, 0)
    atomic_map = calculate_polar_strength_map(
        rdkit_molecule, kallisto_molecule, atoms_and_nbrs, charges
    )
    mw = Chem.RWMol(rdkit_molecule)

    # acceptor strengths
    sa_threshold = 2.0
    sdc_threshold = 2.0
    sdx_threshold = 2.0
    ignore_sa = False
    ignore_sdc = False
    ignore_sdx = False
    exclude_hydrogens = _exclude_hydrogens(
        mw,
        atomic_map,
        sa_threshold,
        sdc_threshold,
        sdx_threshold,
        ignore_sa,
        ignore_sdc,
        ignore_sdx,
    )
    assert len(exclude_hydrogens) == 2
    assert exclude_hydrogens[0] == 2
    assert exclude_hydrogens[1] == 3

    # donor strengths
    sa_threshold = 0.5
    sdc_threshold = 0.5
    sdx_threshold = 0.5
    ignore_sa = True
    ignore_sdc = False
    ignore_sdx = False
    exclude_hydrogens = _exclude_hydrogens(
        mw,
        atomic_map,
        sa_threshold,
        sdc_threshold,
        sdx_threshold,
        ignore_sa,
        ignore_sdc,
        ignore_sdx,
    )

    # extract hydrogen donor strength
    hydrogen_donor_strength = float(mw.GetAtomWithIdx(2).GetProp("sd"))
    assert len(exclude_hydrogens) == 0
    assert np.isclose(hydrogen_donor_strength, 1.0593)


def test_set_donor_props():
    """It correctly sets donor properties."""
    # RDKit carbon molecule
    m = rdkit_molecule_from_smiles("C")
    atom = m.GetAtomWithIdx(0)
    sd = 10.0
    sd_threshold = 7.0
    ignore_sd = False
    condition = _set_donor_props(atom, sd, sd_threshold, ignore_sd)
    assert atom.GetProp("atomNote") == str(sd)
    assert atom.GetProp("sd") == str(sd)
    assert condition == True

    atom2 = m.GetAtomWithIdx(0)
    sd = 0.0
    sd_threshold = 7.0
    ignore_sd = False
    condition = _set_donor_props(atom, sd, sd_threshold, ignore_sd)
    assert condition == False

    # RDKit oxygen molecule
    m = rdkit_molecule_from_smiles("O")
    atom = m.GetAtomWithIdx(0)
    sd = 10.0
    sd_threshold = 7.0
    ignore_sd = False
    condition = _set_donor_props(atom, sd, sd_threshold, ignore_sd)
    assert atom.GetProp("atomNote") == str(sd)
    assert atom.GetProp("sd") == str(sd)
    assert condition == True

    atom2 = m.GetAtomWithIdx(0)
    sd = 0.0
    sd_threshold = 7.0
    ignore_sd = False
    condition = _set_donor_props(atom, sd, sd_threshold, ignore_sd)
    assert condition == False


def test_create_color_scale():
    """It creates successfully RGB mappings."""
    # only valid modes (donor, acceptor)
    with pytest.raises(Exception):
        mode = "invalid"
        idx2strength = dict()
        _create_color_scale(idx2strength, mode)

    idx2strength = {1: 10}
    reds = _create_color_scale(idx2strength, mode="donor")
    assert reds[1] == (1.0, 9.1, 9.1)

    idx2strength = {1: 10}
    blues = _create_color_scale(idx2strength, mode="acceptor")
    assert blues[1] == (9.1, 0.7, 1.0)

    idx2strength = {}
    blues = _create_color_scale(idx2strength, mode="acceptor")
    assert len(blues) == 0


def test_remove_strong_acceptor_hydrogens() -> None:
    """It removes hydrogens that are bonded to strong acceptors."""
    # RDKit molecule
    m = rdkit_molecule_from_smiles("CC(=O)C=CC=C")
    mw = Chem.RWMol(m)
    # index 8 is a hydrogen atom within the embedded molecule
    idx = 8
    hs_to_remove = [idx]
    updated_hs_to_remove = _remove_strong_acceptor_hydrogens(mw, hs_to_remove)
    assert hs_to_remove == updated_hs_to_remove
    # extract neighbor and set as stronh acceptor
    hs_to_remove = [idx]
    neighbor_atom = mw.GetAtomWithIdx(idx).GetNeighbors()[0]
    # strength is 10.0, threshold is 7.0, we do not ignore acceptors (False)
    _set_acceptor_props(neighbor_atom, 10.0, 7.0, False)
    updated_hs_to_remove = _remove_strong_acceptor_hydrogens(mw, hs_to_remove)
    assert updated_hs_to_remove == []


def test_increase_explicit_hydrogens_for_bond_atom() -> None:
    """It increases the number of explicit hydrogens for atom with index `idx`."""
    # RDKit molecule
    m = Chem.MolFromSmiles("CC(=O)C=CC=C")
    # bond start atom index
    bidx = 0
    # bond end atom index
    eidx = 1

    # remove_bidx: True && remove_edix: True
    mw = Chem.RWMol(m)
    ai_to_remove = list()  # type: ignore
    remove_bidx = True
    remove_eidx = True
    batom = mw.GetAtomWithIdx(bidx)
    eatom = mw.GetAtomWithIdx(eidx)
    assert batom.GetNumExplicitHs() == 0
    assert eatom.GetNumExplicitHs() == 0
    # increase explicit number of hydrogens
    mw, ai_to_remove = _increase_explicit_hydrogen_for_bond_atom(
        mw, remove_bidx, bidx, remove_eidx, eidx, ai_to_remove
    )
    # both indices should be in ai_to_remove list
    assert len(ai_to_remove) == 2
    assert ai_to_remove[0] == bidx
    assert ai_to_remove[1] == eidx
    bwant = 1
    ewant = 1
    batom = mw.GetAtomWithIdx(bidx)
    eatom = mw.GetAtomWithIdx(eidx)
    got = batom.GetNumExplicitHs()
    assert got == bwant
    got = eatom.GetNumExplicitHs()
    assert got == ewant

    # remove_bidx: False && remove_eidx: True
    mw = Chem.RWMol(m)
    ai_to_remove = list()
    remove_bidx = False
    remove_eidx = True
    batom = mw.GetAtomWithIdx(bidx)
    eatom = mw.GetAtomWithIdx(eidx)
    assert batom.GetNumExplicitHs() == 0
    assert eatom.GetNumExplicitHs() == 0
    # increase explicit number of hydrogens
    mw, ai_to_remove = _increase_explicit_hydrogen_for_bond_atom(
        mw, remove_bidx, bidx, remove_eidx, eidx, ai_to_remove
    )
    # only eidx index should be in ai_to_remove list
    assert len(ai_to_remove) == 1
    assert ai_to_remove[0] == eidx
    bwant = 1
    ewant = 0
    batom = mw.GetAtomWithIdx(bidx)
    eatom = mw.GetAtomWithIdx(eidx)
    got = batom.GetNumExplicitHs()
    assert got == bwant
    got = eatom.GetNumExplicitHs()
    assert got == ewant

    # remove_bidx: True && remove_eidx: False
    mw = Chem.RWMol(m)
    ai_to_remove = list()
    remove_bidx = True
    remove_eidx = False
    batom = mw.GetAtomWithIdx(bidx)
    eatom = mw.GetAtomWithIdx(eidx)
    assert batom.GetNumExplicitHs() == 0
    assert eatom.GetNumExplicitHs() == 0
    # increase explicit number of hydrogens
    mw, ai_to_remove = _increase_explicit_hydrogen_for_bond_atom(
        mw, remove_bidx, bidx, remove_eidx, eidx, ai_to_remove
    )
    # only bidx index should be in ai_to_remove list
    assert len(ai_to_remove) == 1
    assert ai_to_remove[0] == bidx
    bwant = 0
    ewant = 1
    batom = mw.GetAtomWithIdx(bidx)
    eatom = mw.GetAtomWithIdx(eidx)
    got = batom.GetNumExplicitHs()
    assert got == bwant
    got = eatom.GetNumExplicitHs()
    assert got == ewant

    # remove_bidx: False && remove_eidx: False
    mw = Chem.RWMol(m)
    ai_to_remove = list()
    remove_bidx = False
    remove_eidx = False
    batom = mw.GetAtomWithIdx(bidx)
    eatom = mw.GetAtomWithIdx(eidx)
    assert batom.GetNumExplicitHs() == 0
    assert eatom.GetNumExplicitHs() == 0
    # increase explicit number of hydrogens
    mw, ai_to_remove = _increase_explicit_hydrogen_for_bond_atom(
        mw, remove_bidx, bidx, remove_eidx, eidx, ai_to_remove
    )
    # no index should be in ai_to_remove list
    assert len(ai_to_remove) == 0
    bwant = 0
    ewant = 0
    batom = mw.GetAtomWithIdx(bidx)
    eatom = mw.GetAtomWithIdx(eidx)
    got = batom.GetNumExplicitHs()
    assert got == bwant
    got = eatom.GetNumExplicitHs()
    assert got == ewant


def test_remove_excluded_hydrogens() -> None:
    """It increases the number of explicit hydrogens for list of implicit hydrogens."""
    m = Chem.MolFromSmiles("CC(=O)C=CC=CN")
    mw = Chem.RWMol(m)
    idx = 0
    exclude = [idx]
    _remove_excluded_hydrogens(mw, exclude)
    atom = mw.GetAtomWithIdx(idx)
    want = 1
    got = atom.GetNumExplicitHs()
    assert got == want


def test_increase_explicit_hydrogens() -> None:
    """It increases the number of explicit hydrogens for atom with index `idx`."""
    # increase for non-hydrogen atoms
    m = Chem.MolFromSmiles("CC(=O)C=CC=C")
    mw = Chem.RWMol(m)
    idx = 0
    atom = mw.GetAtomWithIdx(idx)
    _increase_explicit_hydrogens(mw, idx)
    want = 1
    got = atom.GetNumExplicitHs()
    assert got == want
    # do not increase for hydrogen
    m = Chem.MolFromSmiles("[H]")
    mw = Chem.RWMol(m)
    idx = 0
    atom = mw.GetAtomWithIdx(idx)
    _increase_explicit_hydrogens(mw, idx)
    want = 0
    got = atom.GetNumExplicitHs()
    assert got == want


def test_function_fails_for_negative_input() -> None:
    """It exits with a ValueError when a negative value is entered."""
    with pytest.raises(Exception):
        value = -1
        _zero_positive_value_check(value)
