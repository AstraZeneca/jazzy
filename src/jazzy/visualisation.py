"""RDKit functions to convert Jazzy data into rendering."""
# src/jazzy/visualisation.py
from typing import Tuple

import rdkit
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from jazzy.config import ANNOTATION_FONT_SCALE
from jazzy.config import ROUNDING_DIGITS


def _zero_positive_value_check(v: float):
    """Assert that the input is zero or positive.

    Raises:
    ValueError: An error if the input is negative

    """
    if v < 0:
        raise ValueError("The value should be zero or positive (got {})".format(v))


def _increase_explicit_hydrogens(rwmol: Chem.rdchem.RWMol, idx: int):
    """Increase number of explicit hydrogens for atom with index `idx`.

    Args:
    rwmol: An RDKit RWmolecule (rdkit.Chem.rdchem.RWMol)
    idx: index of atom (int)

    """
    oa = rwmol.GetAtomWithIdx(idx)
    if oa.GetAtomicNum() > 1:
        oa.SetNumExplicitHs(oa.GetNumExplicitHs() + 1)


def _increase_explicit_hydrogen_for_bond_atom(
    rwmol: Chem.rdchem.RWMol,
    remove_bidx: bool,
    bidx: int,
    remove_eidx: bool,
    eidx: int,
    ai_to_remove: list,
) -> Tuple[Chem.rdchem.RWMol, list]:
    """Increase number of explicit hydrogens for atom in a bond.

    Args:
    rwmol: An RDKit RWmolecule (rdkit.Chem.rdchem.RWMol)
    remove_bidx: Begin atom in bond will increase explicit hydrogens (bool)
    remove_eidx: End atom in bond will increase explicit hydrogens (bool)

    Returns:
    Tuple with an RDKit RWmolecule and an updated list to remove
    (rdkit.Chem.rdchem.RWMol, list).

    """
    if remove_bidx or remove_eidx:
        if remove_bidx:
            ai_to_remove.append(bidx)
            _increase_explicit_hydrogens(rwmol, eidx)
        if remove_eidx:
            ai_to_remove.append(eidx)
            _increase_explicit_hydrogens(rwmol, bidx)
        rwmol.RemoveBond(bidx, eidx)
    return rwmol, ai_to_remove


def _remove_excluded_hydrogens(rwmol: Chem.rdchem.RWMol, excluded_hydrogens: list):
    """Convert list of hydrogen atoms into explicit hydrogens.

    Convert implicit hydrogen atoms into explicit hydrogen on the atoms where
    they are bonded. E.g., H-N-H is converted to NH2.

    Args:
    rwmol: An RDKit RWmolecule (rdkit.Chem.rdchem.RWMol)
    excluded_hydrogens: list of implicit Hydrogens to convert into explicit ones.

    Returns:
    rdkit_molecule: An RDKit molecule (rdkit.Chem.rdchem.Mol)

    """
    if not excluded_hydrogens:
        return rwmol

    ai_to_remove = list()  # type: ignore
    for bond_idx in reversed(range(rwmol.GetNumBonds())):
        b = rwmol.GetBondWithIdx(bond_idx)
        bidx = b.GetBeginAtomIdx()
        remove_bidx = bidx in excluded_hydrogens
        eidx = b.GetEndAtomIdx()
        remove_eidx = eidx in excluded_hydrogens
        rwmol, ai_to_remove = _increase_explicit_hydrogen_for_bond_atom(
            rwmol,
            remove_bidx,
            bidx,
            remove_eidx,
            eidx,
            ai_to_remove,
        )
    for atom_idx in sorted(ai_to_remove, reverse=True):
        rwmol.RemoveAtom(atom_idx)
    rdkit.Chem.SanitizeMol(rwmol)
    return rwmol.GetMol()


def _remove_strong_acceptor_hydrogens(rwmol: Chem.rdchem.RWMol, hs_to_remove: list):
    """Remove Hydrogens that are bonded to strong acceptors.

    Goes through a list of hydrogens and checks wheather any of their neighbors
    has been annotated with acceptor strength. If so, it removes the hydrogen
    from the removal list to avoid it to collapse with the acceptor atom.

    Args:
    rwmol: An RDKit RWmolecule (rdkit.Chem.rdchem.RWMol)
    hs_to_remove: list of implicit Hydrogens to convert into explicit ones.

    Returns:
    Updated list for hydrogen atoms to be removed.

    """
    updated_hs_to_remove = list()
    for idx in hs_to_remove:
        hs = rwmol.GetAtomWithIdx(idx)
        nbrs = hs.GetNeighbors()
        for nbr in nbrs:
            if not nbr.HasProp("sa"):
                updated_hs_to_remove.append(idx)
                break
    return updated_hs_to_remove


def _create_color_scale(idx2strength: dict, mode: str):
    """Create RGB mapping.

    Normalises values in a dictionary and creates an RGB red or blue scale
    mapping for them. `abs(round(0.9-alpha, ROUNDING_DIGITS))` is a bit hacky
    as it avoids the creation of whites by sacrificing a bit of alpha.
    Also RGB for acceptor has been tweaked a bit to avoid generating scales
    that might hide atoms rendered with blue fonts (e.g., Nitrogen atoms).

    Args:
    idx2strength: mapping between indices and strengths (dict)
    mode: `acceptor` or `donor` (str)

    Returns:
    Mapping between indices and RGB colors (dict).

    """
    valid_modes = ["donor", "acceptor"]
    if mode not in valid_modes:
        raise ValueError("Select valid mode. Either 'donor' or 'acceptor'.")

    strengths = idx2strength.values()
    if len(strengths) == 0:
        return idx2strength
    amin, amax = min(strengths), max(strengths)
    idx2rgb = dict()
    for idx, v in idx2strength.items():
        try:
            alpha = round((v - amin) / (amax - amin), ROUNDING_DIGITS)
        except ZeroDivisionError:
            alpha = round(v, ROUNDING_DIGITS)
        invalpha = abs(round(0.9 - alpha, ROUNDING_DIGITS))
        if mode == "donor":
            rgb = (1.0, invalpha, invalpha)
        else:
            # acceptor
            rgb = (invalpha, 0.7, 1.0)
        idx2rgb[idx] = rgb
    return idx2rgb


def _draw_molecule(
    rwmol: Chem.rdchem.RWMol,
    fig_size: Tuple[int, int],
    atoms_to_highlight: list,
    idx2rgb: dict,
):
    """Draw an RDKit molecule.

    Args:
    rwmol: An RDKit RWmolecule (rdkit.Chem.rdchem.RWMol)

    Returns:
    Text of an SVG image.

    """
    d2d = rdMolDraw2D.MolDraw2DSVG(fig_size[0], fig_size[1])
    d2d.drawOptions().annotationFontScale = ANNOTATION_FONT_SCALE
    d2d.drawOptions().padding = 0.15
    d2d.DrawMolecule(
        rwmol,
        highlightAtoms=atoms_to_highlight,
        highlightAtomColors=idx2rgb,
        highlightBonds=None,
    )
    d2d.FinishDrawing()
    return d2d.GetDrawingText()


def _set_acceptor_props(
    atom: rdkit.Chem.rdchem.Atom, sa: float, sa_threshold: float, ignore_sa: bool
) -> bool:
    """Setting two props to acceptor to avoid mistakes in the highlighting function.

    Args:
    atom: An RDKit atom (rdkit.Chem.rdchem.Atom)
    sa: acceptor strength (float)
    sa_threshold: threshold to show acceptor contribution (float)
    ignore_sa: ignore acceptor contributions (bool)

    Returns:
    If properties have been set to atom (bool).

    """
    condition = False
    if sa != 0 and sa > sa_threshold and not ignore_sa:
        atom.SetProp("atomNote", str(sa))
        atom.SetProp("sa", str(sa))
        condition = True
    return condition


def _set_donor_props(
    atom: rdkit.Chem.rdchem.Atom, sd: float, sd_threshold: float, ignore_sd: bool
) -> bool:
    """Setting two props to donor to avoid mistakes in the highlighting function.

    Args:
    atom: An RDKit atom (rdkit.Chem.rdchem.Atom)
    sd: donor strength (float)
    sd_threshold: threshold to show donor contribution (float)
    ignore_sd: ignore donor contributions (bool)

    Returns:
    If properties have been set to atom (bool).

    """
    condition = False
    if sd != 0 and sd > sd_threshold and not ignore_sd:
        atom.SetProp("atomNote", str(sd))
        atom.SetProp("sd", str(sd))
        condition = True
    return condition


def _exclude_hydrogens(
    rwmol: Chem.rdchem.RWMol,
    atomic_map: dict,
    sa_threshold: float,
    sdc_threshold: float,
    sdx_threshold: float,
    ignore_sa: bool,
    ignore_sdc: bool,
    ignore_sdx: bool,
    rounding_digits=ROUNDING_DIGITS,
):
    """Create `excluded_hydrogens` list from RDKit molecule.

    Args:
    rwmol: An RDKit RWmolecule (rdkit.Chem.rdchem.RWMol)
    atomic_map: molecular map of polar properties (dict)
    sa_threshold: threshold to show acceptor contribution (float)
    sdc_threshold: threshold to show donor contribution on Carbon (float)
    sdx_threshold: threshold to show donor contribution on non-Carbon (float)
    ignore_sa: ignore acceptor contributions (bool)
    ignore_sdc: ignore donor contributions on Carbon (bool)
    ignore_sdx: ignore donor contributions on non-Carbon (bool)

    Returns:
    `excluded_hydrogens` list.

    """
    excluded_hydrogens = list()
    for idx, atom in enumerate(rwmol.GetAtoms()):
        atom_props = atomic_map[idx]

        # not involved in hydrogen bonding
        sa = round(atom_props["sa"], rounding_digits)
        sdc = round(atom_props["sdc"], rounding_digits)
        sdx = round(atom_props["sdx"], rounding_digits)
        if sa == 0.0 and sdc == 0.0 and sdx == 0.0:
            continue

        # acceptor logic
        if _set_acceptor_props(atom, sa, sa_threshold, ignore_sa):
            continue

        # donor logic for carbons
        if _set_donor_props(atom, sdc, sdc_threshold, ignore_sdc):
            continue

        # donor logic for non-carbons
        if _set_donor_props(atom, sdx, sdx_threshold, ignore_sdx):
            continue

        # hydrogen does not fall into above conditions
        # then not strong enough and mark for removal
        if atom_props["z"] == 1:
            excluded_hydrogens.append(idx)

    # remove strong acceptor hydrogens
    if not ignore_sa:
        excluded_hydrogens = _remove_strong_acceptor_hydrogens(
            rwmol, excluded_hydrogens
        )

    return excluded_hydrogens


def _get_highlighted_atoms_and_strength_colors(
    rdkit_molecule: Chem.rdchem.Mol,
    highlight_atoms: bool,
):
    """Get donor and acceptor strengths and create RGB colors from them.

    Args:
    rdkit_molecule: An RDKit RWmolecule (rdkit.Chem.rdchem.RWMol)
    highlight_atoms: show donors in red and acceptors in blue (bool)

    Returns:
    `atoms_to_highlight`, `idx2sa`, `idx2sd`, and `idx2rgb` lists.

    """
    idx2sd = dict()
    idx2sa = dict()
    idx2rgb = dict()
    atoms_to_highlight = list()
    if highlight_atoms:
        atoms = [x for x in rdkit_molecule.GetAtoms()]
        for a in atoms:
            if a.HasProp("sd"):
                idx = a.GetIdx()
                atoms_to_highlight.append(idx)
                idx2sd[idx] = float(a.GetProp("sd"))
            if a.HasProp("sa"):
                idx = a.GetIdx()
                atoms_to_highlight.append(idx)
                idx2sa[idx] = float(a.GetProp("sa"))

        # convert strengths into relative RGB scale and combine the scales
        idx2reds = _create_color_scale(idx2sd, mode="donor")
        idx2blues = _create_color_scale(idx2sa, mode="acceptor")
        idx2rgb = {**idx2reds, **idx2blues}  # type: ignore
    return atoms_to_highlight, idx2rgb


def depict_strengths(
    rdkit_molecule: Chem.rdchem.Mol,
    atomic_map: dict,
    fig_size=(500, 500),
    flatten_molecule=False,
    highlight_atoms=False,
    ignore_sdc=False,
    ignore_sdx=False,
    ignore_sa=False,
    sdc_threshold=0.0,
    sdx_threshold=0.0,
    sa_threshold=0.0,
    rounding_digits=ROUNDING_DIGITS,
):
    """Create an SVG image text from an RDKit molecule and its atomic map.

    The default configuration simply produces a depiction of the input molecule
    and its strengths. `highlight_atoms` highlights atoms in red (donors) and
    blue (acceptors). `flatten_molecule` produces a 2-dimensional depiction of
    the molecule. Any `threshold` parameter allows to set a numeric threshold
    under which strengths are not included in the output depiction.

    Args:
    rdkit_molecule: An RDKit molecule (rdkit.Chem.rdchem.Mol)
    atomic_map: molecular map of polar properties (dict)
    fig_size: size of the depiction in pixels (tuple)
    flatten_molecule: boolean to create 2D depiction (bool)
    highlight_atoms: show donors in red and acceptors in blue (bool)
    ignore_sdc: ignore donor contributions on Carbon (bool)
    ignore_sdx: ignore donor contributions on non-Carbon (bool)
    ignore_sa: ignore acceptor contributions (bool)
    sdc_threshold: threshold to show donor contribution on Carbon (float)
    sdx_threshold: threshold to show donor contribution on non-Carbon (float)
    sa_threshold: threshold to show acceptor contribution (float)

    Returns:
    SVG depiction of given RDKit molecule and its atomic strength map.

    """
    # assert that thresholds are zero or positive
    for v in [sdc_threshold, sdx_threshold, sa_threshold]:
        _zero_positive_value_check(v)

    # copy RDKit molecule otherwise the input would be affected by processing
    rwmol = rdkit.Chem.RWMol(rdkit_molecule)

    # flatten the molecule if required
    # Note: never create SMILES then again RDKit 2D molecule since the mapping
    # could be shuffled! Only recompute 2D coordinates on the 3D molecule.
    if flatten_molecule:
        rdkit.Chem.rdDepictor.Compute2DCoords(rwmol)

    # create list of excluded hydrogens
    excluded_hydrogens = _exclude_hydrogens(
        rwmol,
        atomic_map,
        sa_threshold,
        sdc_threshold,
        sdx_threshold,
        ignore_sa,
        ignore_sdc,
        ignore_sdx,
        rounding_digits,
    )

    # remove excluded hydrogen atoms
    rwmol = _remove_excluded_hydrogens(rwmol, excluded_hydrogens)

    # find atoms with properties appended as all indices have changed now
    (
        atoms_to_highlight,
        idx2rgb,
    ) = _get_highlighted_atoms_and_strength_colors(rwmol, highlight_atoms)

    # draw the molecule and return render
    img_txt = _draw_molecule(
        rwmol=rwmol,
        fig_size=fig_size,
        atoms_to_highlight=atoms_to_highlight,
        idx2rgb=idx2rgb,
    )
    return img_txt
