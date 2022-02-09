"""RDKit functions to convert Jazzy data into rendering."""
# src/jazzy/visualisation.py
from typing import Tuple

import rdkit
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


def fix_explicit_hs(rwmol: Chem.rdchem.RWMol, idx: int):
    """Increase number of explicit Hydrogens for atom with index `idx`."""
    oa = rwmol.GetAtomWithIdx(idx)
    if oa.GetAtomicNum() > 1:
        oa.SetNumExplicitHs(oa.GetNumExplicitHs() + 1)


def remove_selected_hs(rwmol: Chem.rdchem.RWMol, hs_to_remove: list):
    """Convert list of Hydrogen atoms into explicit Hydrogens.

    Convert implicit Hydrogen atoms into explicit Hydrogen on the atoms where
    they are bonded. E.g., H-N-H is converted to NH2.

    Args:
    rwmol: An RDKit RWmolecule (rdkit.Chem.rdchem.RWMol)
    hs_to_remove: list of implicit Hydrogens to convert into explicit ones.

    Returns:
    rdkit_molecule: An RDKit molecule (rdkit.Chem.rdchem.Mol)

    """
    if not hs_to_remove:
        return rwmol

    ai_to_remove = list()
    for bond_idx in reversed(range(rwmol.GetNumBonds())):
        b = rwmol.GetBondWithIdx(bond_idx)
        bidx = b.GetBeginAtomIdx()
        remove_bidx = bidx in hs_to_remove
        eidx = b.GetEndAtomIdx()
        remove_eidx = eidx in hs_to_remove
        if remove_bidx or remove_eidx:
            if remove_bidx:
                ai_to_remove.append(bidx)
                fix_explicit_hs(rwmol, eidx)
            if remove_eidx:
                ai_to_remove.append(eidx)
                fix_explicit_hs(rwmol, bidx)
            rwmol.RemoveBond(bidx, eidx)
    for atom_idx in sorted(ai_to_remove, reverse=True):
        rwmol.RemoveAtom(atom_idx)
    rdkit.Chem.SanitizeMol(rwmol)
    return rwmol.GetMol()


def remove_hs_bonded_to_strong_acceptor(rwmol: Chem.rdchem.RWMol, hs_to_remove: list):
    """Remove Hydrogens that are bonded to strong acceptors.

    Goes through a list of Hydrogens and checks wheather any of their neighbors
    has been annotated with acceptor strength. If so, it removes the Hydrogen
    from the removal list to avoid it to collapse with the acceptor atom.

    Args:
    rwmol: An RDKit RWmolecule (rdkit.Chem.rdchem.RWMol)
    hs_to_remove: list of implicit Hydrogens to convert into explicit ones.

    Returns:
    Updated list for Hydrogen atoms to be removed.

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


def create_color_scale(idx2strength: dict, mode):
    """Create RGB mapping.

    Normalises values in a dictionary and creates an RGB red or blue scale
    mapping for them. `abs(round(0.9-alpha, 3))` is a bit hacky as it avoids
    the creation of whites by sacrificing a bit of alpha. Also RGB for acceptor
    has been tweaked a bit to avoid generating scales that might hide atoms
    rendered with blue fonts (e.g., Nitrogen atoms).

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
            alpha = round((v - amin) / (amax - amin), 3)
        except ZeroDivisionError:
            alpha = round(v, 3)
        invalpha = abs(round(0.9 - alpha, 3))
        if mode == "donor":
            rgb = (1.0, invalpha, invalpha)
        else:
            # acceptor
            rgb = (invalpha, 0.7, 1.0)
        idx2rgb[idx] = rgb
    return idx2rgb


def draw_molecule(
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
    d2d.DrawMolecule(
        rwmol,
        highlightAtoms=atoms_to_highlight,
        highlightAtomColors=idx2rgb,
        highlightBonds=None,
    )
    d2d.FinishDrawing()
    return d2d.GetDrawingText()


def create_hs_to_remove(
    rwmol: Chem.rdchem.RWMol,
    atomic_map: dict,
    sa_treshold: float,
    sdc_treshold: float,
    sdx_treshold: float,
    ignore_sa: bool,
    ignore_sdc: bool,
    ignore_sdx: bool,
):
    """Create `hs_to_remove` list from RDKit molecule.

    Args:
    rwmol: An RDKit RWmolecule (rdkit.Chem.rdchem.RWMol)
    atomic_map: molecular map of polar properties (dict)
    sa_treshold: treshold to show acceptor contribution (float)
    sdc_treshold: treshold to show donor contribution on Carbon (float)
    sdx_treshold: treshold to show donor contribution on non-Carbon (float)
    ignore_sa: ignore acceptor contributions (bool)
    ignore_sdc: ignore donor contributions on Carbon (bool)
    ignore_sdx: ignore donor contributions on non-Carbon (bool)

    Returns:
    `hs_to_remove` list.

    """
    hs_to_remove = list()
    for idx, atom in enumerate(rwmol.GetAtoms()):
        atom_props = atomic_map[idx]

        # not involved in Hydrogen bonding
        sa = round(atom_props["sa"], 3)
        sdc = round(atom_props["sdc"], 3)
        sdx = round(atom_props["sdx"], 3)
        if sa == 0.0 and sdc == 0.0 and sdx == 0.0:
            continue

        # acceptor logic
        # setting two props to avoid mistakes in the highlighting function
        if sa != 0 and sa > sa_treshold and not ignore_sa:
            atom.SetProp("atomNote", str(sa))
            atom.SetProp("sa", str(sa))
            continue

        # donor logic
        # setting two props to avoid mistakes in the highlighting function
        if sdc != 0 and sdc > sdc_treshold and not ignore_sdc:
            atom.SetProp("atomNote", str(sdc))
            atom.SetProp("sd", str(sdc))
            continue

        if sdx != 0 and sdx > sdx_treshold and not ignore_sdx:
            atom.SetProp("atomNote", str(sdx))
            atom.SetProp("sd", str(sdx))
            continue

        # Hydrogen does not fall into above conditions
        # then not strong enough and mark for removal
        if atom_props["z"] == 1:
            hs_to_remove.append(idx)

    # remove hs from hs_to_remove if connected with a strong acceptor
    if not ignore_sa:
        hs_to_remove = remove_hs_bonded_to_strong_acceptor(rwmol, hs_to_remove)

    return hs_to_remove


def depict_strength(
    rdkit_molecule: Chem.rdchem.Mol,
    atomic_map: dict,
    fig_size=(500, 500),
    flatten_molecule=False,
    highlight_atoms=False,
    ignore_sdc=False,
    ignore_sdx=False,
    ignore_sa=False,
    sdc_treshold=0.0,
    sdx_treshold=0.0,
    sa_treshold=0.0,
):
    """Create an SVG image text from an RDKit molecule and its atomic map.

    The default configuration simply produces a depiction of the input molecule
    and its strengths. `highlight_atoms` highlights atoms in red (donors) and
    blue (acceptors). `flatten_molecule` produces a 2-dimensional depiction of
    the molecule. Any `treshold` parameter allows to set a numeric treshold
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
    sdc_treshold: treshold to show donor contribution on Carbon (float)
    sdx_treshold: treshold to show donor contribution on non-Carbon (float)
    sa_treshold: treshold to show acceptor contribution (float)

    Returns:
    SVG depiction of given RDKit molecule and its atomic strength map.

    """
    # copy RDKit molecule otherwise the input would be affected by processing
    rwmol = rdkit.Chem.RWMol(rdkit_molecule)

    # flatten the molecule if required
    # Note: never create SMILES then again RDKit 2D molecule since the mapping
    # could be shuffled! Only recompute 2D coordinates on the 3D molecule.
    if flatten_molecule:
        rdkit.Chem.rdDepictor.Compute2DCoords(rwmol)

    hs_to_remove = create_hs_to_remove(
        rwmol,
        atomic_map,
        sa_treshold,
        sdc_treshold,
        sdx_treshold,
        ignore_sa,
        ignore_sdc,
        ignore_sdx,
    )

    # remove Hydrogen atoms that are not strong enough
    rwmol = remove_selected_hs(rwmol, hs_to_remove)

    # find atoms with properties appended as all indices have changed now
    idx2sd = dict()
    idx2sa = dict()
    idx2rgb = dict()
    atoms_to_highlight = list()
    if highlight_atoms:
        atoms = [x for x in rwmol.GetAtoms()]
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
        idx2reds = create_color_scale(idx2sd, mode="donor")
        idx2blues = create_color_scale(idx2sa, mode="acceptor")
        idx2rgb = {**idx2reds, **idx2blues}  # type: ignore

    # draw the molecule and return render
    img_txt = draw_molecule(
        rwmol=rwmol,
        fig_size=fig_size,
        atoms_to_highlight=atoms_to_highlight,
        idx2rgb=idx2rgb,
    )
    return img_txt
