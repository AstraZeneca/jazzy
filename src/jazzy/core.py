"""Core functions of the jazzy package."""
# src/jazzy/core.py
from typing import Optional

import numpy as np
from kallisto.atom import Atom
from kallisto.methods import getPolarizabilities
from kallisto.methods import getVanDerWaalsRadii
from kallisto.molecule import Molecule
from kallisto.units import Bohr
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import GetPeriodicTable
from rdkit.Chem import PeriodicTable
from rdkit.Chem import rdchem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import Mol

from jazzy.config import EMBEDDING_TYPES
from jazzy.config import ROUNDING_DIGITS
from jazzy.exception import EmbeddingError
from jazzy.exception import KallistoError
from jazzy.exception import NegativeLonePairsError
from jazzy.logging import logger
from jazzy.model import ChargeKey
from jazzy.model import IdxKey


def rdkit_molecule_from_smiles(
    smiles: str, minimisation_method=None, **kwargs
) -> Optional[Mol]:
    """Molecule preparation: Parse SMILES, add hydrogens, and does energy minimisation.

    Args:
        smiles: A molecule SMILES string representation (default '')
        minimisation_method: One of the conformer energy minimisation methods
            as available in RDKit (available is 'MMFF94', 'MMFF94s', or 'UFF')
            (default None)
        kwargs: Keyword arguments

    Keyword Args:
        embedding_type: Molecule embedding method (available as '2D' or '3D')
            (default '3D')
        embedding_seed: Integer seed for the embedding process (default 11)
        embedding_max_iterations: Maximum number of iterations for the embedding

    Returns:
        An RDKit molecule (rdkit.Chem.rdchem.Mol) or None if the process fails

    """
    # create molecule, add hydrogens
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        logger.error("The RDKit SMILES parsing has failed for the molecule: %s", smiles)
        return None
    mh = Chem.AddHs(m)

    # generate embedding
    try:
        _embed_molecule(mh, **kwargs)
    except EmbeddingError:
        return None

    # energy minimisation
    _minimise_molecule(mh, minimisation_method, **kwargs)
    return mh


def _embed_molecule(rdkit_molecule: Mol, **kwargs) -> None:
    """Molecule embedding using either the 2D or 3D method."""
    embedding_type = kwargs.get("embedding_type", "3D")
    embedding_seed = kwargs.get("embedding_seed", 11)
    if embedding_type not in EMBEDDING_TYPES:
        raise ValueError(f"Select a valid embedding type {EMBEDDING_TYPES}")
    if embedding_type == EMBEDDING_TYPES[0]:
        AllChem.Compute2DCoords(rdkit_molecule, sampleSeed=embedding_seed)
    if embedding_type == EMBEDDING_TYPES[1]:
        params = (
            rdDistGeom.ETDG()
        )  # ETDG over others as it seems produce fewer failures
        params.randomSeed = embedding_seed
        params.maxIterations = kwargs.get("embedding_max_iterations", 0)
        emb_code = AllChem.EmbedMolecule(rdkit_molecule, params)
        if emb_code == -1:
            error_msg = (
                "The RDKit embedding has failed for the molecule: "
                f"{Chem.MolToSmiles(rdkit_molecule)} (Seed: {embedding_seed}, "
                f"Max Iterations: {params.maxIterations})"
            )
            logger.error(error_msg)
            raise EmbeddingError(error_msg)


def _minimise_molecule(rdkit_molecule: Mol, minimisation_method: str, **kwargs) -> None:
    """Minimise the energy of a molecule. Only works for 3D molecules."""
    valid_methods = (None, "MMFF94", "MMFF94s", "UFF")
    if not minimisation_method:
        return
    embedding_type = kwargs.get("embedding_type", "3D")
    if embedding_type == EMBEDDING_TYPES[0]:
        logger.warning("Energy minimisation is not supported for 2D molecules.")
        return
    if minimisation_method == valid_methods[1]:
        AllChem.MMFFOptimizeMolecule(rdkit_molecule)
    elif minimisation_method == valid_methods[2]:
        AllChem.MMFFOptimizeMolecule(rdkit_molecule, mmffVariant="MMFF94s")
    elif minimisation_method == valid_methods[3]:
        AllChem.UFFOptimizeMolecule(rdkit_molecule)
    else:
        raise ValueError(f"Select a valid minimisation method {valid_methods}")


def get_covalent_atom_idxs(rdkit_molecule: Mol) -> list:
    """Get covalent indices for atom_idx.

    Creates a list of lists, where indices are atom indices and the content
    of the lists are indices of binding atoms (e.g., [[1],[0]] means that atom
    0 is bound to atom 1 ([1]) and vice versa ([0]).

    """
    atoms_and_idxs = list()
    for _, atom in enumerate(rdkit_molecule.GetAtoms()):
        nbrs = atom.GetNeighbors()
        atoms_and_idxs.append(list(map(lambda x: x.GetIdx(), nbrs)))
    return atoms_and_idxs


def get_all_neighbours(rdkit_molecule: Mol, molecule_covalent_nbrs: list):
    """Get all alpha, beta, and gamma neighbours.

    Create dictionaries for covalent bonding partner (alpha), all nearest
    neighbours (beta), and nearest-nearest (gamma) neighbours. Every dictionary
    contains the atomic index within the molecule as key and a list of
    neighbours as the value.

    Args:
        rdkit_molecule: RDKit molecule
        molecule_covalent_nbrs: List of lists containing covalent neighbours

    Returns:
        Dictionaries for alpha, beta, and gamma neighbours

    """
    alpha = dict()
    beta = dict()
    gamma = dict()
    for i, _ in enumerate(rdkit_molecule.GetAtoms()):
        atom_idxs_dict = get_atom_and_nbrs_idxs_dict(i, molecule_covalent_nbrs)
        alpha[i] = atom_idxs_dict[IdxKey.Alpha]
        beta[i] = atom_idxs_dict[IdxKey.Beta]
        gamma[i] = atom_idxs_dict[IdxKey.Gamma]
    return alpha, beta, gamma


def kallisto_molecule_from_rdkit_molecule(rdkit_molecule: Mol) -> Molecule:
    """Create a kallisto molecule from RDKit molecule.

    Args:
        rdkit_molecule: RDKit molecule

    Returns:
        A kallisto molecule (kallisto.molecule.Molecule)

    Raises:
        KallistoError: An error if the kallisto molecule cannot be created

    """
    # get the name of the molecule if it comes from SDF
    name = ""
    if rdkit_molecule.HasProp("_Name"):
        name = rdkit_molecule.GetProp("_Name")
    # get all xyz coordinates and split into list of lines
    xyz = Chem.rdmolfiles.MolToXYZBlock(rdkit_molecule).split("\n")
    # remove empty lines or molecule name from list
    xyz = [string for string in xyz if string != "" and string != name]
    # remove number of atoms as given in xmol files (first line)
    xyz = xyz[1:]

    # setup periodic table
    pt = GetPeriodicTable()
    # create list of atoms
    atoms = []
    # create kallisto molecule
    for coord in xyz:
        elem, x, y, z = coord.split()[:4]

        # convert atomic coordinates from Angstrom to Bohr
        position = [float(x) / Bohr, float(y) / Bohr, float(z) / Bohr]
        atom = Atom(symbol=pt.GetAtomicNumber(elem), position=position)
        atoms.append(atom)
    kallisto_mol = Molecule(symbols=atoms)
    if "numbers" not in kallisto_mol.arrays.keys():
        raise KallistoError(
            "The kallisto molecule was not created for the input '{}'".format(
                Chem.MolToSmiles(rdkit_molecule)
            )
        )
    return kallisto_mol


def get_charges_from_kallisto_molecule(
    kallisto_molecule: Molecule, charge: int
) -> list:
    """Calculate electronegativity equilibration (EEQ) atomic partial charges.

    Args:
        kallisto_molecule: kallisto molecule
        charge: molecular charge (int)

    Returns:
        List of electronegativity equilibration atomic partial charges

    Raises:
        KallistoError: If 'charges' contains any NaNs.
    """
    charges = list(kallisto_molecule.get_eeq(charge=charge))
    if np.isnan(charges).any():
        raise KallistoError("EEQ calculation resulted in NaNs")
    return charges


def calculate_polar_strength_map(
    rdkit_molecule: Mol,
    kallisto_molecule: Molecule,
    atoms_and_nbrs: list,
    charges: list,
    d=6.1475,
    a=-2.2316,
    t=0.274,
) -> dict:
    """Calculate the polar strength map.

    Generates a molecular dictionary where keys correspond to atom indices,
    and values are dictionaries of calculated properties. These properties
    include polar acceptor strength (sa), polar donor strength of Hydrogens
    bonded to Carbon (sdc) or other non-Hydrogen atoms (sdx), number of lone
    pairs (num_lp), atomic number (z), and formal charge (q).

    Args:
        rdkit_molecule: RDKit molecule
        kallisto_molecule: kallisto molecule
        atoms_and_nbrs: list of lists of atoms and their bonded atoms
        charges: list of atomic partial charges
        d: parameter that ensures that the donor strength (sdx) of
            Hydrogen in H2O is equal to 1.000
        a: parameter that ensures that the acceptor strength (sa) of
            one lone pair on Oxygen in H2O is equal to 1.000
        t: bond reduction factor (fixed to ensure conditions above)

    Returns:
        Map of molecular polar properties (plus some extras)

    """
    # atomic numbers (e.g., C=6, H=1)
    ats = [atom.GetAtomicNum() for atom in rdkit_molecule.GetAtoms()]
    # atomic covalent coordination numbers
    cns = kallisto_molecule.get_cns(cntype="cov")
    # atomic static polarizabilities
    alps = getPolarizabilities(ats, cns, charges, 0)  # type: ignore

    mol_map = dict()
    for idx, atom in enumerate(rdkit_molecule.GetAtoms()):
        sa = 0.0
        sdx = 0.0
        sdc = 0.0
        nlps = 0
        z = ats[idx]
        q = atom.GetFormalCharge()
        eeq = charges[idx]
        hyb = atom.GetHybridization().name.lower()
        alp = alps[idx]

        # if Hydrogen -> evaluate donor strength
        if z == 1:
            sd = get_donor_atom_strength(idx, atoms_and_nbrs, charges, d, t)

            # distinguish between Carbon and other non-Hydrogen atoms
            nbr = atom.GetNeighbors()
            if nbr:
                # Carbon case
                if nbr[0].GetAtomicNum() == 6:
                    sdc = sd
                # other non-Hydrogen case
                else:
                    sdx = sd
        else:
            nlps = get_lone_pairs(atom)
            if nlps > 0:
                sa = get_acceptor_atom_strength(idx, atoms_and_nbrs, charges, a, t)

        # create dict where key is atom index and values are properties
        atom_dict = {
            "z": z,
            "q": q,
            "eeq": round(eeq, ROUNDING_DIGITS),
            "alp": round(alp, ROUNDING_DIGITS),
            "hyb": hyb,
            "num_lp": nlps,
            "sdc": round(sdc, ROUNDING_DIGITS),
            "sdx": round(sdx, ROUNDING_DIGITS),
            "sa": round(sa, ROUNDING_DIGITS),
        }

        mol_map[idx] = atom_dict
    return mol_map


def get_lone_pairs(atom) -> int:
    """Get the number of lone pairs for an atom.

    The method is similar to that in
    https://github.com/rdkit/blob/master/Code/GraphMol/Aromaticity.cpp with
    the exception that it calculates the explicit valence via atom functions
    (not PeriodicTable). This is because the original method produces
    miscalculations for some systems (e.g., in molecules with SO2 groups).

    """
    # set up a periodic table
    pt = Chem.GetPeriodicTable()
    symbol = atom.GetSymbol()
    valence_electrons = PeriodicTable.GetNOuterElecs(pt, symbol)
    unavailable_electrons = atom.GetExplicitValence()
    charge = atom.GetFormalCharge()
    free_electrons = valence_electrons - unavailable_electrons - charge
    return int(free_electrons / 2)


def get_atom_and_nbrs_idxs_dict(atom_idx: int, molecule_covalent_nbrs: list):
    """Get all neighbours for atom_idx.

    Extract all covalent bonding partner (alpha), all nearest neighbours
    (beta), and all nearest-nearest neighbours (gamma) for atom with index
    'atom_idx'.

    Args:
        atom_idx: index of atom to extract neighbours for
        molecule_covalent_nbrs: List of lists containing covalent neighbours

    Returns:
        alpha: list of all covalent bonding atom indices of atom_idx
        beta: list of nearest neighbour atom indices of atom_idx
        gamma: list of nearest-nearest neighbour atom indices of atom_idx

    """
    # extract alpha neighbours
    alpha_idxs = molecule_covalent_nbrs[atom_idx]

    # extract beta neighbours
    beta_idxs = list()
    for _, a in enumerate(alpha_idxs):
        b = molecule_covalent_nbrs[a]
        diff = list(set([atom_idx]) ^ set(b))
        if len(diff) > 0:
            beta_idxs.extend(diff)

    # extract gamma neighbours
    gamma_idxs = list()
    for _, b in enumerate(beta_idxs):
        c = molecule_covalent_nbrs[b]
        inter = list(set(alpha_idxs).intersection(set(c)))
        diff = list(set(inter) ^ set(c))
        gamma_idxs.extend(diff)
    gamma_idxs = list(dict.fromkeys(gamma_idxs))
    return {
        IdxKey.Atom: atom_idx,
        IdxKey.Alpha: alpha_idxs,
        IdxKey.Beta: beta_idxs,
        IdxKey.Gamma: gamma_idxs,
    }


def get_charges_from_atom_list(atom_idxs: list, charges: list) -> list:
    """List-based function for charge retrieval."""
    return _get_value_from_list(atom_idxs, charges)


def _get_value_from_list(idxs: list, values: list) -> list:
    """Retrieve values from list based on indices."""
    return [values[idx] for idx in idxs]


def get_donor_atom_strength(
    atom_idx: int, atoms_and_nbrs: list, charges: list, d=63.7, t=0.274
) -> float:
    """Donor strength calculation - equation 10."""
    q, q_delta = calculate_q_and_delta_q(atom_idx, atoms_and_nbrs, charges, t)
    return d * (q + q_delta)


def get_acceptor_atom_strength(
    atom_idx: int, atoms_and_nbrs: list, charges: list, a=-4.4362, t=0.274
) -> float:
    """Acceptor strength calculation - equation 12 and 13."""
    q, q_delta = calculate_q_and_delta_q(atom_idx, atoms_and_nbrs, charges, t)
    return a * (q + q_delta)


def calculate_q_and_delta_q(
    atom_idx: int, atoms_and_nbrs: list, charges: list, t=0.274
):
    """Calculates charge and delta charge.

    Suitable for both donor and acceptor strength calculations.

    """
    # get idxs and charges
    atom_idxs_dict, q_dict = _get_idxs_and_q_dicts(atom_idx, atoms_and_nbrs, charges)

    # calculate q delta
    q_delta = _calculate_q_delta(
        q_dict[ChargeKey.Alpha], q_dict[ChargeKey.Beta], q_dict[ChargeKey.Gamma], t
    )
    return q_dict[ChargeKey.Atom], q_delta


def _get_idxs_and_q_dicts(atom_idx: int, atoms_and_nbrs: list, charges: list):
    """Wraps get_atom_and_nbrs_idxs_dict and _get_q_dict."""
    # get indices
    atom_idxs_dict = get_atom_and_nbrs_idxs_dict(atom_idx, atoms_and_nbrs)

    # get charges
    q_dict = _get_q_dict(atom_idxs_dict, charges)
    return atom_idxs_dict, q_dict


def _get_q_dict(atom_idxs_dict, charges):
    """Generates a dictionary of list of charges for each group of atoms."""
    q_dict = dict()
    q_dict[ChargeKey.Atom] = get_charges_from_atom_list(
        [atom_idxs_dict[IdxKey.Atom]], charges
    )[0]
    q_dict[ChargeKey.Alpha] = get_charges_from_atom_list(
        atom_idxs_dict[IdxKey.Alpha], charges
    )
    q_dict[ChargeKey.Beta] = get_charges_from_atom_list(
        atom_idxs_dict[IdxKey.Beta], charges
    )
    q_dict[ChargeKey.Gamma] = get_charges_from_atom_list(
        atom_idxs_dict[IdxKey.Gamma], charges
    )
    return q_dict


def _calculate_q_delta(q_alpha, q_beta, q_gamma, t):
    """Calculates q delta as per Equation 9 in Gerber's paper."""
    return t * sum(q_alpha) + (t**2) * sum(q_beta) + (t**3) * sum(q_gamma)


def calculate_delta_apolar(
    rdkit_molecule: Mol,
    mol_map: dict,
    g0: float,
    gs: float,
    gr: float,
    gpi1: float,
    gpi2: float,
) -> float:
    """Calculate apolar free energy contribution.

    Calculate the apolar contribution to the free energy.
    Equation 7-8 - parameter names are adapted.

    Args:
        rdkit_molecule: RDKit molecule
        mol_map: molecular map of polar properties
        g0: zeroth order parameter
        gs: surface parameter
        gr: ring parameter
        gpi1: first pi parameter
        gpi2: second pi parameter

    Returns:
        Apolar free energy contribution (float)

    """
    # get atoms from rdkit molecule
    atoms = rdkit_molecule.GetAtoms()
    # atomic numbers (e.g., C=6, H=1)
    at = [atom.GetAtomicNum() for atom in atoms]
    # hybridization states
    hi = [atom.GetHybridization() for atom in atoms]
    # ring count
    ring = rdMolDescriptors.CalcNumRings(rdkit_molecule)
    # non-hydrogen count
    idxs = [idx for idx, nuclear_charge in enumerate(at) if nuclear_charge != 1]
    # atomic dynamic polarizabilities
    alp = [atom["alp"] for atom in mol_map.values()]
    # van der Waals radii in Angstrom
    vdws = getVanDerWaalsRadii(
        len(at),
        at,  # type: ignore
        alp,  # type: ignore
        vdwtype="rahm",
        scale=Bohr,
    )

    ns = 0
    nl = dict()
    for _, idx in enumerate(idxs):
        # non-hydrogen partner count
        nl[idx] = 0
        vdw = vdws[idx]
        atom = atoms[idx]
        nbrs = atom.GetNeighbors()
        for nbr in nbrs:
            if nbr.GetAtomicNum() != 1:
                nl[idx] += 1

        # sphere term ns
        bracket = 1 - (nl[idx]) / (hi[idx] + 1)
        ns += 4 * np.pi * np.power(vdw, 2) * bracket

    # zeroth order term
    dga = g0
    # surface term
    dga += gs * ns
    # ring term
    dga += gr * ring

    # identify pi orbitals in hybridization states
    # SP  : 2 p-orbitals
    # SP2 : 1 p-orbital
    sp, sp2 = 0, 0
    for _, v in enumerate(hi):
        if v == rdchem.HybridizationType.SP:
            sp += 2
        elif v == rdchem.HybridizationType.SP2:
            sp2 += 1

    # pi-orbital terms
    dga += sp * gpi1
    dga += sp2 * gpi2
    return float(dga)


def calculate_delta_polar(
    mol_map: dict,
    atoms_and_nbrs: list,
    gd: float,
    ga: float,
    expd: float,
    expa: float,
) -> float:
    """Calculate polar contribution to free energy.

    Calculate the polar contribution to the free energy.
    Equation 14 - parameter names adapted.

    Args:
        mol_map: molecular map of polar properties as obtained by
            the `calculate_polar_strength_map()` function
        atoms_and_nbrs: list of lists of atoms and their bonded atoms
        gd: hydrogen bond donor parameter
        ga: hydrogen bond acceptor parameter
        expd: exponent for hydrogen bond donor
        expa: exponent for hydrogen bond acceptor

    Returns:
        Polar free energy contribution (float)

    Raises:
        NegativeLonePairsError: if the input compound contains atoms with
            negative number of lone pairs (lone pairs < 0)

    """
    don = 0.0
    acc = 0.0

    for _, idx in enumerate(mol_map):
        atom = mol_map[idx]
        nh = 0
        sdi = 0.0
        nlps = 0
        sak = 0.0

        # if atom is Hydrogen -> skip
        # get contribution from heavy partner atom
        if atom["z"] == 1:
            continue

        # donor contribution
        nbrs_idxs = atoms_and_nbrs[idx]
        for _, nbr_idx in enumerate(nbrs_idxs):
            nbr = mol_map[nbr_idx]

            # if neighbor is Hydrogen
            # get strength and add to total
            if nbr["z"] == 1:
                sd = nbr["sdx"] + nbr["sdc"]
                sdi += sd
                nh += 1

        # mean donor strength
        if nh > 0:
            sdi = sdi / nh

        # acceptor contribution
        sak = atom["sa"]
        nlps = atom["num_lp"]
        if nlps < 0:
            raise NegativeLonePairsError(
                "The input compound contains atoms with "
                f"negative lone pairs (got {nlps} for atom at idx {idx})"
            )

        don += sdi * (nh**expd)
        acc += sak * (nlps**expa)
    # equation 14
    return float(gd * don + ga * acc)


def any_hydrogen_neighbors(rdkit_atom: Chem.rdchem.Atom):
    """Returns True is Hydrogen is a neighbor else False.

    Written because `rdkit.Chem.rdchem.Atom.GetTotalHs()` does not work on
    embedded molecules.

    Args:
        rdkit_atom: rdkit.Chem.rdchem.Atom object

    Returns:
        Boolean

    """
    for nbr in rdkit_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 1:
            return True
    return False


def interaction_strength(idx: int, mol_map: dict, acceptor_exp: float) -> float:
    """Calculate interaction strength for atom with index `idx`."""
    acceptor_strength = mol_map[idx]["sa"]
    num_lp = mol_map[idx]["num_lp"]
    if num_lp != 0:
        return acceptor_strength * (num_lp**acceptor_exp)
    return 0.0


def calculate_delta_interaction(
    rdkit_molecule: Mol,
    mol_map: dict,
    atoms_and_nbrs: list,
    gi: float,
    expa: float,
    f: float,
):
    """Calculate interaction contribution to free energy.

    Calculate the interation contribution to the free energy.
    Equation 15 - parameter names adapted.

    Args:
        rdkit_molecule: RDKit molecule
        mol_map: molecular map of polar properties as obtained by
        atoms_and_nbrs: list of lists of atoms and their bonded atoms
        gi: interaction parameter
        expa: exponent hydrogen bond acceptor
        f: correction parameter

    Returns:
        Interaction free energy contribution (float)

    """
    # extract neighbors
    alpha, beta, gamma = get_all_neighbours(rdkit_molecule, atoms_and_nbrs)

    dgi = 0.0
    atoms = rdkit_molecule.GetAtoms()
    for idx, _ in enumerate(atoms):
        # strength of origin
        s_origin = interaction_strength(idx, mol_map, expa)

        # if any Hydrogen is available continue
        if any_hydrogen_neighbors(atoms[idx]):
            continue

        # alpha neighbor contributions
        for _, va in enumerate(alpha[idx]):
            if any_hydrogen_neighbors(atoms[va]):
                continue

            # strength of partner
            s_partner = interaction_strength(va, mol_map, expa)
            dgi += s_origin * s_partner

        # beta neighbor contribution
        for _, vb in enumerate(beta[idx]):
            if any_hydrogen_neighbors(atoms[vb]):
                continue

            # strength of partner
            s_partner = interaction_strength(vb, mol_map, expa)
            dgi += s_origin * s_partner

        # gamma neighbor contribution
        for _, vc in enumerate(gamma[idx]):
            if any_hydrogen_neighbors(atoms[vc]):
                continue

            # strength of partner
            s_partner = interaction_strength(vc, mol_map, expa)
            dgi += f * s_origin * s_partner
    return gi * dgi
