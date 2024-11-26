"""Application programming interface for the jazzy package."""
# src/jazzy/api.py
import base64

from jazzy.config import Config
from jazzy.config import MINIMISATION_METHOD
from jazzy.config import ROUNDING_DIGITS
from jazzy.core import calculate_delta_apolar
from jazzy.core import calculate_delta_interaction
from jazzy.core import calculate_delta_polar
from jazzy.core import calculate_polar_strength_map
from jazzy.core import get_charges_from_kallisto_molecule
from jazzy.core import get_covalent_atom_idxs
from jazzy.core import kallisto_molecule_from_rdkit_molecule
from jazzy.core import rdkit_molecule_from_smiles
from jazzy.exception import exception_handling
from jazzy.exception import JazzyError
from jazzy.helpers import condense_atomic_map
from jazzy.helpers import convert_map_to_tuples
from jazzy.helpers import sum_atomic_map
from jazzy.visualisation import depict_strengths


# global jazzy config (parameter)
config = Config()


def __smiles_to_molecule_objects(
    smiles, minimisation_method=MINIMISATION_METHOD, **kwargs
):
    """Private method for converting SMILES into RDKit and kallisto objects."""
    if smiles == "":
        raise JazzyError("An empty SMILES string was passed.")
    rdkit_mol = rdkit_molecule_from_smiles(
        smiles, minimisation_method=minimisation_method, **kwargs
    )
    if rdkit_mol is None:
        raise JazzyError("The SMILES '{}' could not be processed.".format(smiles))
    kallisto_mol = kallisto_molecule_from_rdkit_molecule(rdkit_mol)
    return rdkit_mol, kallisto_mol


@exception_handling
def molecular_vector_from_smiles(
    smiles: str, minimisation_method=MINIMISATION_METHOD, only_strengths=False, **kwargs
):
    """API route to calculate molecular free energy vector.

    Calculates the apolar (dga), the polar (dgp), and the interaction (dgi)
    contribution to the free energy.

    Args:
        smiles: A molecule SMILES string representation (default '')
        minimisation_method: One of the conformer energy minimisation methods
            as available in RDKit (available as 'MMFF94', 'MMFF94s', or 'UFF')
             (default None)
        only_strengths: Boolean value that determines wheather to calculate only
            strengts or even more
        kwargs: Keyword arguments

    Keyword Args:
        embedding_type: Molecule embedding method (available as '2D' or '3D')
            (default '3D')
        embedding_seed: Integer seed for the embedding process (default 11)
        embedding_max_iterations: Maximum number of iterations for the embedding


    Returns:
        Molecular strength vector with or without free energy contributions

    """
    mol_vector = dict()
    # generate an RDKit molecule
    rdkit_molecule, kallisto_molecule = __smiles_to_molecule_objects(
        smiles, minimisation_method, **kwargs
    )
    atoms_and_nbrs = get_covalent_atom_idxs(rdkit_molecule)
    kallisto_charges = get_charges_from_kallisto_molecule(kallisto_molecule, 0)
    atomic_map = calculate_polar_strength_map(
        rdkit_molecule, kallisto_molecule, atoms_and_nbrs, kallisto_charges
    )
    mol_vector = sum_atomic_map(atomic_map)

    # add free energy contributions
    if not only_strengths:
        dg = dict()
        dg["dga"] = calculate_delta_apolar(
            rdkit_molecule,
            atomic_map,
            config.g0,
            config.gs,
            config.gr,
            config.gpi1,
            config.gpi2,
        )
        dg["dgp"] = calculate_delta_polar(
            atomic_map, atoms_and_nbrs, config.gd, config.ga, config.expd, config.expa
        )
        dgi = calculate_delta_interaction(
            rdkit_molecule, atomic_map, atoms_and_nbrs, config.gi, config.expa, config.f
        )
        dg["dgtot"] = sum(dg.values()) + dgi
        dg = {k: round(dg[k], ROUNDING_DIGITS) for k in dg}
        mol_vector = {**mol_vector, **dg}  # type: ignore
    return mol_vector


@exception_handling
def deltag_from_smiles(smiles: str, minimisation_method=MINIMISATION_METHOD, **kwargs):
    """API route to calculate molecular free energy scalar.

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
        Free energy as scalar rounded

    """
    # generate basic descriptors
    rdkit_molecule, kallisto_molecule = __smiles_to_molecule_objects(
        smiles, minimisation_method, **kwargs
    )
    atoms_and_nbrs = get_covalent_atom_idxs(rdkit_molecule)
    kallisto_charges = get_charges_from_kallisto_molecule(kallisto_molecule, 0)
    atomic_map = calculate_polar_strength_map(
        rdkit_molecule, kallisto_molecule, atoms_and_nbrs, kallisto_charges
    )

    # generate free energy scalar
    dg = dict()
    dg["dga"] = calculate_delta_apolar(
        rdkit_molecule,
        atomic_map,
        config.g0,
        config.gs,
        config.gr,
        config.gpi1,
        config.gpi2,
    )
    dg["dgp"] = calculate_delta_polar(
        atomic_map, atoms_and_nbrs, config.gd, config.ga, config.expd, config.expa
    )
    dg["dgi"] = calculate_delta_interaction(
        rdkit_molecule, atomic_map, atoms_and_nbrs, config.gi, config.expa, config.f
    )
    return round(sum(dg.values()), ROUNDING_DIGITS)


@exception_handling
def atomic_tuples_from_smiles(
    smiles: str, minimisation_method=MINIMISATION_METHOD, **kwargs
):
    """API route to generate a tuple representation on the atomic map.

    Not recommended if serialization is needed.

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
        Tuple representation of the atomic map.

    """
    # generate basic descriptors
    rdkit_molecule, kallisto_molecule = __smiles_to_molecule_objects(
        smiles, minimisation_method, **kwargs
    )
    atoms_and_nbrs = get_covalent_atom_idxs(rdkit_molecule)
    kallisto_charges = get_charges_from_kallisto_molecule(kallisto_molecule, 0)
    atomic_map = calculate_polar_strength_map(
        rdkit_molecule, kallisto_molecule, atoms_and_nbrs, kallisto_charges
    )
    return convert_map_to_tuples(atomic_map)


@exception_handling
def atomic_map_from_smiles(
    smiles: str, minimisation_method=MINIMISATION_METHOD, **kwargs
):
    """API route to generate a condensed representation on the atomic map.

    Recommended if serialization is needed.

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
        Condensed representation of the atomic map.

    """
    # generate basic descriptors
    rdkit_molecule, kallisto_molecule = __smiles_to_molecule_objects(
        smiles, minimisation_method, **kwargs
    )
    atoms_and_nbrs = get_covalent_atom_idxs(rdkit_molecule)
    kallisto_charges = get_charges_from_kallisto_molecule(kallisto_molecule, 0)
    atomic_map = calculate_polar_strength_map(
        rdkit_molecule, kallisto_molecule, atoms_and_nbrs, kallisto_charges
    )
    return condense_atomic_map(atomic_map)


@exception_handling
def atomic_strength_vis_from_smiles(
    smiles: str,
    minimisation_method=MINIMISATION_METHOD,
    encode=False,
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
    **kwargs
):
    """API route to generate an SVG image from SMILES string.

    Args:
        smiles: A molecule SMILES string representation (default '')
        minimisation_method: One of the conformer energy minimisation methods
            as available in RDKit (available is 'MMFF94', 'MMFF94s', or 'UFF')
             (default None)
        encode: If True, returns the base64-encoded SVG image (default False)
        fig_size: Tuple representing the size of the generated image
             (default (500, 500))
        flatten_molecule: If True, flattens the molecule (default False)
        highlight_atoms: If True, highlights the atoms (default False)
        ignore_sa: Ignore acceptor contributions (default False)
        ignore_sdc: Ignore Carbon-donor contributions (default False)
        ignore_sdx: Ignore Heteroatom-donor contributions (default False)
        rounding_digits: Number of digits to round the strengths
        sa_threshold: Acceptor threshold (default 0.0)
        sdc_threshold: Carbon-donor threshold (default 0.0)
        sdx_threshold: Heteroatom-donor threshold (default 0.0)
        kwargs: Keyword arguments


    Keyword Args:
        embedding_type: Molecule embedding method (available as '2D' or '3D')
            (default '3D')
        embedding_seed: Integer seed for the embedding process (default 11)
        embedding_max_iterations: Maximum number of iterations
             for the embedding


    Returns:
        SVG image either 2D or 3D.

    """
    # generate basic descriptors
    rdkit_molecule, kallisto_molecule = __smiles_to_molecule_objects(
        smiles, minimisation_method, **kwargs
    )
    atoms_and_nbrs = get_covalent_atom_idxs(rdkit_molecule)
    kallisto_charges = get_charges_from_kallisto_molecule(kallisto_molecule, 0)
    atomic_map = calculate_polar_strength_map(
        rdkit_molecule, kallisto_molecule, atoms_and_nbrs, kallisto_charges
    )
    img_txt = depict_strengths(
        rdkit_molecule=rdkit_molecule,
        atomic_map=atomic_map,
        fig_size=fig_size,
        flatten_molecule=flatten_molecule,
        highlight_atoms=highlight_atoms,
        ignore_sdc=ignore_sdc,
        ignore_sdx=ignore_sdx,
        ignore_sa=ignore_sa,
        sdc_threshold=sdc_threshold,
        sdx_threshold=sdx_threshold,
        sa_threshold=sa_threshold,
        rounding_digits=rounding_digits,
    )
    if encode:
        img_txt = base64.b64encode(img_txt.encode("utf-8"))
    return img_txt
