"""Main jazzy program."""
from typing import Tuple

import click

from jazzy.api import atomic_strength_vis_from_smiles
from jazzy.api import molecular_vector_from_smiles
from jazzy.config import Config
from jazzy.utils import JazzyError

pass_config = click.make_pass_decorator(Config, ensure=True)


@click.group()
def cli():
    """Command line interface of Jazzy."""


# CLI subcommands


@cli.command("vec")
@pass_config
@click.option(
    "--opt",
    default=None,
    type=str,
    show_default=True,
    help="Optimisation method (None, MMFF94, MMFF94s, or UFF).",
)
@click.option("--strength_only", is_flag=True)
@click.argument("smiles", type=str, default=None, required=True)
def vec(config, smiles: str, opt: str, strength_only: bool) -> None:
    """Calculate molecular descriptors."""
    # raise ValueError if minimisation method is not valid
    valid_methods = [None, "MMFF94", "MMFF94s", "UFF"]
    if opt not in valid_methods:
        raise ValueError(
            f"Minimisation method within '--opt' flag is not valid {valid_methods}"
        )

    # write SMILE to global config
    config.smiles = smiles

    # obtain molecular vector via API
    mol_vector = dict()
    try:
        mol_vector = molecular_vector_from_smiles(
            smiles=smiles,
            minimisation_method=opt,
            only_strengths=strength_only,
        )
        mol_vector["__status"] = "success"
        mol_vector["smiles"] = config.smiles
    except JazzyError:
        mol_vector["__status"] = "error"
        mol_vector["smiles"] = config.smiles
    print(mol_vector)


@cli.command("vis")
@pass_config
@click.option(
    "--opt",
    default=None,
    type=str,
    show_default=True,
    help="Optimisation method (None, MMFF94, MMFF94s, or UFF).",
)
@click.option(
    "--fig_size",
    default=(500, 500),
    type=(int, int),
    show_default=True,
    help="Size of SVG image in pixels.",
)
@click.option(
    "--sdc_threshold",
    default=0.0,
    type=float,
    show_default=True,
    help="Treshold strength to depic Carbon donors.",
)
@click.option(
    "--sdx_threshold",
    default=0.0,
    type=float,
    show_default=True,
    help="Treshold strength to depic non-Carbon donors.",
)
@click.option(
    "--sa_threshold",
    default=0.0,
    type=float,
    show_default=True,
    help="Treshold strength to depic acceptors.",
)
@click.option("--base64", is_flag=True)
@click.option("--flatten_molecule", is_flag=True)
@click.option("--highlight_atoms", is_flag=True)
@click.option("--ignore_sdc", is_flag=True)
@click.option("--ignore_sdx", is_flag=True)
@click.option("--ignore_sa", is_flag=True)
@click.argument("smiles", type=str, default=None, required=True)
def vis(
    config,
    smiles: str,
    opt: str,
    fig_size: Tuple[int, int],
    base64: bool,
    flatten_molecule: bool,
    highlight_atoms: bool,
    ignore_sdc: bool,
    ignore_sdx: bool,
    ignore_sa: bool,
    sdc_threshold: float,
    sdx_threshold: float,
    sa_threshold: float,
):
    """Create SVG image."""
    # raise ValueError if minimisation method is not valid
    valid_methods = [None, "MMFF94", "MMFF94s", "UFF"]
    if opt not in valid_methods:
        raise ValueError(
            f"Minimisation method within '--opt' flag is not valid {valid_methods}"
        )

    # write SMILE to global config
    config.smiles = smiles

    # obtain visualisation vector via API
    vis_vector = dict()
    try:
        svg = atomic_strength_vis_from_smiles(
            smiles=smiles,
            minimisation_method=opt,
            encode=base64,
            fig_size=fig_size,
            flatten_molecule=flatten_molecule,
            highlight_atoms=highlight_atoms,
            ignore_sdc=ignore_sdc,
            ignore_sdx=ignore_sdx,
            ignore_sa=ignore_sa,
            sdc_threshold=sdc_threshold,
            sdx_threshold=sdx_threshold,
            sa_threshold=sa_threshold,
        )
        vis_vector["svg"] = svg
        vis_vector["smiles"] = config.smiles
        vis_vector["__status"] = "success"
    except JazzyError:
        vis_vector["__status"] = "error"
        vis_vector["smiles"] = config.smiles
    print(vis_vector)
