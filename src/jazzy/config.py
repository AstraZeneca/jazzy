"""Global config file of jazzy package."""
# src/jazzy/config.py

MINIMISATION_METHOD = None
ROUNDING_DIGITS = 4
ANNOTATION_FONT_SCALE = 0.7


class Config(object):
    """Global config object.

    Global jazzy config file for command-line interface (CLI).
    Passable to CLI methods via '@pass_config' decorator.
    Init method initializes the global parameters that are neccessary within
    the calculation of free energy contributions.

    Parameter in config:
    self.g0: apolar zeroth order parameter
    self.gs: apolar surface parameter
    self.gr: apolar ring parameter
    self.gpi1: apolar first pi parameter
    self.gpi2: apolar second pi parameter
    self.gd: hydrogen bond donor parameter
    self.ga: hydrogen bond acceptor parameter
    self.expd: exponent hydrogen bond donor
    self.expa: exponent hydrogen bond acceptor
    self.gi: interaction parameter
    self.smiles: SMILES string of system

    """

    def __init__(
        self,
        g0=1.884,
        gs=0.0467,
        gr=-3.643,
        gpi1=-1.602,
        gpi2=-1.174,
        gd=-0.908,
        ga=-16.131,
        expd=0.50,
        expa=0.34,
        gi=4.9996,
        f=0.514,
    ):
        """Initialize global config."""
        self.g0 = g0
        self.gs = gs
        self.gr = gr
        self.gpi1 = gpi1
        self.gpi2 = gpi2
        self.gd = gd
        self.ga = ga
        self.expd = expd
        self.expa = expa
        self.gi = gi
        self.f = f
        self.smiles = None
