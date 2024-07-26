"""Model classes for attribute definition."""
# src/jazzy/model.py
from enum import Enum


class ChargeKey(str, Enum):
    """Attributes of the charge dictionary."""

    Atom: str = "q"
    Alpha: str = "q_alpha"
    Beta: str = "q_beta"
    Gamma: str = "q_gamma"


class IdxKey(str, Enum):
    """Attributes of the atom index dictionary."""

    Atom: str = "atom"
    Alpha: str = "alpha"
    Beta: str = "beta"
    Gamma: str = "gamma"
