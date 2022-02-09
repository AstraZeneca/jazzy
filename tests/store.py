"""Global molecular structure store for unit tests."""
# tests/store.py
from kallisto.atom import Atom
from kallisto.molecule import Molecule
from kallisto.units import Bohr


def pyridine():
    """Create a pyridine molecule and return as Molecule."""
    coords = [
        [1.3603, 0.0256, 0.0],
        [0.6971, -1.202, 0.0],
        [-0.6944, -1.2184, 0.0],
        [-1.3895, -0.0129, 0.0],
        [-0.6712, 1.1834, 0.0],
        [0.6816, 1.196, 0.0],
        [2.453, 0.1083, 0.0],
        [1.2665, -2.1365, 0.0],
        [-1.2365, -2.1696, 0.0],
        [-2.4837, 0.0011, 0.0],
        [-1.1569, 2.1657, 0.0],
    ]
    coords = [[float(j) / Bohr for j in i] for i in coords]

    symbols = [
        "C",
        "C",
        "C",
        "C",
        "C",
        "N",
        "H",
        "H",
        "H",
        "H",
        "H",
    ]

    atoms = []
    for i, _ in enumerate(coords):
        atoms.append(Atom(symbols[i], position=coords[i]))
    return Molecule(symbols=atoms)


def neopentane():
    """Create a neopentane molecule and return as Molecule."""
    coords = [
        [0.000000, 0.0, 0.0],
        [0.881905, 0.881905, 0.881905],
        [-0.881905, -0.881905, 0.881905],
        [0.881905, -0.881905, -0.881905],
        [-0.881905, 0.881905, -0.881905],
        [-1.524077, 0.276170, -1.524077],
        [1.524077, 1.524077, 0.276170],
        [1.524077, -0.276170, -1.524077],
        [1.524077, 0.276170, 1.524077],
        [-1.524077, -0.276170, 1.524077],
        [1.524077, -1.524077, -0.276170],
        [-0.276170, 1.524077, -1.524077],
        [0.276170, 1.524077, 1.524077],
        [0.276170, -1.524077, -1.524077],
        [-0.276170, -1.524077, 1.524077],
        [-1.524077, 1.524077, -0.276170],
        [-1.524077, -1.524077, 0.276170],
    ]
    coords = [[float(j) / Bohr for j in i] for i in coords]

    symbols = [
        "C",
        "C",
        "C",
        "C",
        "C",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
    ]

    atoms = []
    for i, _ in enumerate(coords):
        atoms.append(Atom(symbols[i], position=coords[i]))
    return Molecule(symbols=atoms)
