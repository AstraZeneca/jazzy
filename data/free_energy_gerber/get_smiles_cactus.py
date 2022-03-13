"""Converts synonyms into SMILES for the data from Gerber's paper."""
# data/free_energy_gerber/get_smiles_cactus.py
import re
from io import BytesIO

import pandas as pd
import pycurl


def getsmiles_cactus(name):
    """Converts synonyms into SMILES strings.

    A function to use the public cactus (National Institutes of Cancer Research)
    webservice to retrieve a smiles string from a synonym.

    Args:
    name: any trivial or IUPAC name for a molecule

    Returns:
    Canonical smiles string for that molecule.

    """
    # hard-coded rules for some unusual names
    name = re.sub(" ", "%20", name)
    name = re.sub("Propine", "Propyne", name)
    name = re.sub("Butine", "Butyne", name)
    name = re.sub("Pentine", "Pentyne", name)
    name = re.sub("Hexine", "Hexyne", name)
    name = re.sub("Heptine", "Heptyne", name)
    name = re.sub("Octine", "Octyne", name)
    name = re.sub("Nonine", "Nonyne", name)
    name = re.sub("1-Buten-3-ine", "1-Buten-3-yne", name)
    name = re.sub("methoate", "methanoate", name)
    url = "https://cactus.nci.nih.gov/chemical/structure/" + name + "/smiles"
    buffer = BytesIO()
    c = pycurl.Curl()
    c.setopt(c.URL, url)
    c.setopt(c.WRITEDATA, buffer)
    c.perform()
    c.close()
    smiles = buffer.getvalue().decode("UTF-8")
    print(name, smiles)
    return smiles


def main():
    """Runs a batch of name conversions into SMILES."""
    data = "01-free_energies_gerber.txt"
    df = pd.read_csv(data, sep="\t")
    df["SMILES"] = df.apply(lambda row: getsmiles_cactus(row["Structure"]), axis=1)
    df.to_csv("02-free_energies_gerber_smiles.txt", sep="\t")


main()
