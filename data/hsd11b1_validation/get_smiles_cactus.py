"""Converts synonyms into SMILES for the data from Gerber's paper."""
# data/hsd11b1_validation/get_smiles_cactus.py
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
    data = "01-robb_data.txt"
    df = pd.read_csv(data, sep="\t")
    df["SMILES"] = df.apply(lambda row: getsmiles_cactus(row["Iupac"]), axis=1)
    df.to_csv("02-robb_data_smiles.txt", sep="\t")


main()
