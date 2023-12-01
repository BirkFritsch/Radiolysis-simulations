# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 17:28:05 2023

@author: Birk Fritsch

This script is based on

Fritsch et al., "Precision of Radiation Chemistry Networks: Playing Jenga with Kinetic Models for Liquid-Phase Electron Microscopy",
Precision Chemistry, 2023, doi:10.1021/prechem.3c00078

and is to be regarded as addition to AuRaCh - Automated Radiation Chemistry simulations.

Dependencies:
    
Python 3.10.6+ (requires f-string compatibility, written and tested with 3.10.6, should also work with 3.6+ (not tested))
pandas 1.5.1+ (may work with other versions, but written and tested with 1.5.1)
"""
import pandas as pd


def sparsen_reaction_set(reaction_file, excluded_species):
    """
    Function to write a new reaction txt file without any reactions containing excluded species.
    The file is stored is current working directory as "<reaction_file>_sparsened.txt".
    In addition, a csv table listing the number of all initial reactions with a flag whether they are kept or
    thrown out is stored as "<reaction_file>_legend.csv". The latter table is also returned as
    pandas DataFrame.

    Parameters
    ----------
    reaction_file : str
        Name of the to-be-sparsened reaction file.
    excluded_species : list
        A list containing reactants that should be excluded from reaction_file.

    Returns
    -------
    df : pandas.DataFrame
        Index: Number (int) of reaction in reaction_file, starting from 1. The index is named "Reaction".
        Columns: Single column named "Used" containing boolean values.

    """
    new_file = []
    used_reactions = []

    with open(reaction_file) as file:

        for n, line in enumerate(file.readlines()):
            # crop comments
            reaction_column = line.split(";")[0]
            line_fragments = reaction_column.split()
            if any([any([sp == s for s in line_fragments]) for sp in excluded_species]):
                continue

            new_file.append(line)
            used_reactions.append(n)

    # write new file
    # new filename
    new_file_name = reaction_file.split(".txt")[0] + "_sparsened.txt"

    with open(new_file_name, "w") as file:

        for reaction in new_file:
            file.write(f"{reaction}")

    # store used reactions as csv file
    reaction_label = [True if i in used_reactions else False for i in range(n+1)]
    df = pd.DataFrame(reaction_label, columns=["Used"], index=range(1, len(reaction_label)+1))
    df.index.name="Reaction"
    df.to_csv(new_file_name.split(".txt")[0] + "_legend.csv")

    return df


if __name__ == "__main__":
    reaction_file = "Reactions_HAuCl4_extended.txt"
    excluded_species = (
        "Cl2 Cl2O Cl2O2 Cl2O3 Cl2O4 Cl3- ClO ClO- ClO2 ClO2- ClO3 ClO3- HClO HClO2 HO3 O O- O3 O3- O4".split()
    )
    df = sparsen_reaction_set(reaction_file, excluded_species)
