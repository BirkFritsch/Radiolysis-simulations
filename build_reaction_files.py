# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 17:28:05 2023

@author: b.fritsch
"""
import numpy as np
import pandas as pd
import AuRaCh as au


def _label_used_reactions(used_reactions_list, number_of_all_reactions, file_name):
    """
    Produces a pd.DataFrame labelling all reactions, when producing a sparse reaction set.
    The DataFrame will be stored as Excel (.xlsx) file in the working directory.

    Parameters
    ----------
    used_reactions_list : list, or iterable
         container (list or similar) with a boolean entry for each reaction in .
    number_of_all_reactions : int
        Number of reactions in the original reaction file.
    file_name : str or path-like
        String for storing the output DataFrame as excel. It will automatically appended by '_legend.xlsx'.

    Returns
    -------
    df : pd.DataFrame
        DataFrame containing a column named 'Used' of boolean entries and an index named 
        'Reaction Line Number'. 'Reaction Line Number' starts at 1.
    """
    reaction_label = [True if i in used_reactions_list else False 
                      for i in range(number_of_all_reactions)]
    #store used reactions as excel file       
    df = pd.DataFrame(reaction_label, columns=['Used'])
    df['Reaction Line Number'] = np.arange(1, len(reaction_label) + 1)
    df.set_index('Reaction Line Number', inplace=True)
    df.to_excel(file_name.split('.txt')[0] + '_legend.xlsx')
    
    return df


def write_reaction_file(file_name, new_file_list):
    """
    Writes a text file named file_name with every line being a string in new_file_list.

    Parameters
    ----------
    file_name : str
        File name including the chosen ending
    new_file_list : list or iterable
        List or similar containing strings.

    Returns
    -------
    None.
    """
    with open(file_name, 'w') as file:
        for reaction in new_file_list:
            file.write(f'{reaction}')


def sparsen_reaction_set(reaction_file, excluded_species):
    """
    Creates a new reaction file being a subset of reaction_file  that excludes all reactions using any of the
    reactants listed in excluded_species.

    Parameters
    ----------
    reaction_file : str or path-like
        AuRaCh-ready filename of the original reaction file.
    excluded_species : list
        list or iterable containing species names as strings matching entries 
        in reactions of reaction_file.

    Returns
    -------
    df : pd.DataFrame
        DataFrame containing a column named 'Used' of boolean entries and an index named 
        'Reaction Line Number'. 'Reaction Line Number' starts at 1. This is meant as a legend
        to track which reactions of the original file were used.
    """
    new_file = []
    used_reactions = []
    
    with open(reaction_file) as file:
        
        for n, line in enumerate(file.readlines()):
            
            #crop comments
            reaction_column = line.split(';')[0]
            
            line_fragments = reaction_column.split()
            if any([any([sp == s for s in line_fragments]) for sp in excluded_species]): 
                continue
        
            new_file.append(line)
            used_reactions.append(n)
        
    #write new file
    #new filename
    new_file_name = reaction_file.split('.txt')[0] + '_sparsened.txt'
    
    write_reaction_file(new_file_name, new_file)
            
    df = _label_used_reactions(used_reactions, n, new_file_name)
    
    return df


def build_reaction_subset(reaction_file, included_species):
    """
    Creates a new reaction file being a subset of reaction_file that includes all 
    reactions using any of the reactants listed in included_species.


    Parameters
    ----------
    reaction_file : str or path-like
        AuRaCh-ready filename of the original reaction file.
    included_species :  list
        list or iterable containing species names as strings matching entries 
        in reactions of reaction_file.
    Returns
    -------
    df : pd.DataFrame
        DataFrame containing a column named 'Used' of boolean entries and an index named 
        'Reaction Line Number'. 'Reaction Line Number' starts at 1. This is meant as a legend
        to track which reactions of the original file were used.
    """
    new_file = []
    used_reactions = []
    with open(reaction_file) as file:
        
        for n, line in enumerate(file.readlines()):
            #check if line is reaction to exclude header
            if '-->' not in line:
                continue
            
            #crop to reactions
            reaction_column = line.split(';')[0]
            line_fragments = reaction_column.split()
            #remove functional symbols
            reaction_operators = ['+', '-->']
            line_fragments = [s for s in line_fragments if not s.isdigit()]
            reactants = [s for s in line_fragments 
                              if not any([s.strip() == o for o in reaction_operators])
                              ]
            
            if all([any([sp == incl_sp for incl_sp in included_species]) 
                    for sp in reactants]): 
                
                new_file.append(line)
                used_reactions.append(n)
    
    name_tag = " ".join(included_species)
    if len(included_species) > 20:
        name_tag = " ".join(included_species[:10]) + ' and more'
    new_file_name = f'Reactions of {name_tag}.txt'
    
    write_reaction_file(new_file_name, new_file)
    
    df = _label_used_reactions(used_reactions, n, new_file_name)
    
    return df


def automated_reaction_detection(reaction_file, concentrations_gval_file):
    """
    Detects all participating reactions in reaction_file based on simulation starting conditions.
    This means, it finds all reactions based on all species that have non-zero concentrations or 
    G-values in the AuRaCh-ready concentrations_gval_file.
    The function then iterates over all reactions in reaction_file to identify all secondary
    products to subsequently build a reaction set that contains all reactions of reaction_file
    where all initial and secondary reactants are regarded.


    Parameters
    ----------
    reaction_file : str or path-like
        AuRaCh-ready filename of the original reaction file.
    concentrations_gval_file : str or path-like
        AuRaCh-ready filename of the concentrations/g-value file.

    Returns
    -------
    relevant_reactants : list
        list containing all identified reactants in alphabetical order.

    """
    __, rmap = au.create_reactions(reaction_file)
    reactants = au._make_reactant_set(rmap)
    # link initial concentrations to reactants and assign 0 if no other value is given
    dict_reactants, dict_gvals, *__ = au.get_initial_conditions(
        reactants, conditions=concentrations_gval_file
    )
    #retreive all species with non-zero concentrations or g-values
    nonzero_g_vals = [r for r in reactants if dict_gvals[r] != 0]
    nonzero_concentrations = [r for r in reactants if dict_reactants[r] != 0]
    used_reactants = set()
    used_reactants.update(nonzero_g_vals)
    used_reactants.update(nonzero_concentrations)
    
    """now iterate over all educts and look at all products. Add products to used_starting_reactants 
    if it is not in there already"""
    
    initial_number_of_species = 0
    new_number_of_species = len(used_reactants)
    while initial_number_of_species < new_number_of_species:
        for n in rmap:
            if all(s in used_reactants for s in rmap[n]['educts']):
                for product in rmap[n]['products']:
                    used_reactants.add(product)
        initial_number_of_species = new_number_of_species
        new_number_of_species= len(used_reactants)
        
    relevant_reactants = sorted(used_reactants)
        
    return relevant_reactants
        

def build_reaction_set_from_initial_conditions(reaction_file, concentrations_gval_file):
    """
    Creates a new reaction file being a subset of an AuRaCh-ready, data base-like reaction_file
    that includes all reactions based on simulation starting conditions.
    This means, the new reaction file includes all reactions based on all species that have a 
    non-zero concentrations or G-values in the AuRaCh-ready concentrations_gval_file.
    
    Moreover, an Excel (.xlsx) sheet will be stored that references whether a reaction 
    from the original reaction_file was used.

    Parameters
    ----------
    reaction_file : str or path-like
        AuRaCh-ready filename of the original reaction file.
    concentrations_gval_file : str or path-like
        AuRaCh-ready filename of the concentrations/g-value file.
    Returns
    -------
    included_species :  list
        list containing all identified reactants in alphabetical order.
    label_df :  pd.DataFrame
        DataFrame containing a column named 'Used' of boolean entries and an index named 
        'Reaction Line Number'. 'Reaction Line Number' starts at 1. This is meant as a legend
        to track which reactions of the original reaction_file were used.
    """
    included_species = automated_reaction_detection(reaction_file, concentrations_gval_file)
    label_df = build_reaction_subset(reaction_file, included_species)
    
    return included_species, label_df


if __name__ == '__main__':
    reaction_file = 'Reaction_Database.txt'
    gval_file = 'conditions_H2O.txt'
    
    #first example
    excluded_species = 'Cl2 Cl2O Cl2O2 Cl2O3 Cl2O4 Cl3- ClO ClO- ClO2 ClO2- ClO3 ClO3- HClO HClO2 HO3 O O- O3 O3- O4'.split()
    df = sparsen_reaction_set(reaction_file, excluded_species)
    
    #second example
    included_species = ['H2O', 'H+', 'OH-', 'OH', 'N']
    df = build_reaction_subset(reaction_file, included_species)

    #third example
    used_reactants, df = build_reaction_set_from_initial_conditions(reaction_file, 
                                                                    gval_file)
    
