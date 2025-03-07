set# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 11:12:12 2022

@author: b.fritsch, Guiseppe De Salvo
"""

import re
import numpy as np
import pandas as pd
from .. import AuRaCh as au


def ensure_integer_prefactor(prefactor):
    """
    This is required because COMSOL does only work well with integers as
    stoichiomentric factors

    Parameters
    ----------
    prefactor : TYPE
        DESCRIPTION.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    None.

    """
    if prefactor != int(prefactor):
        
        raise ValueError(f'{prefactor} is not convertible to an integer value')


def calculate_conversion(dct):

    #reaction order
    n = dct['educt prefix'].sum()
    
    #conversion M --> mM
    conversion = 1000**(-n+1)
    
    return conversion


def replace_charge_with_letter(species):
    
    return 'P'.join('n'.join(species.split('-')).split('+'))


def make_species_4_comsol(species_pn):
    
    
    #replace +/- signs with P/n. Capitalization of P is chosen to preserve alphanumeric ordering
    #species_pn = 'P'.join('n'.join(species.split('-')).split('+'))
    #replace ^ with 
    if '^' in species_pn:
        body, charge = species_pn.split('^')
        #get charges
        searcher = re.search('\d+', charge)
        new_charge = charge
        if searcher:
            factor_str = searcher.group()
            _, charge_letter = charge.split(factor_str)
            factor = int(factor_str)    
            new_charge = factor*charge_letter 
        
        species_pn = body + new_charge
        
    for sign in '[]()':    
        if sign in species_pn:
            species_pn = species_pn.replace(sign, '')
        
    return species_pn
            

def write_reaction_file_2_comsol(reaction_file):
     
    rate_constants, reaction_dictionary = au.create_reactions(reaction_file)
    
    output_string = 'Reaction number, formula, rate constant\n'
    
    for reaction in reaction_dictionary:
        
        
        #get reaction number written
        r_string = f'{reaction+1}, '
        
        #write educts:
        educt_str_lst = []
        for educt, prefactor in zip(reaction_dictionary[reaction]['educts'],
                                    reaction_dictionary[reaction]['educt prefix']):
            
            ensure_integer_prefactor(prefactor)
            comsol_educt = make_species_4_comsol(educt)
            educt_str_lst.append(f'{int(prefactor)}{comsol_educt}')
        educt_str = '+'.join(educt_str_lst)
            
        #add educt part to r_string
        r_string += educt_str
        
        #add reaction file
        r_string += '=>'
        
        #write product part
        product_str_lst = []
        for product, prefactor in zip(reaction_dictionary[reaction]['products'],
                                    reaction_dictionary[reaction]['product prefix']):
            
            ensure_integer_prefactor(prefactor)
            comsol_product = make_species_4_comsol(product)
            product_str_lst.append(f'{int(prefactor)}{comsol_product}')
        product_str = '+'.join(product_str_lst)    
                    
        #add product part to r_string
        r_string += product_str
        
        output_string += r_string
        output_string += ', '
        
        #add rate constant
        rconst = rate_constants[reaction]
        #convert to mM
        rconst *= calculate_conversion(reaction_dictionary[reaction])
        
        output_string += f'{rconst}\n'
    
    filename = reaction_file.split('.txt')[0] + '_COMSOL.txt'
    
    with open(filename, 'w') as textfile:
        textfile.write(output_string)


def write_concentration_file_2_comsol(reaction_file,
                                      concentration_file):

    rate_constants, reaction_dictionary = au.create_reactions(reaction_file)
    
    #create COMSOL-readable reaction file    
    reactants = au._make_reactant_set(reaction_dictionary)   
    
    dict_reactants, dict_gvals = au.get_initial_concentrations(
        reactants, concentrations=concentration_file
    )
    
    c0_list = []
    g_val_list = []
    reactants_COMSOL = []
    reactants_reaction_set = []
    
    for species in reactants:
        
        c0 = dict_reactants[species]
        c0_list.append(f'{c0}[M]')
        
        g_val = dict_gvals[species]
        g_val_list.append(f'{g_val}[1/eV]')
        
        species_reaction_set = make_species_4_comsol(species)
        #replace +/- signs with P/n. Capitalization of P is chosen to preserve alphanumeric ordering
        species_pn = replace_charge_with_letter(species_reaction_set)
        reactants_COMSOL.append(species_pn)
        #keep species as written in the reaction file for sorting 
        reactants_reaction_set.append(species_reaction_set)
    
    df=pd.DataFrame(np.array([c0_list,
                              g_val_list,
                              reactants_COMSOL,
                              reactants_reaction_set]).T,
                    columns=['InitialConcentration',
                             'g-Value',
                             'Species',
                             'Sorting column'])
    df.sort_values(by=['Sorting column'],
                   inplace=True)
    df.drop(columns='Sorting column',
            inplace=True)
    df.set_index('Species',
                 inplace=True)
    
    filename = concentration_file.split('.txt')[0] + '_COMSOL.txt'
    df.to_csv(filename)
    

def main(reaction_file,
         concentration_file):
    
    write_reaction_file_2_comsol(reaction_file)
    write_concentration_file_2_comsol(reaction_file,
                                      concentration_file)
    
    


#%%

if __name__ == "__main__":
    
    reaction_file = 'Reactions_HAuCl4_sparsened.txt'
    concentration_file = 'C0_HAuCl4.txt'
    rate_constants, reaction_dictionary = au.create_reactions(reaction_file)
    main(reaction_file,
         concentration_file)
        