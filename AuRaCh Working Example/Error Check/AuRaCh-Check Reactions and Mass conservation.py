# -*- coding: utf-8 -*-
"""
This is the supplementary code to the publication:

Fritsch et al., Radiolysis‐Driven Evolution of Gold Nanostructures –
 Model Verification by Scale Bridging In Situ Liquid‐Phase Transmission Electron
 Microscopy and X‐Ray Diffraction, Advanced Science, 2022, DOI:10.1002/advs.202202803

If you find this tool helpful, please cite our work.

This code is published under a GPL-3.0 licence.
"""

import re
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


def extract_prefix(string, re_object):
    """
    Function to extract the prefix of a reactant based on regular expressions

    Parameters
    ----------
    string : str
        A string containing the reactant and an optional prefix in integer syntax.
    re_object : re.Pattern object
        

    Returns
    -------
    prefix : int
        
    reactant : str

    """
    #match regular expression to string
    re_result = re_object.search(string)
    #if re_result is not None, this will be truthy
    if re_result:
        #split prefix and reactant
        try:
            prefix, reactant = string.split()
        except ValueError:
            print(string)
            raise ValueError
        #convert prefix to int
        prefix = int(prefix)
    
    else:
        #dummy prefix for later calculation
        prefix = 1
        #remove spaces
        reactant = string.split()[0]
    
    return prefix, reactant



def remove_pattern(reactant, re_object):
    """
    Removes a predefined re-pattern from a string
    
    INPUT:
        reactant (str): string to be operated on
        re_object (re.Pattern): re.Pattern object with defined pattern
        
    RETURN:
        reactant with pattern removed
    """
    
    parts = re_object.findall(reactant)
    
    for p in parts:
        #remove parenthesis part:
        reactant = ''.join(reactant.split(p))
    
    return reactant



def get_charge(reactant):
    """
    Extract charge from reactant
    
    INPUT:
        reactant(str)
        
    RETURN:
        num (int): Charge
    """
    reactant = '^'.join(reactant.split('**'))
    
    find_negative = re.compile('(\^\d*-|\-)')
    find_positive = re.compile('(\^\d*\+|\+)')
    find_num = re.compile('\d+')
            ##count charges
    num = 0
    #negative:
    charge_match = find_negative.findall(reactant)
    if charge_match:
        if charge_match[-1] == '-':
            num = -1
        elif '-' in charge_match[-1]:
            #extract number
            try:
                num = find_num.search(charge_match[-1]).group()
            except AttributeError:
                num = 1
            num = -int(num)
    
    #positive:
    charge_match = find_positive.findall(reactant)
    if charge_match:
        
        assert num == 0, f'Positive and negative charges found in {reactant}'
        
        if charge_match[-1] == '+':
            num = 1
            
        elif '+' in charge_match[-1]:
            #extract number
            try:
                num = find_num.search(charge_match[-1]).group()
            except AttributeError:
                num = 1
            num = int(num)            
    

    return num
    


def analyze_reaction_site(site):
    """
    Takes educt/product site of a reaction and extracts reactants and the net charge
    """
    #prepare prefix extraction
    find_prefix = re.compile('^\s?\d+\.?\d*\s+')
    #check educts
    elements = dict()
    
    #prepare charge finding
    charges = 0
    
    for agent in site:
        #separate prefix
        prefix, reactant = extract_prefix(agent, find_prefix)
        
        #get atoms in reactant
        elems = extract_atoms([reactant])
        
        #combine atoms and reactant-prefix
        for atom in elems:
            
            elems[atom] *= prefix
            
            try:
                elements[atom] += elems[atom]
            except KeyError:
                elements[atom] = elems[atom]
                
        charge = get_charge(reactant)
        
        charges += prefix * charge
        
    return elements, charges
    
    
    
def extract_atoms(reactant_list):
    """
    Counts the atoms for reactants in a list of reactants.
    """
    #prepare atom extraction
    find_atoms = re.compile('[A-Z][a-z]*\d*')
    #prepare parenthesis extraction:
    find_parenthesis = re.compile('\(.*\)')
    #prepare number extraction
    find_num = re.compile('\d+')
    
    #store elements
    elements = dict()
    
    #iterate over reactants:
    for reactant in reactant_list:
  
        #remove parenthesis
        reactant = remove_pattern(reactant, find_parenthesis)
        
        #get elements
        for atom in find_atoms.findall(reactant):
            
            #get number of a single atom:
            try:
                num = int(find_num.findall(atom)[0])
            except IndexError:
                num = 1
            
            #strip number
            raw_atom = remove_pattern(atom, find_num)
            
            try:
                elements[raw_atom] += num
            except KeyError:
                elements[raw_atom] = num
        
        
    return elements



def extract_species(reaction_string):
    """
    Extract individual reactants out of a reaction string
    """
    #replace tabs by spaces:
    reaction_string = ' '.join(reaction_string.split('\t'))
    #split educts and products into separate lists
    educts, products = reaction_string.split('-->')
    educts = educts.split(' + ')
    products = products.split(' + ')
    
    #remove all spaces at end or beginning of string
    educts = [e.strip() for e in educts]
    products = [p.strip() for p in products]
    
    return educts, products



def analyze_reaction(reaction_string):
    """
    Get element and charge count for educt and product site of a reaction string
    """
    educts, products = extract_species(reaction_string)
    
    educt_elements, educt_charge = analyze_reaction_site(educts)
    product_elements, product_charge = analyze_reaction_site(products)
    
    return educt_elements, educt_charge, product_elements, product_charge


def _extract_agent(agents, find_prefix = re.compile('^\s?\d+\.?\d*\s+')):
    
    return [extract_prefix(a, find_prefix)[1] for a in agents] 



def _extract_reactants(reaction):
    
        educts, products = extract_species(reaction)
        educts = _extract_agent(educts)
        products = _extract_agent(products)
        
        return educts, products
    


def find_identical_reactions(reaction_list):
    """
    Searches for pairs of reactions with matching educt and product species.
    Hence, if a reaction appears more than twice, the number of matches scales
    accordingly.
    """
    pairs = []
    #iterate over list of strings skipping the last one
    for n, string in enumerate(reaction_list[:-1]):
        
        educts, products = _extract_reactants(string)       
        
        #iterate over residual strings not checked yet, starting with n+1 to avoid diagonals (n=n)
        for m, string2 in enumerate(reaction_list[n+1:], n+1):
            educts2, products2 = _extract_reactants(string2)
        
            #now, check for equality using sets:
            if all([set(educts) == set(educts2),
                    set(products) == set(products2)]):
                #add indices starting from 1 to be in accordance with DataFrame
                pairs.append([n+1,m+1])
                
    return pairs
                

def check_reactions(reaction_list):
    """
    Checks a list of reactions for consistency.
    """
    doc = {}
    
    for n, reaction in enumerate(reaction_list, 1):
        #extract elements and charges for product/educt site
        try:
            educt_e, educt_c, product_e, product_c = analyze_reaction(reaction)
        except Exception as e:
            print(f'Error in Reaction number {n+1}.\n\t{reaction}')
            raise e
        #check charges
        charge_test = educt_c == product_c
        
        #check elements
        elem_test = educt_e == product_e
        
        doc[n] = [charge_test, elem_test]
        
        #sort dictionarys for better readability
        educt_e = {k:educt_e[k] for k in sorted(educt_e)}
        product_e = {k:product_e[k] for k in sorted(product_e)}

        
        if not all(doc[n]):
            
            message = f'Error in Reaction {n}:\n'
            
            if not charge_test:
                message += f'\tEduct charge = {educt_c}\n'
                message += f'\tProduct charge = {product_c}\n'
        
            if not elem_test:
                message += f'\tEduct elements = {educt_e}\n'
                message += f'\tProduct elements = {product_e}\n'
                
            message += '\t' + reaction
            
            print(message)
            
        #Charge test, Element test, Educt charges, Product charges, Educt Elements, Product Elements, Reaction
        doc[n].extend([educt_c,  product_c, product_e, product_e, reaction])
        
    return doc



def read_reaction_file(file_str):
    """
    Reads a reaction file, checkes for elementary and charging consistency,
    stores this result as csv file and returns the corresponding output as
    DataFrame.
    """
    with open(file_str) as file:
        reaction_list = file.readlines()
    
    #strip other columns:
    reaction_list = [line.split(';')[0] for line in reaction_list]
    
    out = check_reactions(reaction_list)
    
    #store data in DataFrame
    df = pd.DataFrame.from_dict(out, orient='index',
                                columns=['Charge test',
                                         'Element test',
                                         'Educt charges',
                                         'Product charges',
                                         'Educt Elements',
                                         'Product Elements',
                                         'Reaction'])
    
    df['Reaction'] = [' '.join(r.split('\t')) for r in df['Reaction']]
    df.to_csv(file_str.split('.')[0]+'_checked.csv')
    
    return df



def check_concentration(file):
    """
    Plots the evolution of elements of a simulation output file and returns the DataFrame.
    """
    df = pd.read_csv(file)
    #set time column as index
    df.set_index('time / s', inplace=True)

    #remove unit if necessary:
    species = [s.split(' ')[0] for s in df.columns] 
    
    #elements as list:
    elements = [*extract_atoms(species)]
    
    #do fancy numpy stuff to get a quick result
    out = np.zeros((df.shape[0], len(elements)))
    
    #iterate over species
    for s in species:

        #get dictionary with atoms in s
        atoms = extract_atoms([s])
        #iterate over each atom
        for a in atoms:
            
            #get position in element list
            idx = elements.index(a)
            #add concentration times the atom for each column
            out[:,idx] += df[s] * atoms[a]
    
    #summarize in DataFrame to store column name
    df_out = pd.DataFrame(out, columns=elements)
    #match index with simulation
    df_out['time / s'] = df.index
    df_out.set_index('time / s', inplace=True)
    
    #show & store results:
    
    savename = ''.join(file.split('.')[:-1])
    
    f, ax = plt.subplots(figsize=(5,5), dpi=100)
    for s in df_out:
        
        ax.loglog(df_out[s], label=s)
    
    ax.tick_params(direction='in', which = 'both')
    ax.legend(loc=0)
    ax.set(xlabel= 'time / s' , 
           ylabel='Concentration / M',
           xlim=(1e-12, df_out.index.max())
           )
    plt.show()
    
    f.savefig(savename+'_checked.png', dpi=600, bbox_inches='tight')
    plt.close('all')
    
    df_out.to_csv(savename+'_checked.csv')
    
    return df, df_out



def display_identical_reactions(filename):
    """
    INPUT: (path or str) Name of the reaction file to check, being passed to read_reaction_file
    
    RETURN:
        duplicates: (pd.DataFrame) Contains 2 columns: 'Reaction', 'Match'.
    Its index matches with rows in the reaction file.
    'Reaction' displays the respective reaction string.
    'Match' displays the first match found. Consequently, if a reaction is duplicated
    more than once, it will be found more often, until all combinations are shown.
    Permutations are avoided.
    """
    #read file
    out=read_reaction_file(filename)
    #extract pairs of doubled reactions
    pairs = find_identical_reactions(out['Reaction'].values)
    #define containers for later indexing
    partner_idx = []
    indexes = []
    #flatten pairs tuples and append to containers in respective order
    for pair in pairs:
        idx, idx2 = pair
        indexes.append(idx)
        partner_idx.append(idx2)
        indexes.append(idx2)
        partner_idx.append(idx)
    
    #select relevant rows out of out-DataFrame
    duplicates = out.loc[indexes]    
    #assign respective partner indices
    duplicates['Match'] = partner_idx
    #return relevant columns
    return duplicates[['Reaction', 'Match']]
#%%

filename = 'Reactions_HAuCl4_errors.txt'
duplicates = display_identical_reactions(filename)
print(duplicates)

#%%
conc_file = 'Reactions_HAuCl4_C0_HAuCl4_1000.0_enm2s_Output.csv'
df, df_out = check_concentration(conc_file)

