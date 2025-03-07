# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 11:47:09 2025

@author: b.fritsch
"""

import re
import numpy as np
import pandas as pd
from decimal import Decimal
import AuRaCh as au

def make_nice_exp(s):
    
    if str(s).strip() == '':
        return ''
    
    elif s == int('0'):
        return '$0$'
    try:    
        s = f'{Decimal(s):.3e}'
    except Exception as exception:
        print(s)
        raise exception
    num, exp = s.split('e')
    
    if float(exp) == 0:
                
        return f'${int(float(num))}$'
        
    for i in range(len(num)-1,0, -1):
        if i < 1:
            break
        if num[i] == '0':
            num = num[:i]
    
    if num.endswith('.'):
        num = num[:-1]
        
    if exp.startswith('+'):
        exp = exp.split('+')[1]

    return f'${num}\cdot10^' + '{' + str(exp) + '}$'


def make_nice_string(raw):

    store = []
    for agent in raw:
        s_out = ''
        prefix, reactant = au._extract_prefix(agent, re.compile('^\s?\d+\.?\d*\s+'))

        if prefix != 1:
            s_out += f'{prefix} '

        s_out += au.make_nice_label(reactant)#'\\cf{' + reactant + '}'

        store.append(s_out)

    return store


def make_nice_equation(s, only_half='educt'):
    
    s = ' '.join(s.split('\t'))
      
    #adjust reactants:
    educt_site, product_site = s.split('-->')
    #split product and educt sites into individual reactants
    raw_educts = educt_site.split(' + ')
    raw_products = product_site.split(' + ')

    
    #introduce mhchem syntax   
    educt_list = make_nice_string(raw_educts)
    product_list = make_nice_string(raw_products)
    
    if only_half:
    #reaction arrow
        if only_half=='educt':
            s_designed = ' + '.join(educt_list)
        elif only_half == 'product':
            s_designed = ' + '.join(product_list)
        else:
            raise ValueError(f'Unknown only_half value: {only_half}')
    else:
        s_designed = ' + '.join(educt_list) + ' $\\longrightarrow$ ' + ' + '.join(product_list)
    
    #remove additional spaces
    s_designed = ' '.join(s_designed.split())
    
    return s_designed    
    

def convert_reaction_set_to_LaTeX(reaction_file, clear_comments=False):
    """
    Using pandas to export a reaction_file to be implemented in LaTeX documents 
    via booktabs and longtable.
    This requires 
    
    \\usepackage{booktabs}
    
    \\usepackage{longtable}
    
    in your LaTeX preamble.
    
    See https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_latex.html
    
    
    
    Parameters
    ----------
    reaction_file : str
        csv file containing reaction sets for AuRaCh
    clear_comments : bool, optional
        Flag whether to clear the third column from the right. The default is False.

    Returns
    -------
    data : pd.DataFrame
        DataFrame built from reaction_file containing LaTeX notations.
    """
    
    #load data
    data = pd.read_csv(reaction_file, sep=';', header=0)
    if data.shape[1] == 3:
        data = pd.read_csv(reaction_file, sep=';', names=['Reaction', 'Rate Constant', 'Comment'])
    elif data.shape[1] == 2:
        data = pd.read_csv(reaction_file, sep=';', names=['Reaction', 'Rate Constant'])
    column_mapper = {col:col.strip().title() for col in data.columns}
    data.rename(columns=column_mapper, inplace=True)
    data.set_index(np.arange(1,data.shape[0]+1), inplace=True)
    
    #convert strings
    data['$~$'] = data['Reaction'].apply(lambda s:make_nice_equation(s, only_half='product'))
    data['Reaction'] = data['Reaction'].apply(lambda s:make_nice_equation(s, only_half='educt'))
    data['\\,'] = '$\\longrightarrow$'
    for col in data.columns:
        if col in ['$~$', '\\,', 'Reaction', 'Comment']:
            continue
        data[col] = data[col].apply(make_nice_exp)
    
    #restructure column order
    final_column_order = ['Reaction', '\\,', '$~$'] #
    for col in data.columns:
        if (col in final_column_order) or (col == 'Comment'):
            continue
        if data[col].any():
            final_column_order.append(col)
    #comment column should go last
    if 'Comment' in data.columns:
        final_column_order.append('Comment')    
    
    data = data[final_column_order]
    
    if clear_comments:
        data['Comment'] = ''
    
    #export using pandas
    save_name = reaction_file
    if '.' in save_name:
        save_name = '.'.join(save_name.split('.')[:-1])
    save_name += '.tex'
    
    #predefine first four columns (index + reaction)
    column_formatter = 'rrcl'
    #now, put all column formats to 'r' for every non-comment column, which goes last
    idx_of_preset_columns = len(column_formatter) - 1
    for col in range(len(data.columns[idx_of_preset_columns:-1])):
        column_formatter +='r'
    #make comment column leftbound
    column_formatter += 'l'
    
    with pd.option_context("max_colwidth", 100000):
        data.to_latex(save_name, longtable=True, escape=False, column_format=column_formatter)

    return data


if __name__ == '__main__':
    reaction_file = 'CO2_upconversion_on_Cu.txt'
    data = convert_reaction_set_to_LaTeX(reaction_file, clear_comments=True)
    
    