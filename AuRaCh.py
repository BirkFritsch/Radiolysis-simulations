# -*- coding: utf-8 -*-
"""
This is the supplementary code to the publication:

Fritsch et al., Radiolysis‐Driven Evolution of Gold Nanostructures –
 Model Verification by Scale Bridging In Situ Liquid‐Phase Transmission Electron
 Microscopy and X‐Ray Diffraction, Advanced Science, 2022, 9, 2202803, DOI:10.1002/advs.202202803

If you find this tool helpful, please cite our work.

This code is published under a GPL-3.0 licence.
"""

import os
import re
import logging
import traceback
import numpy as np
import pandas as pd
from scipy.constants import N_A, e, gas_constant, h
from scipy.constants import k as kB
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt
from datetime import datetime
from itertools import cycle, product
from decimal import Decimal


logging.basicConfig(level=logging.INFO)

def convert_decimal_to_latex(decimal_str):
    
    prefactor, exponent = decimal_str.split('e')
    
    return f'${prefactor}{{\\cdot}}10^{{{int(exponent)}}}$'


def sci_num_tex(value):
    
    if isinstance(value,str):
        latex_string = convert_decimal_to_latex(value)    
    else:    
        value_str = f"{value:.2e}"
        latex_string = convert_decimal_to_latex(value_str)
        
    return latex_string


def e_flux_to_dose_rate_angstrom(
    flux, thickness, stopping_power=2.36, mean_free_path=3.8e-7, t=1
):
    """
    INPUT:
    
    flux (float or array-like) in e/(Angstroem²s)
    Thickness (float or array-like) https://github.com/BirkFritsch/Radiolysis-simulations/blob/main/AuRaCh%20Working%20Example%20GUI/AuRaCh_GUI.pyin m
    Mean free path (float or array-like) in m
    stopping_power (float or array-like) in MeV(cm)²/g
    
    RETURN:
    
    dose rate (float or array-like) in Gy/s 
    """

    dose_rate = 1e5 * e* (stopping_power/(1e-20*t)) * flux * (1 + thickness/mean_free_path)
    
    return dose_rate


def e_flux_to_dose_rate_nm(flux, thickness, stopping_power = 2.36, mean_free_path=3.8e-7,t=1):
    """
    INPUT:
    
    flux (float or array-like) in e/(nm²s)
    Thickness (float or array-like) in m
    Mean free path (float or array-like) in m
    stopping_power (float or array-like) in MeV(cm)²/g
    
    RETURN:
    
    dose rate (float or array-like) in Gy/s 
    """

    dose_rate = 1e5 * e* (stopping_power/(1e-18*t)) * flux * (1 + thickness/mean_free_path)
    
    return dose_rate



def e_flux_to_dose_rate(flux, thickness, stopping_power = 2.36, mean_free_path=3.8e-7,t=1):
    """
    Wrapper around e_flux_to_dose_rate_angstrom. Deprecated.
    Use  e_flux_to_dose_rate_angstrom or e_flux_to_dose_rate_nm instead.
    
    INPUT:
    

    flux (float or array-like) in e/(Angstroem2s)
    Thickness (float or array-like) in m
    Mean free path (float or array-like) in m
    stopping_power (float or array-like) in MeV(cm)2/g
    
    RETURN:
    
    dose rate (float or array-like) in Gy/s 
    """

    
    return e_flux_to_dose_rate_angstrom(flux, thickness, stopping_power, mean_free_path, t)



def _extract_prefix(string, re_object):
    """
    Function to extract the prefix of a reactant based on regular expressions

    Parameters
    ----------
    string : str
        A string containing the reactant and an optional prefix in float syntax.
    re_object : re.Pattern object


    Returns
    -------
    prefix : float

    reactant : str

    """
    # match regular expression to string
    re_result = re_object.search(string)
    # if re_result is not None, this will be truthy
    if re_result:
        # split prefix and reactant
        try:
            prefix, reactant = string.split()
        except ValueError as v:
            logging.critical(
                f"ValueError for {string}. This means, no patterb are found, although expected by the regular expression."
            )
            raise v
        # convert prefix to float
        prefix = float(prefix)

    else:
        # dummy prefix for later calculation
        prefix = 1.0
        # remove spaces
        reactant = string.split()[0]

    return prefix, reactant


def _create_reactions(filename,
                     settings=None):
    """
    Dechifre a reaction file. Every line in the file specified is interpretated
    as reaction. Educt and product site are separated by '-->'. The rate constant
    is separated from the reaction by a semi colon.
    Example:

        2 A + 4 B --> 2 AB- + 2 B+ ; 0.35

        this would return the following:
            np.array([0.35]), {0:{'educts':['A', 'B'], 'products':['AB-', 'B+'], 'educt prefix':np.array([2.,4.]), 'product prefix':np.array([2.,2.]),}}

    Parameters
    ----------
    filename : str
        Name of the file containing reaction equations and rate constants. The default is 'reactions.txt'.
    settints : dict or None
        settings dictionary or None. Default is None

    Returns
    -------
    rate_constants : numpy.ndarray
        Extracted rate constant in order of appearance.
    reaction_map : dict
        Nested dictionary:
            For every reaction a dictionary with the keys 'educts, products, educt prefix, product prefix' is created.
            educts, products -> lists containing strings
            educt prefix, product prefix -> np.ndarrays containing floats

            The key for this sub_dictionary is an integer based on the order of apperance.

    """
    # read reaction file
    
    logging.info(f'Read {filename}')

    with open(filename) as reaction_file:
        # create a list with every line in file as string
        reaction_list = reaction_file.readlines()

    # create output containers
    rate_constants = []
    reaction_map = {}

    # define re.Pattern object for prefix extraction
    find_prefix = re.compile("^\s?\d+\.?\d*\s+")

    #define adjustable counter
    n = 0
    # do the work, j is for debugging
    for j, reaction in enumerate(reaction_list):

        # replace tabs by spaces:
        reaction = " ".join(reaction.split("\t"))
        

        #Header handling. If no reaction or ; is found, the line is assumed to be a header
        try:
            # split the reaction off 
            reaction_col, *other_cols = reaction.split(";")#
            # divide reaction into two sides based on reaction marker
            educt_site, product_site = reaction_col.split("-->")
            #ensure that other cols are there. This will raise a ValueError if ; is missing
            other_cols[0]
            
        except Exception as e:
            logging.warning(
                f"""Could not handle the reaction file in line {n}: {reaction}.
                    I skip line {n} and continue."""
            )
            logging.info(type(e).__name__)
            logging.info(traceback.format_exc())
            
            continue
        
        # deal with the other columns based on the amount of entries
        num_of_other_cols = len(other_cols)
        
        if num_of_other_cols <= 2:
            
            logging.info(f'Assume given rate constant in reaction {j}')

            #assume old format reaction ; rate constant k; comment
            k = [other_cols[0]]
        
        elif num_of_other_cols == 3:
            
            activation_energy = float(other_cols[1]) 
            arrhenius_constant = float(other_cols[0]) 
            
            k = [calculate_rate_constant(activation_energy,
                                        temperature=settings["temperature"],
                                        a1=arrhenius_constant)]
            
            logging.info(f"line {j}: Arrhenius writing assumed")

        elif num_of_other_cols > 5:
            #calculate rate constant based on other values
            #if first other column is a float, continue:
            try:
                k = [float(other_cols[0])]
            except ValueError:
                #Check for Arrhenius constant
                try:
                    arrhenius_constant = float(other_cols[1])
                    logging.info(f"line {j}: Arrhenius constant detected")
                except ValueError:
                    arrhenius_constant = None
                #read rest of required information
                e0, symmetry_factor, amount_of_charges = (float(other_cols[2]),
                                                          float(other_cols[3]),
                                                          float(other_cols[4]))
                
                activation_energy = potential_dependent_activation_energy(
                                                        e0,
                                                        settings["potential"],
                                                        symmetry_factor,
                                                        amount_of_charges,
                                                        )
                #calculate k based on electrochemistry
                k = [calculate_rate_constant(activation_energy,
                                            temperature=settings["temperature"],
                                            a1=arrhenius_constant)]
        else:
            raise ValueError(f"Unknown amount of columns in line {j}")

        rate_constants.append(float(k[0]))
        # split product and educt sites into individual reactants
        raw_educts = educt_site.split(" + ")
        raw_products = product_site.split(" + ")
        # create output containers for sub dictionary
        educts = []
        products = []
        educt_prefix = np.ones(len(raw_educts))
        product_prefix = np.ones(len(raw_products))
        # get prefixes for every reactant
        for m, educt in enumerate(raw_educts):
            prefix, reactant = _extract_prefix(educt, find_prefix)
            educts.append(reactant)
            educt_prefix[m] = prefix

        for m, product in enumerate(raw_products):
            prefix, reactant = _extract_prefix(product, find_prefix)
            products.append(reactant)
            product_prefix[m] = prefix
        
        # store results
        reaction_map[n] = {
            "educts": educts,
            "products": products,
            "educt prefix": educt_prefix,
            "product prefix": product_prefix,
        }
        #increase n 
        n += 1
    #transform rate_constants to array
    rate_constants = np.array(rate_constants)
    
    return rate_constants, reaction_map


def create_reactions(filename,
                     settings=None):
    """
    Dechifre a reaction file. Every line in the file specified is interpretated
    as reaction. Educt and product site are separated by '-->'. The rate constant
    is separated from the reaction by a semi colon.
    Example:

        2 A + 4 B --> 2 AB- + 2 B+ ; 0.35

        this would return the following:
            np.array([0.35]), {0:{'educts':['A', 'B'], 'products':['AB-', 'B+'], 'educt prefix':np.array([2.,4.]), 'product prefix':np.array([2.,2.]),}}

    Parameters
    ----------
    filename : str
        Name of the file containing reaction equations and rate constants. The default is 'reactions.txt'.
    settints : dict or None
        settings dictionary or None. Default is None

    Returns
    -------
    rate_constants : numpy.ndarray
        Extracted rate constant in order of appearance.
    reaction_map : dict
        Nested dictionary:
            For every reaction a dictionary with the keys 'educts, products, educt prefix, product prefix' is created.
            educts, products -> lists containing strings
            educt prefix, product prefix -> np.ndarrays containing floats

            The key for this sub_dictionary is an integer based on the order of apperance.

    """
    # read reaction file
    
    logging.info(f'Read {filename}')

    with open(filename) as reaction_file:
        # create a list with every line in file as string
        reaction_list = reaction_file.readlines()
    
    #Header handling. If no reaction is found, the line is assumed to be a header
    if '-->' in reaction_list[0]:
        return _create_reactions(filename,
                                 settings)
    
    else:
        reaction_file_df = pd.read_csv(filename, header=0, sep=';\s+')
    # create output containers
        rate_constants = []
        reaction_map = {}
    
        # define re.Pattern object for prefix extraction
        find_prefix = re.compile("^\s?\d+\.?\d*\s+")
    
        #define adjustable counter
        n = 0
        # do the work, j is for debugging
        for j, reaction in enumerate(reaction_file_df["Reaction"]):
    
            #get rate constant:
            k = reaction_file_df.loc[j, 'Rate constant']
            
            if ~np.isnan(k):
                logging.debug(f"line {j}: k is used directly")
            else:
       
                activation_energy = reaction_file_df.loc[j, 'Activation energy']
                arrhenius_constant = reaction_file_df.loc[j, 'Arrhenius constant']
                
                if ~np.isnan(arrhenius_constant):
                
                    k = calculate_rate_constant(activation_energy,
                                                temperature=settings["temperature"],
                                                a1=arrhenius_constant)
                
                    logging.debug(f"line {j}: Arrhenius calculation used for k")
                else:
                    logging.debug(f"line {j}: Arrhenius constant calculated using transition state theory")
                    #read rest of required information
                    e0, symmetry_factor, amount_of_charges = (reaction_file_df.loc[j, 'Activation energy'],
                                                              reaction_file_df.loc[j, 'Symmetry factor'],
                                                              reaction_file_df.loc[j, 'Amount of charges'])
                    
                    activation_energy = potential_dependent_activation_energy(
                                                            e0,
                                                            settings["potential"],
                                                            symmetry_factor,
                                                            amount_of_charges,
                                                            )
                    #calculate k based on electrochemistry
                    k = calculate_rate_constant(activation_energy,
                                                temperature=settings["temperature"],
                                                a1=arrhenius_constant)
                
            rate_constants.append(k)
            
            #handle reactions
            # replace tabs by spaces:
            reaction = " ".join(reaction.split("\t"))
            # divide reaction into two sides based on reaction marker
            educt_site, product_site = reaction.split("-->")
            
            # split product and educt sites into individual reactants
            raw_educts = educt_site.split(" + ")
            raw_products = product_site.split(" + ")
            # create output containers for sub dictionary
            educts = []
            products = []
            educt_prefix = np.ones(len(raw_educts))
            product_prefix = np.ones(len(raw_products))
            # get prefixes for every reactant
            for m, educt in enumerate(raw_educts):
                prefix, reactant = _extract_prefix(educt, find_prefix)
                educts.append(reactant)
                educt_prefix[m] = prefix
    
            for m, product in enumerate(raw_products):
                prefix, reactant = _extract_prefix(product, find_prefix)
                products.append(reactant)
                product_prefix[m] = prefix
            
            # store results
            reaction_map[n] = {
                "educts": educts,
                "products": products,
                "educt prefix": educt_prefix,
                "product prefix": product_prefix,
            }
            #increase n 
            n += 1
        #transform rate_constants to array
        rate_constants = np.array(rate_constants)
    
    return rate_constants, reaction_map


def potential_dependent_activation_energy(e0,
                                          potential,
                                          symmetry_factor,
                                          amount_of_charges=1,
                                          ):
    """
    Calculate the potential-dependent shift of the activation barrier based on
    a linear interpolation.
    If multiple charges are transferred, it is assumed that they all behave
    independently and thus the potential dependency scales again linearly.
    The calculation is inspired by Eq.(5) and (6) of Zijlstra et al.,
    Electrochimica Acta 335 (2020), 135665, https://doi.org/10.1016/j.electacta.2020.135665
    
    Note that if more than one array-like input parameter is given, they must
    be of the same shape.

    Parameters
    ----------
    e0 : int, float, or array-like
        Activation energy at 0 V in J
    potential : int, float, or array-like with n entries
        The applied potential in forward (reductive) direction in V.
    symmetry_factor : int, float, or array-like
        Value from 0-1 describing the transition rate in the reaction.
        If the value is 0, no potential-dependency will be calculated.
        If it is between 0 and 1, its absolute value must add up to one with
        the corresponding backwards reaction.
        Note that the sign of the symmetry factor depends on the direction of
        the reaction (reduction -> positive, oxidation -> negative)
    amount_of_charges : int, float, or array-like, optional
        Amount of electrons transferred in the direction. In case of positive
        charges, use negative values!
        The default is +1.

    Returns
    -------
    E_potential: float, or array-like 
        The potential-depedent activation barrier(s) in J/mol.

    """ 
    #translate e0 from kJ/mol to J/mol
    e0 *= 1000
    #calculate potential
    E_potential = e0 + symmetry_factor * amount_of_charges * e * N_A * potential

    return E_potential


def calculate_rate_constant(activation_energy,
                            temperature=298.15,
                            a1=None):
    """
    Calculates a rate constant based on temperature and activation energy.
    If an empirical prefactor a1 is given, an Arrhenius-like calculation is performed.
    If not, transition state theory is used assuming a transition coefficient of 1.

    Parameters
    ----------
    activation_energy : int, float or array-like,
        Activation energy in units of J/mol.
    temperature : int, float or array-like, optional
        Temperature in Kelvin. The default is 298.15 K.
    a1 : int, float, array-like or None, optional
        Amplitude of the exponential function. This can be an empirical Arrhenius
        prefactor, or None. In the latter case, transition state theory is used
        to calculate it as a function of temperature. Here, the ideal case of a
        transition coefficient of 1 is assumed.
        The default is None.

    Returns
    -------
    rate_constant : float, or array-like
        Calculated temperature-dependent rate constant.

    """
    #use transition state theory to calculate a1 if no better value is given
    if a1 is None:
        #use transition state theory to calculate a1 if no better value is given.
        a1 = kB * temperature / h
    #calculate exponential term
    energy_ratio = - activation_energy / (gas_constant * temperature)
    #calculate rate constant
    rate_constant = a1 * np.exp(energy_ratio)
       
    return rate_constant


def _make_reactant_set(reaction_map):
    """
    Create a set containing every species as string

    Parameters
    ----------
    reaction_map : dict
        nested dictionary, as returned by create_reactions


    Returns
    -------
    reactant_set : unique list

    """
    reactant_set = set()

    for n in reaction_map:
        for site in ("educts", "products"):
            for reactant in reaction_map[n][site]:
                reactant_set.add(reactant)

    return sorted(list(reactant_set))


def get_initial_concentrations(
    concentration_set,
    concentrations="Water_concentrations and gvalues.txt",
    sep=None,
    default_concentration=0.0,
    default_gval=0.0,
):
    """
    Create a dictionary with the initial concentrations and g values given in a file

    Parameters
    ----------
    concentration_set : iterable
        Iterable object containing a preferably unique list of the required species names.

    concentrations : str, optional
        Name of the concentration file. The default is 'Concentrations.txt'.

    sep : str or None, optional.
        Seperator of the columns in conenctrations. If None (default), a space or tab is assumed.

    default_concentration : float, optional
        Value to assign to species that are not found in concentrations. The default is 0.

    default_gval : float, optional
        Value to assign to species that are not found in concentrations. The default is 0.


    Returns
    -------
    dict_reactants : dict
        Dictionary containing the species name as key and the concentrations as values.

    dict_gvals : dict
        Dictionary containing the species name as key and the g value as values.

    """
    # create a concenctration dictionary:
    dict_reactants = {}
    # create a gvalue dictionary:
    dict_gvals = {}
    # read concentration file
    with open(concentrations) as concentration_file:
        concentration_list = concentration_file.read().split("\n")

    # assign concentration to dict_reactants. Skip first line as header
    for n, species in enumerate(concentration_list[1:], start=2):
        # try to split the values into two. Empty lines and header line will fail
        splitresult = species.split(sep)
        try:
            species, c0, gval = splitresult
        except ValueError as v:
            if len(splitresult) == 2:
                (species, c0), gval = splitresult, default_gval
                logging.warning(
                    f"No g-value found for {species}. Thus, it is set to {default_gval}."
                )
            elif len(splitresult) == 0:
                logging.error(f"Nothing found in line {n}. Thus, it is ignored.")
                continue
            else:
                logging.critical(f"Problems with {species}.")
                raise v

        # assign species as key and the concentration as value
        dict_reactants[species] = float(c0)

        # convert the generation value per single eV instead of per 100 ev
        gval = float(gval) / 100
        # assign species as key and generation value as value
        dict_gvals[species] = gval

    # check whether all species have been dealt with.
    for species in concentration_set:
        if species not in dict_reactants:
            logging.info(
                f"{species} is not defined in {concentrations}. {default_concentration} M is assigned as initial concentration and {default_gval} is set as G-value."
            )
            # Assign default_concentration otherwise
            dict_reactants[species] = default_concentration
            dict_gvals[species] = default_gval

    return dict_reactants, dict_gvals


def get_initial_conditions(
    concentration_set,
    conditions="Water_concentrations and gvalues.txt",
    sep=None,
    default_concentration=0.0,
    default_gval=0.0,
    default_gval_alt=0.0,
    default_solubility=100.0,
    default_diff_coeff=2.3e-9,
):
    """
    Create a dictionary with the initial concentrations, solubility parameters, and (up to two sets of) g values given in a file.

    Parameters
    ----------
    concentration_set : iterable
        Iterable object containing a preferably unique list of the required species names.

    conditions : str, optional
        Name of the concentration file. The default is 'Concentrations.txt'.

    sep : str or None, optional.
        Seperator of the columns in conenctrations. If None (default), a space or tab is assumed.

    default_concentration : float, optional
        Value to assign to species that are not found in conditions. The default is 0.

    default_solubility : float, optional
        Value to assign to species that are not found in conditions in percent. The default is 100.
        
    default_gval : float, optional
        Value to assign to species that are not found in conditions. The default is 0.
        
    default_diff_coeff : float, optional
        Value to assign to species that are not found in conditions. The default is 2.3e-9.


    Returns
    -------
    dict_reactants : dict
        Dictionary containing the species name as key and the concentrations as values.
        
    dict_gvals : dict
        Dictionary containing the species name as key and the g value as values.

    dict_gvals_alt : dict
        Dictionary containing the species name as key and the g value as values.

    dict_diff_coeff : dict
        Dictionary containing the species name as key and the solubility as values.
        
    dict_solubility : dict
        Dictionary containing the species name as key and the solubility as values.

    """
    
    cond_df = pd.read_csv(conditions, sep="\t")
    
    # check whether there are parameters for all species.
    for species in concentration_set:
        # add species if species is not yet in dataframe
        if species not in cond_df['species'].values:
            cond_df.loc[len(cond_df)] = {'species': species}
            
            
    column_dct = {'species':                                ['None', 'Default','species'],    # unit, assigned default value, variable name
                  'c0 / M':                                 ['M', default_concentration, 'reactants'], 
                  'g values / (molecules / 100 eV)':        ['(molecules / 100 eV)', default_gval, 'gvals'],
                  'g values alt / (molecules / 100 eV)':    ['(molecules / 100 eV)', default_gval_alt, 'gvals_alt'],
                  'solubility / %':                         ['%', default_solubility, 'solubility'],
                  'diffusioncoefficient / m^2/s':           ['m^2/s', default_diff_coeff, 'diff_coeff'],
                  }
    
    
    # Replace nan values in dataframe with default values
    # Extract default values from the dictionary (2nd entry in each list)
    default_values = {key: value[1] for key, value in column_dct.items()}
    
    # notify user about replaced values
    nan_indices = {col: cond_df[cond_df[col].isna()].index.tolist() for col in cond_df.columns}
    for key, nan_index_list in nan_indices.items():
        if len(nan_index_list) > 0:
            logging.warning(
                f"No value for {key} found for {[cond_df['species'][i] for i in nan_index_list]}. \n \
                    Thus they are set to {column_dct[key][1]} {column_dct[key][0]}."
                    )
    
    # Fill NaN values with the default values
    cond_df.fillna(value=default_values, inplace=True)
    
    
    # g values as per eV instead of per 100 eV
    cond_df['g values / (molecules / 100 eV)'] = cond_df['g values / (molecules / 100 eV)'] / 100
    cond_df['g values alt / (molecules / 100 eV)'] = cond_df['g values alt / (molecules / 100 eV)'] / 100
    # solubility as float instead of %
    cond_df['solubility / %'] = cond_df['solubility / %'] / 100
    
    
    # Rename dataframe columns
    # Create a mapping of old column names to display names (3rd entry in the list)
    rename_map = {key: value[2] for key, value in column_dct.items()}
    
    # Rename the columns in the DataFrame
    cond_df.rename(columns=rename_map, inplace=True)
    
    
    column_dicts = {'dict_' + col: cond_df.set_index("species")[col].to_dict() for col in cond_df.columns if col != "species"}
    
    dict_reactants  = column_dicts['dict_reactants']
    dict_gvals      = column_dicts['dict_gvals']
    dict_gvals_alt  = column_dicts['dict_gvals_alt']
    dict_diff_coeff = column_dicts['dict_diff_coeff']
    dict_solubility = column_dicts['dict_solubility']
       
    ###!!!
    # check whether all species have been dealt with.
    for species in concentration_set:
        if species not in dict_reactants:
            logging.warning(
                f"{species} is not defined in {conditions}. Default values are set."
            )
            # Assign default_concentration otherwise
            dict_reactants[species] = default_concentration
            dict_gvals[species]     = default_gval
            
    return dict_reactants, dict_gvals, dict_gvals_alt, dict_diff_coeff, dict_solubility



def _make_reaction_array(rate_constants, reaction_map, reactants, reactants_arr):
    """
    Create numpy array out of dictionary for the reaction rates in uM as input for ODE solver

    Parameters
    ----------
    rate_constants :
    concentrations : str, optional
        Name of the concentration file. The default is 'Concentrations.txt'.

    reaction_map :  dict
        nested dictionary, as returned by create_reactions

    reactants : iterable
        unique list of reactants

    reactants_arr :

    Returns
    -------
    r : numpy array containing concentrations

    """
    # map current concentration with species
    dict_reactants = {sp: val for sp, val in zip(reactants, reactants_arr)}

    # set up numpy arrays for the reaction rates and to store the products
    r = np.zeros(len(reaction_map))

    # iterate over reaction map
    for n in reaction_map:

        # educt concentration array:
        educt_c = np.zeros(len(reaction_map[n]["educts"]))

        # iterate over all educts in reaction
        for i, educt in enumerate(reaction_map[n]["educts"]):
            educt_c[i] = dict_reactants[educt]

        # multiply rate constants * PRODUCT(educt_concentration ** prefix)
        r[n] = rate_constants[n] * np.prod(educt_c ** reaction_map[n]["educt prefix"])

    return r


def _make_product_array(reactions, reactants, linked_reactions):
    """
    Calculates the concentration change of every equation per step

    Parameters
    ----------
    reactants : iterable
        unique list of reactants

    linked_reactions : dict
        Dictionary containing the stochiometric information for a species with
        every reaction involved


    Returns
    -------
    products : numpy array containing the new concentrations

    """

    # define output array
    products = np.zeros(len(reactants))

    # fill output array by iterating over all species. linked_reactions is used as a roadmap for connections
    for i, species in enumerate(reactants):

        # get coupled equation numbers and stochiometric infos from linked_reactions
        for eq, prefix in zip(
            linked_reactions[species]["educts"],
            linked_reactions[species]["educt prefix"],
        ):

            # calculate the respespective change
            products[i] = products[i] - prefix * reactions[eq]

        # get coupled equation numbers and stochiometric infos from linked_reactions
        for eq, prefix in zip(
            linked_reactions[species]["products"],
            linked_reactions[species]["product prefix"],
        ):
            # calculate the respespective change
            products[i] = products[i] + prefix * reactions[eq]

    return products


def _get_product_connections(reactants, reaction_map):
    """
    Creates all links for a species between the reactions

    Parameters
    ----------
    reactants : iterable
        unique list of reactants

    reaction_map : dict
        nested dictionary, as returned by create_reactions


    Returns
    -------
    linked_reactions : dict
        Dictionary containing the stochiometric information for a species with
        every reaction involved

    """

    # create output dictionary
    linked_reactions = {}

    # iterate over all existing species
    for species in reactants:

        # create a nested dictionary to store the equation number where the species occurs as educt or product
        linked_reactions[species] = {
            "educts": [],
            "educt prefix": [],
            "products": [],
            "product prefix": [],
        }

        # iterate over all reactions
        for i in reaction_map:

            # check if species occurs as educt in equation number i
            for agent in reaction_map[i]["educts"]:
                if species == agent:
                    # store equation number
                    linked_reactions[species]["educts"].append(i)
                    # get stochiometric info
                    stoch = reaction_map[i]["educt prefix"][
                        reaction_map[i]["educts"].index(species)
                    ]
                    # store stochiometric info
                    linked_reactions[species]["educt prefix"].append(stoch)

            # check if species occurs as educt in equation number i
            for agent in reaction_map[i]["products"]:
                if species == agent:
                    # store equation number if this is the case
                    linked_reactions[species]["products"].append(i)
                    # get stochiometric info
                    stoch = reaction_map[i]["product prefix"][
                        reaction_map[i]["products"].index(species)
                    ]
                    # store stochionetric info
                    linked_reactions[species]["product prefix"].append(stoch)

    return linked_reactions


def _make_species_arrays(reactants, dict_reactants, dict_gvals, dict_gvals_alt={}):
    """
    Converts the initial concentrations and g values into numpy arrays

    Parameters
    ----------
    reactants :  iterable
        unique list of all reactants

    dict_reactants : dict
        Dictionary containing the stochiometric information for a species with
        every reaction involved

    dict_gvals : dict
        dictionary containing the g values for every species
        
    dict_gvals_alt : dict
        dictionary containing alternative g values for every species

    Returns
    -------
    reactants_arr : numpy array

    gval_arr : numpy array
    
    gval_alt_arr : numpy array
    """

    reactants_arr = np.zeros(len(reactants))
    gval_arr = np.zeros(reactants_arr.shape)
    gval_alt_arr = np.zeros(reactants_arr.shape)

    for n, species in enumerate(reactants):

        reactants_arr[n] = dict_reactants[species]
        gval_arr[n] = dict_gvals[species]
        gval_alt_arr[n] = dict_gvals_alt[species]

    return reactants_arr, gval_arr, gval_alt_arr


def _make_species_arrays_mobility(reactants, dict_reactants, dict_solubility, dict_diff_coeff):
    """
    Converts the initial concentrations and g values into numpy arrays

    Parameters
    ----------
    reactants :  iterable
        unique list of all reactants

    dict_reactants : dict
        Dictionary containing the stochiometric information for a species with
        every reaction involved

    dict_solubility : dict
        dictionary containing the solubility value for every species
        
    dict_diff_coeff : dict
        dictionary containing the diffusion coefficient value for every species

    Returns
    -------
    reactants_arr : numpy array

    solubility_arr : numpy array
    """

    reactants_arr = np.zeros(len(reactants))
    solubility_arr = np.ones(reactants_arr.shape)
    diff_coeff_arr = np.zeros(len(reactants))

    for n, species in enumerate(reactants):

        reactants_arr[n] = dict_reactants[species]
        solubility_arr[n] = dict_solubility[species]
        diff_coeff_arr[n] = dict_diff_coeff[species]

    return reactants_arr, solubility_arr, diff_coeff_arr


def radiolysis(
    t,
    reactants_arr,
    gval_arr,
    reactants,
    linked_reactions,
    rmap,
    rate_constants,
    dose_rate=2e7,
    density=1.0,
    flow_velocity=0,
    inflow_concentrations=0, 
    length=1e-9,
    solubility_arr=1,
    diffusioncoefficient=0
):
    """
    Converts the initial concentrations and g values into numpy arrays

    Parameters
    ----------
    t: float,
        dummy variable to feed to ODE solver. As all reaction rates are
        normalized to 1 s into account, this input is automatically dealt with
        by the solver

    reactants_arr: numpy array
        array containing the concentration prior to calculation

    gval_arr: numpy array
        array containing the g values

    reactants :  iterable
        unique list of all reactants

    linked_reactions : dict
        Dictionary containing the stochiometric information for a species with
        every reaction involved

    rmap : dict
        dictionary containing the g values for every species

    rate_constants

    dose_rate: float, optional
        dose rate in Gy/s. Default is 2e7.

    density: float, optional
        Liquid density in kg/l. Default is 1 for water (assumed by N. Schneider).
        
    flow_velocity: float, optional
        Liquid flow velocity in m/s. Default is 0.
    
    inflow_concentrations: 0 or numpy array, optional; mandatory if species are mobile (i.e. length parameter supplied)
        Concentration of surrounding (or flowing) medium. Default 0
    
    length: float, optional
        Length of irradiated area in flow direction. default 1e-7.
        
    solubility_arr: float or array containing 0 and 1
        array defining for each species whether it is soluble (influenced by the flow), represented by a 1, or insoluble (not influenced by the flow), represented by a zero
        
    diffusioncoefficient: float (or array?), optional
        Diffusioncoefficient in m^2/s. Default 0. A reasonable value for small molecules/atoms in water is ~2.1e-9
        

    Returns
    -------
    products : numpy array containing the updated concentration
    """

    # create reaction rate array
    reactions = _make_reaction_array(rate_constants, rmap, reactants, reactants_arr)

    # create product array based on reaction set parameters
    products = _make_product_array(reactions, reactants, linked_reactions)  
    # products is a change of concentration (mol/s); here from the reactions happening

    # conversion from Gy/s to eV/(m³s)
    dose_rate = dose_rate / (e * 0.001)
    
    """
    # add to radiolysis --> 1e3: Conversion from m³ to L
    products = products + density * dose_rate * gval_arr / (N_A * 1e3) \
        + (flow_velocity / beamlength) * (inflow_concentrations-reactants_arr) \
        + (diffusioncoefficient / ((beamlength/2)**2)) * (inflow_concentrations-reactants_arr)
        # the second line adds flow
        # the third line adds diffusion
    """
    
    # add to radiolysis --> 1e3: Conversion from m³ to L
    products = products + density * dose_rate * gval_arr / (N_A * 1e3)
    # products is a change of concentration (mol/s); here from the reactions happening + radiolytically created species
    
    # add flow
    # add diffusion
    
    if length!=0.0:
        products += (flow_velocity / length) * (inflow_concentrations-reactants_arr)*solubility_arr \
            + (diffusioncoefficient / length**2) * (inflow_concentrations-reactants_arr)*solubility_arr
        # the first line adds flow; unit: m/s / m * (mol-mol)*1 --> mol/s
        # the second line adds diffusion; unit: m^2/s / m^2 * (mol-mol)*1 --> mol/s
    
    return products


def _string_there(string, dct):
    """
    Checks if string is used as a key in dict and prints a message if so.

    Parameters
    ----------
    string : str

    dct: dict
    """
    if string in dct:

        logging.warning("{} is redefined".format(string))


def read_settings(setting_file):
    """
    input: Filename as string

    return: dictionary containing settings
    """

    def read_out_standard_line(line):

        # split line at equal sign
        __, val = line.split("=")
        # convert to float and return value
        return float(val)

    # create storage dictionary:
    setting_dct = {}

    # create keys to look for in the lines
    dose_rate_key = "dose rate"
    liquid_thickness_key = "liquid thickness"
    t_end_key = "duration"
    mean_free_path_key = "mean free path"
    stopping_power_key = "stopping power"
    rtol_key = "rtol"
    atol_key = "atol"
    density_key = "density"
    dwell_time_key = "dwell time"
    return_time_key = "return time"
    flow_velocity_key = "flow velocity"
    length_key = "length"
    diffusioncoefficient_key = 'diffusioncoefficient'
    temperature = "temperature"
    potential = "potential"
    g2_g1_ratio = "g2 to g1 ratio"

    # get settings file
    with open(setting_file) as file:
        # create a list with every line in file as string
        setting_list = file.readlines()

    ##iterate over lines and extract features:
    for line in setting_list:

        # skip if starts with comment:
        if line.startswith("#"):
            continue

        # extract dose rate line
        if line.lower().startswith(dose_rate_key):
            # inform user if dose rate is used twice
            _string_there(dose_rate_key, setting_dct)
            # line looks like "dose rate = 1 e/A2s" -> extract unit and value
            __, e_flux = line.split("=")
            e_flux, flux_unit = e_flux.split()
            # convert value to float and store
            setting_dct["e Flux"] = float(e_flux)
            # store unit for documentation
            setting_dct["Flux unit"] = flux_unit
            # skip key-loop

        else:
            # check every other key
            for key in [
                liquid_thickness_key,
                t_end_key,
                mean_free_path_key,
                stopping_power_key,
                rtol_key,
                atol_key,
                density_key,
                dwell_time_key,
                return_time_key,
                flow_velocity_key,
                length_key,
                diffusioncoefficient_key,
                temperature,
                potential,
                g2_g1_ratio
            ]:

                if line.lower().startswith(key):
                    # inform user that key is used twice
                    _string_there(key, setting_dct)
                    # convert value to float and store
                    setting_dct[key] = read_out_standard_line(line)
                    # break loop if key is found
                    # break

    # deal with dose_rate unit conversion:
    if flux_unit == "Gy/s":
        dose_rate = setting_dct["e Flux"]

    elif flux_unit == "e/A2s":
        dose_rate = e_flux_to_dose_rate_angstrom(
            setting_dct["e Flux"],
            setting_dct[liquid_thickness_key],
            stopping_power=setting_dct["stopping power"],
            mean_free_path=setting_dct["mean free path"],
        )
    elif flux_unit == "e/nm2s":
        dose_rate = (
            e_flux_to_dose_rate_nm(
                setting_dct["e Flux"],
                setting_dct[liquid_thickness_key],
                stopping_power=setting_dct["stopping power"],
                mean_free_path=setting_dct["mean free path"],
            )
        )
    else:
        raise ValueError(
            'Dose rate unit "{}" not recognized. It must be one of the following: {}'.format(
                flux_unit, ["Gy/s", "e/A2s", "e/nm2s"]
            )
        )
    print('Read from file: E-flux:', str(e_flux), ' Dose rate :', str(dose_rate))

    # store calculated dose rate
    setting_dct[dose_rate_key] = dose_rate

    return setting_dct


interesting_species = []
def display_results(times, concentrations, setting_dct, reactants, elements=None, 
                    interesting_species=[], file_type='xlsx'):

    f, ax = plt.subplots(figsize=(5, 5), dpi=150)

    relevant_c = setting_dct["atol"] * 2
    ##filter used species
    # create output container
    used_species = []
    # loop over solver results. Range is used because the species position is relevant
    for i in range(0, concentrations.shape[0]):
        # store species if its max. concentration reaches RELEVANT_C

        if concentrations[i].max() >= relevant_c:
            used_species.append(i)

    # define line layout
    amount_of_used_species = len(used_species)
    color = plt.cm.nipy_spectral(np.linspace(0, 1, amount_of_used_species + 2))[1:-1]
    symbols = cycle(
        [
            "-",
            "--",
            "-.",
            ":",
        ]
    )
    # title layouting

    liquid_thickness_str = setting_dct["liquid thickness"] * 1e9
    if liquid_thickness_str == int(liquid_thickness_str):
        liquid_thickness_str = int(liquid_thickness_str)
    else:
        liquid_thickness_str == round(liquid_thickness_str, 1)

    for i, n in enumerate(used_species):
        
        species = reactants[n]
        
        if elements is not None:
            #skip non-selected elements
            cond = [elem in species for elem in elements]
            if interesting_species:
                cond.append(species in interesting_species)
            if not any(cond):
                continue  
        
        elif interesting_species:
            if species not in interesting_species:
                continue
        
        s = next(symbols)
        ax.loglog(times, concentrations[n], ls=s, c=color[i], label=make_nice_label(species))

                  
                
    # lims = ax.get_ylim()
    
    """
    # define title string based on e flux unit:
    flux_unit = setting_dct["Flux unit"]
    dose_rate = setting_dct["dose rate"]
    e_flux = setting_dct["e Flux"]

    decimal_str = convert_decimal_to_latex(f'{Decimal(dose_rate):.2e}')

    if flux_unit == "e/A2s":
        title_str = f"{e_flux} {angstr} and {liquid_thickness_str} nm H$_2$O $\\Rightarrow$ {decimal_str} {gys}"
    elif flux_unit == "e/nm2s":
        title_str = f"{e_flux} {nms} and {liquid_thickness_str} nm H$_2$O $\\Rightarrow$ {decimal_str} {gys}"
    elif flux_unit == "Gy/s":
        title_str = f"{decimal_str} {gys}"
  
    potential_part = ""
    if "potential" in setting_dct:
        potential_part = f"{setting_dct['potential']}V_"
        title_str += f"; {setting_dct['potential']} V" + "$_\\mathrm{SHE}$"
    
    temperature_part = ""
    if "temperature" in setting_dct:
        temperarture_part = f"{setting_dct['temperature']}K_"
        title_str += f" at {setting_dct['temperature']} K" 
        
        
    gys = "$\\frac{\\mathrm{Gy}}{\\mathrm{s}}$"
    angstr = "$\\frac{\\mathrm{e}^-}{\\mathrm{\\AA}^2\\mathrm{s}}$"
    nms = "$\\frac{\\mathrm{e}^-}{\\mathrm{nm}^2\\mathrm{s}}$"
    
    """
    
    # new title string
    iter_variables = setting_dct["iter variables"]
    
    
    e_flux = setting_dct["e Flux"]    
    flux_unit = setting_dct["Flux unit"]
    dose_rate = setting_dct["dose rate"]
    
    
    decimal_str = convert_decimal_to_latex(f'{Decimal(dose_rate):.2e}')
    
    if len(iter_variables) == 0:
        unit_tex = setting_dct["Flux unit"]
        title_str = f"{decimal_str} {unit_tex}"
        
    else:
        title_str = ""
        for variable in iter_variables:
            unit_tex = variable_to_unit_dct[variable][1]
            title_str += f"{setting_dct[variable]} {unit_tex}; "
            
    ax.set(
        ylim=(relevant_c, concentrations.max() * 2),
        xlim=(1e-12, times.max()),
        xlabel="Time / s",
        ylabel="Concentration / M",
        title=title_str,
    )


    # ax.grid()
    # create amount of columns based on amount of species
    legend_columns = round(amount_of_used_species / 10)
    if legend_columns == 0:
        legend_columns = 1

    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), ncol=legend_columns)
    ax.tick_params(direction="in", which="both")
    plt.show()

    # save results

    flux_unit_save = "".join(flux_unit.split("/"))
    reaction_set_name = setting_dct["reaction set name"].rsplit("\\", 1)[-1]
    if "\\" in reaction_set_name:
        reaction_set_name = reaction_set_name.rsplit("\\", 1)[-1]
    concentrations_name = setting_dct["c0 name"]
    
    #notify displayed elements
    elem_str = ""
    if elements is not None:
        elem_str = "_" + "".join(elements)
    
    
    output_folder = f"{setting_dct['output folder']}"
    
    if len(iter_variables) == 0:
        
        savename = os.path.join(output_folder, 
            f"{reaction_set_name}_{concentrations_name}_{e_flux}_{flux_unit_save}"
        )
    else:
        savename = ""
        for variable in iter_variables:
            unit_plain = variable_to_unit_dct[variable][0]
            
            
            savename += f"{setting_dct[variable]}_{unit_plain}_"
        savename = os.path.join(output_folder, 
                                savename.rsplit("_", 1)[0]
                                )
        
    # save figure
    for end in ['png', 'svg', 'pdf']:
        f.savefig(
            f"{savename}{elem_str}.{end}",
            bbox_inches="tight",
            dpi=300,
            transparent=True
        )
        #print(f"saved {savename}{elem_str}.{end}")

    plt.close("all")

    ##use Pandas to save the data to a human readable file
    # create DataFrame
    df = pd.DataFrame(concentrations.T, columns=reactants)
    # Assign time in seconds as index
        # if there was a STEM iteration calculated before: 
        # add the last timestamp of the previous step to the timestamp
    df["time / s"] = times
    df.set_index("time / s", inplace=True)
    
    
    # save    
    if file_type == "xlsx":
        df.to_excel(savename + ".xlsx")
    else:
        df.to_csv(savename + f".{file_type}")
    
    return df


def make_nice_label(string, custom_label={}):
    """
    Converts a chemical formula to LaTeX notation using its math environment. 
    
    If 'string' ends with '-' or '+', they are automatically interpreted as superscript.
    Numbers within 'string' are subscripted. To superscript those, they must be prefixed by '^' or '**'.
    Radical labels are not supported.
    
    Examples:
        
        make_nice_label('H2SO4-') returns 'H$_{2}$SO$_{4}^-$'
        
        make_nice_label('C60+') returns 'C$_{60}^+$'
        
        make_nice_label('Au^3+') returns 'Au$^{3+}$' but make_nice_label('Au3+') returns 'Au$_{3}^+$'
    
    'eh-' and 'es-' are designed using custom labels that can be overwritten by the 'custom_label' argument.

    Parameters
    ----------
    string : str
        Chemical formula.
    custom_label : dict, optional
        Dictionary that allows custom output by using string as key and the custom output as value. The default is {}.

    Returns
    -------
    out_str : str
        LaTeX-formatted chemical formula.
    """
    ignore ={}
    ignore['eh-'] = 'e$_\\mathrm{h}^-$'
    ignore['es-'] = 'e$_\\mathrm{s}^-$'
    
    #merge with custom_label
    ignore = {**ignore, **custom_label}
    
    if string in ignore:
        
        out_str = ignore[string]
    
    else:
        
        string = '^'.join(string.split('**'))
        
        numbers = re.compile('\d+')
        used_nums = []

                #charge part
        if '^' in string:
            body, charge = string.split('^')
            charge = '$^{' + charge + '}$'
        else:
            body = string
            charge = ''

        if body.endswith('-') or body.endswith('+'):
            charge = f'$^{body[-1]}$'
            body = body[:-1]

        out_str = body

        for i in numbers.finditer(out_str):

            found = i.group()

            if found in used_nums:
                continue

            replacement = '$_{'+found+'}$'
            out_str = replacement.join(out_str.split(found))

            used_nums.append(found)

        out_str += charge
        out_str = ''.join(out_str.split('$$'))
    
    return out_str


def main(reaction_info, initial_info, setting_dct, extended_reaction_info=[],
         interesting_species=[]):
    """
    Main execution program based on the inputs of text files.
    """
    # unpack variable
    gval_arr, reactants,linked_reactions, rmap, rate_constants = reaction_info
    
    
    if len(extended_reaction_info)>0:
        gval_arr_alt, inflow_concentrations, solubility_arr, diff_coeff_arr = extended_reaction_info
    
        g2_g1_ratio = setting_dct["g2 to g1 ratio"]
        print("G2G1: ", g2_g1_ratio)
    
        #if g1_g2_ratio != 1.0: 
            #gval_arr = g1_g2_ratio * gval_arr + (1-g1_g2_ratio) * gval_arr_alt #1.0 --> 100% G1 + 0% G2, 0.0 --> 0% G1 + 100% G2; allows only for INTERpolation
        if g2_g1_ratio != 0.0:
            gval_arr = gval_arr + g2_g1_ratio * (gval_arr_alt - gval_arr)  #0.0 --> 100% G1 + 0% G2, 1.0 --> 0% G1 + 100% G2; allows only for INTER and EXTRApolation
    
    initial_concentrations = initial_info

    # define time interval
    time_span = 0.0, setting_dct["duration"]

    # starting timestamp
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    logging.info(f"Simulation started. Current Time ={current_time}")
    

    # perform ODE solving
    solutions = solve_ivp(
        lambda t, c0: radiolysis(
            t,
            c0,
            gval_arr,
            reactants,
            linked_reactions,
            rmap,
            rate_constants,
            dose_rate=setting_dct["dose rate"],
            density=setting_dct["density"],
            flow_velocity=setting_dct["flow velocity"],
            inflow_concentrations=initial_concentrations,
            length=setting_dct["length"],
            solubility_arr=solubility_arr,
            diffusioncoefficient=diff_coeff_arr,
        ),
        time_span,
        initial_concentrations,
        method="LSODA",
        rtol=setting_dct["rtol"],
        atol=setting_dct["atol"],
    )

    # retreive output parameters of interest
    times = solutions.t
    # convert solution to micro molar
    concentrations = solutions.y

    # end timestamp
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")

    logging.info(f"Simulation ended. Current Time ={current_time}")

    # plot and store
    df = display_results(times, concentrations, setting_dct, reactants, 
                         interesting_species=interesting_species)

    return df


def main_loop(reaction_set, init_conditions, setting_file, infl_concentrations=None, interesting_species=[]):   # XXX check obsolensence; some parameters need to be integrated via extended_reaction_info
    """
    Main execution program based on the inputs of text files.
    """
    
    #define file names
    
    #load settings
    setting_dct = read_settings(setting_file)

    fluxes = np.logspace(11, 15, 3)
    #fluxes = [1e12]
    doses = fluxes
    
    flow_velocities = [0]  # XXX

    lengths = np.logspace(-10, -6, 3)
    #beam_lengths = [1e-5]
    g1_g2_ratios = [1] #[1, 0.5, 0] #np.linspace(0, 1, 3)#[::-1]
    
    for g1_g2_ratio in g1_g2_ratios:
        for length in lengths:
            for flow_v in flow_velocities:
                for flux, dose in zip(fluxes, doses):
                    
                    logging.info(f'Start with: {length} m; {flow_v} m/s; {dose} Gy/s')
            
                    setting_dct['e Flux'] = flux
                    setting_dct['dose rate'] = dose
                    setting_dct['flow velocity'] = flow_v
                    setting_dct['length'] = length
                    setting_dct['g1_g2_ratio'] = g1_g2_ratio
                    
                    #split rate constants and create a dictionary containing educts and products
                    rate_constants, rmap = create_reactions(reaction_set, setting_dct)
                    #create a set of reactants
                    reactants = _make_reactant_set(rmap)
                    #link initial concentrations to reactants and assign 0 if no other value is given
                    #dict_reactants, dict_gvals = get_initial_concentrations(reactants,
                    #                                                        concentrations=init_concentrations)
                    dict_reactants, dict_gvals, dict_gvals_alt, dict_diff_coeff, dict_solubility = get_initial_conditions(reactants,
                                                                                                         conditions=init_conditions)
                    
                               
                    linked_reactions = _get_product_connections(reactants, rmap)
                    initial_concentrations, gval_arr, gval_arr_alt = _make_species_arrays(reactants,
                                                                            dict_reactants,
                                                                            dict_gvals, 
                                                                            dict_gvals_alt)
                    
                    inflow_concentrations, solubility_arr, diff_coeff_arr = _make_species_arrays_mobility(reactants,
                                                                    dict_reactants, dict_solubility, dict_diff_coeff)
                    
                    """
                    ### old implementation relying on separate "IN" file for flow 
                    # assume that inflow_concentrations are same as init_concentrations and all species are soluble (flushable) if not otherwise specified
                    if not infl_concentrations:
                        inflow_concentrations = initial_concentrations
                        insoluble_arr = np.ones(len(reactants))
                    else:
                        dict_reactants_flow, dict_reactants_insoluble = get_initial_concentrations(reactants,
                                                                                concentrations=infl_concentrations,
                                                                                default_gval=1.0)  # default_gval goes to insolubility value now
                                   
                        inflow_concentrations, insoluble_arr = _make_species_arrays_mobility(reactants,
                                                                        dict_reactants_flow, dict_reactants_insoluble)
    
                    """   
                    
                    #define time interval
                    time_span = 0., setting_dct['duration'] 
                    
                    #starting timestamp
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    print("Start Time =", current_time)
                    
                    #perform ODE solving
                    solutions = solve_ivp(lambda t, c0: radiolysis(t,
                                                                   c0,
                                                                   gval_arr,
                                                                   reactants,
                                                                   linked_reactions, 
                                                                   rmap,
                                                                   rate_constants,
                                                                   dose_rate= setting_dct['dose rate'],
                                                                   density=setting_dct["density"],
                                                                   flow_velocity=setting_dct['flow velocity'],
                                                                   inflow_concentrations=inflow_concentrations,
                                                                   length=setting_dct["length"],
                                                                   solubility_arr=solubility_arr,
                                                                   diffusioncoefficient=diff_coeff_arr,
                                                                   g1_g2_ratio=g1_g2_ratio,
                                                                   gval_arr_alt=gval_arr_alt,
                                                                   ),
                                      time_span,
                                      initial_concentrations,
                                      method= 'LSODA',
                                      rtol=setting_dct['rtol'],
                                      atol=setting_dct['atol'],
                                      #max_step=10e-9   #XXX
                                      )
                    
                    #retreive output parameters of interest
                    times = solutions.t
                    concentrations = solutions.y
                    
                    #end timestamp
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    
                    print("End Time =", current_time)                
                    
                    #create names for plotting
                    setting_dct['reaction set name'] = reaction_set.split('.txt')[0]
                    #setting_dct['c0 name']  = init_concentrations.split('.txt')[0]
                    
                    setting_dct['c0 name'] = init_conditions.split(".txt")[0]
                    
                    #plot and store
                    display_results(times, concentrations, setting_dct, reactants, interesting_species=interesting_species)

  
def main_T_loop(reaction_set, init_concentrations, setting_file, temperatures):  # XXX! should be obsolete
    """
    Main execution program based on the inputs of text files.
    """

    # define file names

    # load settings
    setting_dct = read_settings(setting_file)

    # split rate constants and create a dictionary containing educts and products
    rate_constants, rmap = create_reactions(reaction_set, setting_dct)
    # create a set of reactants
    reactants = _make_reactant_set(rmap)
    # link initial concentrations to reactants and assign 0 if no other value is given
    dict_reactants, dict_gvals = get_initial_concentrations(
        reactants, concentrations=init_concentrations
    )

    linked_reactions = _get_product_connections(reactants, rmap)
    initial_concentrations, gval_arr = _make_species_arrays(
        reactants, dict_reactants, dict_gvals
    )

    # define time interval
    time_span = 0.0, setting_dct["duration"]
    for temperature in temperatures:
        
        setting_dct["temperature"] = temperature
        
        #recalculate rate constants with new temperature:
        rate_constants, rmap = create_reactions(reaction_set, setting_dct)
        
        # starting timestamp
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        logging.info(f"Simulation started. Current Time ={current_time}")
    
        # perform ODE solving
        solutions = solve_ivp(
            lambda t, c0: radiolysis(
                t,
                c0,
                gval_arr,
                reactants,
                linked_reactions,
                rmap,
                rate_constants,
                dose_rate=setting_dct["dose rate"],
                density=setting_dct["density"],
            ),
            time_span,
            initial_concentrations,
            method="LSODA",
            rtol=setting_dct["rtol"],
            atol=setting_dct["atol"],
        )
    
        # retreive output parameters of interest
        times = solutions.t
        # convert solution to micro molar
        concentrations = solutions.y
    
        # end timestamp
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
    
        logging.info(f"Simulation ended. Current Time ={current_time}")
    
        # create names for plotting
        setting_dct["reaction set name"] = reaction_set.split(".txt")[0]
        setting_dct["c0 name"] = init_concentrations.split(".txt")[0] + f'_{temperature}K'
    
        # plot and store
        display_results(times, concentrations, setting_dct, reactants)
   
    
def load_ground_state(reaction_set, init_conditions, setting_file):  
    # load settings
    setting_dct = read_settings(setting_file)

    # split rate constants and create a dictionary containing educts and products
    rate_constants, rmap = create_reactions(reaction_set, setting_dct)
    # create a set of reactants
    reactants = _make_reactant_set(rmap)
    # link initial concentrations to reactants and assign 0 if no other value is given
    dict_reactants, dict_gvals, dict_gvals_alt, dict_diff_coeff, dict_solubility = get_initial_conditions(reactants,
                                                                                         conditions=init_conditions)
    linked_reactions = _get_product_connections(reactants, rmap)
    initial_concentrations, gval_arr, gval_arr_alt = _make_species_arrays(reactants,
                                                            dict_reactants,
                                                            dict_gvals, 
                                                            dict_gvals_alt)
    inflow_concentrations, solubility_arr, diff_coeff_arr = _make_species_arrays_mobility(reactants,
                                                    dict_reactants, dict_solubility, dict_diff_coeff)


    # create names for plotting
    setting_dct["reaction set name"] = reaction_set.split(".txt")[0]
    setting_dct["c0 name"] = init_conditions.split(".txt")[0]

    reaction_info = (gval_arr,
                     reactants,
                     linked_reactions,
                     rmap,
                     rate_constants)
    
    extended_reaction_info = (gval_arr_alt,
                              inflow_concentrations,
                              solubility_arr,
                              diff_coeff_arr)

    initial_info = initial_concentrations
    
    return reaction_info, initial_info, setting_dct, extended_reaction_info
    

def looper(reaction_info, reaction_set, initial_conc, setting_dct, extended_reaction_info,
           looped_variables={}, interesting_species=[]):
    
    now = datetime.now()
    datetime_str = f"{now.strftime('%Y-%m-%d_%H-%M')}"
    
    output_folder = "Simulation Results"
    os.makedirs(output_folder, exist_ok=True)
    setting_dct["output folder"] = output_folder
    
    
    if len(looped_variables) == 0:
        df = main(reaction_info, initial_conc, setting_dct, extended_reaction_info,
                  interesting_species)
        
    else:
        #XXX! create subfolder with date, save setting_dct therein in a readable format

        output_folder = os.path.join(output_folder, datetime_str)
        os.makedirs(output_folder, exist_ok=True)
        
        setting_dct["output folder"] = output_folder
        setting_dct["iter variables"] = looped_variables.keys()
        # Save to a text file
        with open(os.path.join(output_folder, "initial_conditions.txt"), "w") as file:
            for key, value in setting_dct.items():
                file.write(f"{key}: {value}\n")
            
            file.write("\n \n" + \
                       "species c_0 \t\tG1 \t\tG2 \tsolubility[100%] \n")
            
            
            for key, c_0, g1, g2, sol in zip(reaction_info[1], 
                                        initial_conc, 
                                        reaction_info[0], 
                                        extended_reaction_info[0],
                                        extended_reaction_info[2]
                                        ):   # species and concentration
                file.write(f"{key}: \t{c_0:.1e} \t{g1:.1e} \t{g2:.1e} \t{sol}\n")

        # Use itertools.product to generate all combinations
        for combination in product(*looped_variables.values()):
            # Map the combination back to the variable names
            variable_mapping = dict(zip(looped_variables.keys(), combination))
            
            for key, value in variable_mapping.items():
                setting_dct[key] = value
                
            if "temperature" in looped_variables.keys():  # when the temperature is varied, the rate constants need to be recalculated
                rate_constants, rmap = create_reactions(reaction_set, setting_dct)
            
            print("using G2G1:", setting_dct["g2 to g1 ratio"])
            df = main(reaction_info, initial_conc, setting_dct, extended_reaction_info,
                      interesting_species)
            
    return df


def calculate_temperature_from_g2_g1_ratio(target_temperature, g1_temperature, g2_temperature):
    """
    Computes the g2_g1_ratio to be used for temperature-dependent g-value generation, e.g. to be fed to looper.
    
    Parameters
    ----------
    target_temperature : int, float or array-like
        Target temperature in K.
    g1_temperature : int or float
        Temperature in K matching the G-values in the g1 column of the init_conditions file.
    g2_temperature : int or float
        Temperature in K matching the G-values in the g2 column of the init_conditions file..

    Returns
    -------
    g2_g1_ratio: int, float, or array
        g2_g1_ratio for calculating g-values at a given target_temperature.

    """
    g2_g1_ratio = (target_temperature - g1_temperature) / (g2_temperature - g1_temperature) 
   
    return g2_g1_ratio

#%%
if __name__ == "__main__":

    ### HUI
    reaction_set = "Reaction Sets\Temperature-dependent_Water_Reactions+SparseHAuCl4.txt"
    init_conditions = "conditions_H2O.txt"
    setting_file = "settings.txt"
    interesting_species = []
    
    doserates = np.logspace(3,5,2)
    lengths = np.logspace(-9, -3, 2)

    
    looped_variables = {"dose rate":[1e9],# doserates,
                        "length": [0],
                        "g2 to g1 ratio": [0.0, 1.0],
                        }
    
    
    ## do temperature-dependent simulations
    temperatures = np.linspace(20, 100, 9) + 273.15
    g2_g1_ratio = calculate_temperature_from_g2_g1_ratio(temperatures,
                                                         293.15, #T of g1
                                                         373.15 #T of G2
                                                         )
    looped_variables["temperature"] = temperatures
    looped_variables["g2 to g1 ratio"] = g2_g1_ratio
    init_conditions = "conditions_HAuCl4.txt"

    
    #note that diffusion coefficients are regarded temperature-independent here.
    
    
    reaction_info, initial_conc, setting_dct, extended_reaction_info = load_ground_state(reaction_set, 
                                                                                         init_conditions, 
                                                                                         setting_file)
    # "$\\frac{\\mathrm{Gy}}{\\mathrm{s}}$"
            
    ## used to generate filename; triggerword in settings_dct: [unit_plain, unit_tex]
    variable_to_unit_dct = {"e Flux":           [setting_dct["Flux unit"].replace("/s", "s-1"), setting_dct["Flux unit"]],
                            "liquid thickness": ["m_thick", r"$\mathrm{m}_t$"],
                            "duration":         ["s", "s"],
                            "flow velocity":    ["ms-1", r"$\frac{m}{s}$"],
                            "length":           ["m_long", r"$\mathrm{m}_l$"],
                            "diffusioncoefficient": ["m2s", r"\frac{m^2}{s}"],
                            "temperature":      ["K", r"$K$"],
                            "potential":        ["V", r"$V$"],
                            "density":          ["kgl-1", r"$\frac{kg}{L^{-1}}$"],
                            "mean free path":   ["m_mfp", r"$\mathrm{m}_\lambda$"],
                            "stopping power":   ["MeV(cm)2/g", "$\\frac{{MeV(cm)^2}{g}$"],
                            "dose rate":        [setting_dct["Flux unit"].replace("/s", "s-1"), setting_dct["Flux unit"]],
                            "g2 to g1 ratio":   ["g2Tg1", "$_{g2tg1}$"]
                            }

    last_df = looper(reaction_info, reaction_set, initial_conc, setting_dct, extended_reaction_info,
                         looped_variables, interesting_species)
       
    #df = main(reaction_info, initial_conc, setting_dct, interesting_species)
    
