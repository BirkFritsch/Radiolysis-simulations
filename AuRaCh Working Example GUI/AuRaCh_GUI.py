"""
This is the supplementary code to the publication:

Fritsch et al., Radiolysis‐Driven Evolution of Gold Nanostructures –
 Model Verification by Scale Bridging In Situ Liquid‐Phase Transmission Electron
 Microscopy and X‐Ray Diffraction, Advanced Science, 2022, DOI:10.1002/advs.202202803

If you find this tool helpful, please cite our work.

This code is published under a GPL-3.0 licence.
"""

import sys
import re
import os 
import numpy as np
import logging
import pandas as pd
from scipy.constants import N_A, e
from scipy.integrate import solve_ivp
import matplotlib
from matplotlib import pyplot as plt
from datetime import datetime
from itertools import cycle
from decimal import Decimal
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
matplotlib.use('Qt5Agg')

from PyQt5 import QtCore,QtGui
from PyQt5.QtCore import QThread,QCoreApplication,Signal, Slot



from PyQt5.QtWidgets import (
    QApplication,
    QComboBox,
    QFormLayout,
    QLineEdit,
    QStackedLayout,
    QVBoxLayout,
    QWidget,
    QDialogButtonBox,
    QHBoxLayout,
    QRadioButton,
    QFileDialog,
    QPushButton,
    QMessageBox,
    QTextEdit,
    QMenuBar,
    QMenu,
    QCheckBox,
    QTableView,
    QMainWindow,
    QWidget,
    QTextBrowser,
    QPlainTextEdit,
    QGridLayout,
    QLabel,
   
)
from PyQt5.QtWidgets import QMessageBox,QStatusBar
from PyQt5 import QtWidgets
from PyQt5.QtGui import QIntValidator,QTextCursor
from PyQt5.QtCore import Qt
from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot


log = logging.getLogger("Log")

logging.basicConfig(level=logging.WARNING)

matplotlib.use('Qt5Agg')



def e_flux_to_dose_rate(flux, thickness, stopping_power = 2.36, mean_free_path=3.8e-7,t=1):
    """
    INPUT:
    
    flux (float or array-like) in e/(Angstroem²s)
    Thickness (float or array-like) in m
    Mean free path (float or array-like) in m
    stopping_power (float or array-like) in MeV(cm)²/g
    
    RETURN:
    
    dose rate (float or array-like) in Gy/s 
    """

    dose_rate = 1e5 * e* (stopping_power/(1e-20*t)) * flux * (1 + thickness/mean_free_path)
    
    return dose_rate



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
    #match regular expression to string
    re_result = re_object.search(string)
    #if re_result is not None, this will be truthy
    if re_result:
        #split prefix and reactant
        try:
            prefix, reactant = string.split()
        except ValueError as v:
            logging.critical(f'ValueError for {string}. This means, no patterb are found, although expected by the regular expression.')
            raise v
        #convert prefix to float
        prefix = float(prefix)
    
    else:
        #dummy prefix for later calculation
        prefix = 1.
        #remove spaces
        reactant = string.split()[0]
    
    return prefix, reactant



def create_reactions(filename):
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
    #read reaction file
    
    with open(filename) as reaction_file:
        #create a list with every line in file as string
        reaction_list = reaction_file.readlines()

    #create output containers
    rate_constants = np.zeros(len(reaction_list))
    reaction_map = {}
    
    #define re.Pattern object for prefix extraction
    find_prefix = re.compile('^\s?\d+\.?\d*\s+')
    
    #do the work
    for n, reaction in enumerate(reaction_list):
        
        #replace tabs by spaces:
        reaction = ' '.join(reaction.split('\t'))
        
        #split the rate constant k off the reaction site. The asterix accounts for comments. Thus, k is a tuple.
        reaction_part, *k = reaction.split(';')
        
        #divide reaction into two sides based on reaction marker
        try:
            educt_site, product_site = reaction_part.split('-->')
        except ValueError as v:
            logging.critical(f'Could not handle the reaction file in line {n}: {reaction}. I skip line {n} and continue.')
            raise v
        #store k
        try:
            #first index in k is the rate constant.
            rate_constants[n] = float(k[0])
        except ValueError:
            logging.error(f'Could read the rate constant in line {n}: {reaction}. It is set to 0, which effectively disables the reaction.')
            rate_constants[n] = 0.
        #split product and educt sites into individual reactants
        raw_educts = educt_site.split(' + ')
        raw_products = product_site.split(' + ')
        
        #create output containers for sub dictionary
        educts = []
        products = []
        educt_prefix = np.ones(len(raw_educts))
        product_prefix = np.ones(len(raw_products))
        
        #get prefixes for every reactant
        for m, educt in enumerate(raw_educts):
            prefix, reactant = _extract_prefix(educt, find_prefix)
            educts.append(reactant)
            educt_prefix[m] = prefix
        
        for m, product in enumerate(raw_products):
            prefix, reactant = _extract_prefix(product, find_prefix)
            products.append(reactant)
            product_prefix[m] = prefix
        
        #store results
        reaction_map[n] = {
                           'educts':educts,
                           'products':products,
                           'educt prefix':educt_prefix,
                           'product prefix':product_prefix,
                           }
        
    
    return rate_constants, reaction_map



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
        for site in ('educts', 'products'):
            for reactant in reaction_map[n][site]:
                reactant_set.add(reactant)
                
    return sorted(list(reactant_set))
        


def get_initial_concentrations(concentration_set, concentrations='Water_concentrations and gvalues.txt',
                               sep=None, default_concentration = 0., default_gval=0.):
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
    #create a concenctration dictionary:
    dict_reactants = {}
    #create a gvalue dictionary:
    dict_gvals = {}
    #read concentration file
    with open(concentrations) as concentration_file:
        concentration_list = concentration_file.read().split('\n')
        
    #assign concentration to dict_reactants. Skip first line as header
    for n, species in enumerate(concentration_list[1:], start=2):
        #try to split the values into two. Empty lines and header line will fail
        splitresult = species.split(sep)
        try:
            species, c0, gval = splitresult
        except ValueError as v:
            if len(splitresult) == 2:
                (species, c0), gval = splitresult, default_gval
                logging.debug(f'No G-value found for {species}. Thus, it is set to {default_gval}.')
            elif len(splitresult) == 0:
                logging.info(f'Nothing found in line {n}. Thus, it is ignored.')
                continue
            else:
                logging.critical(f'Problems with {species}.')
                raise v
            
        #assign species as key and the concentration as value
        dict_reactants[species] = float(c0)
        
        #convert the generation value per single eV instead of per 100 ev
        gval = float(gval) / 100
        #assign species as key and generation value as value
        dict_gvals[species] = gval
        
    #check whether all species have been dealt with.
    for species in concentration_set:
        if species not in dict_reactants:
            logging.warning(f'{species} is not defined in {concentrations}. {default_concentration} M is assigned as initial concentration and {default_gval} is set as G-value.')           
            #Assign default_concentration otherwise   
            dict_reactants[species] = default_concentration
            dict_gvals[species] = default_gval
            
    return dict_reactants, dict_gvals



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
    #map current concentration with species
    dict_reactants = {sp:val for sp, val in zip(reactants, reactants_arr)}
    
    #set up numpy arrays for the reaction rates and to store the products
    r = np.zeros(len(reaction_map))

    #iterate over reaction map
    for n in reaction_map:
        
        #educt concentration array:
        educt_c = np.zeros(len(reaction_map[n]['educts']))
        
        #iterate over all educts in reaction
        for i, educt in enumerate(reaction_map[n]['educts']):
            educt_c[i] = dict_reactants[educt]
        
        #multiply rate constants * PRODUCT(educt_concentration ** prefix)
        r[n] = rate_constants[n] * np.prod(educt_c**reaction_map[n]['educt prefix'])
        
#    #Unit conversion: Molar to micromolar
#    r = r * 1e-6
#        
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
    
    #define output array
    products = np.zeros(len(reactants))    
    
    #fill output array by iterating over all species. linked_reactions is used as a roadmap for connections
    for i, species in enumerate(reactants):
       
        #get coupled equation numbers and stochiometric infos from linked_reactions
        for eq, prefix in zip(linked_reactions[species]['educts'],
                               linked_reactions[species]['educt prefix']):

            #calculate the respespective change
            products[i] = products[i] - prefix * reactions[eq]
            
            
        #get coupled equation numbers and stochiometric infos from linked_reactions
        for eq, prefix in zip(linked_reactions[species]['products'],
                               linked_reactions[species]['product prefix']):
            #calculate the respespective change
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
    
    #create output dictionary
    linked_reactions={}
    
    #iterate over all existing species
    for species in reactants:
        
        #create a nested dictionary to store the equation number where the species occurs as educt or product
        linked_reactions[species] = {'educts':[],
                                     'educt prefix':[],
                                     'products':[],
                                     'product prefix':[]}
        
        #iterate over all reactions
        for i in reaction_map:
            
            #check if species occurs as educt in equation number i
            for agent in reaction_map[i]['educts']:
                if species == agent:
                    #store equation number
                    linked_reactions[species]['educts'].append(i)
                    #get stochiometric info
                    stoch = reaction_map[i]['educt prefix'][reaction_map[i]['educts'].index(species)]
                    #store stochiometric info
                    linked_reactions[species]['educt prefix'].append(stoch)
                    
            #check if species occurs as educt in equation number i
            for agent in reaction_map[i]['products']:
                if species == agent:
                    #store equation number if this is the case
                    linked_reactions[species]['products'].append(i)
                    #get stochiometric info
                    stoch = reaction_map[i]['product prefix'][reaction_map[i]['products'].index(species)]
                    #store stochionetric info
                    linked_reactions[species]['product prefix'].append(stoch)

                
    return linked_reactions



def _make_species_arrays(reactants, dict_reactants, dict_gvals):
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
            
    Returns
    -------
    reactants_arr : numpy array
    
    gval_arr : numpy array 
    """
    
    reactants_arr = np.zeros(len(reactants))
    gval_arr = np.zeros(reactants_arr.shape)
    
    for n, species in enumerate(reactants):
        
        reactants_arr[n] = dict_reactants[species]
        gval_arr[n] = dict_gvals[species]
        
    return reactants_arr, gval_arr



def radiolysis(t, reactants_arr, gval_arr, reactants, linked_reactions, rmap,
               rate_constants, dose_rate=2e7, density=1.):
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
            
    Returns
    -------
    products : numpy array containing the updated concentration
    """
    
    
    #create reaction rate array
    reactions = _make_reaction_array(rate_constants, rmap, reactants, reactants_arr)
    
    #create product array based on reaction set parameters    
    products = _make_product_array(reactions, reactants, linked_reactions)
    
    #conversion from Gy/s to eV/(m³s)
    dose_rate = dose_rate/(e*0.001)
    
    #add to radiolysis --> 1e3: Conversion from m³ to L 
    products = products + density * dose_rate * gval_arr / (N_A * 1e3)

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
        
        logging.warning('{} is redefined'.format(string))



def read_settings(setting_file):
    """
    input: Filename as string
    
    return: dictionary containing settings
    """
    
    def read_out_standard_line(line):
        
        #split line at equal sign
        __, val = line.split('=')
        #convert to float and return value
        return float(val)
    
    #create storage dictionary:
    setting_dct = {}
    
    #create keys to look for in the lines
    dose_rate_key = 'dose rate'
    liquid_thickness_key = 'liquid thickness'
    t_end_key = 'duration'
    mean_free_path_key = 'mean free path'
    stopping_power_key = 'stopping power'
    rtol_key = 'rtol'
    atol_key = 'atol'
    density_key = 'density'
    
    #get settings file
    with open(setting_file) as file:
        #create a list with every line in file as string
        setting_list = file.readlines()
    
    
    ##iterate over lines and extract features:
    for line in setting_list:
        
        #skip if starts with comment:
        if line.startswith('#'):
            continue
       
        #extract dose rate line 
        if line.lower().startswith(dose_rate_key):
            #inform user if dose rate is used twice
            _string_there(dose_rate_key, setting_dct)
            #line looks like "dose rate = 1 e/A2s" -> extract unit and value
            __, e_flux = line.split('=')
            e_flux, flux_unit = e_flux.split()
            #convert value to float and store
            setting_dct['e Flux'] = float(e_flux)
            #store unit for documentation
            setting_dct['Flux unit'] = flux_unit
            #skip key-loop
            
        else:
            #check every other key
            for key in [liquid_thickness_key,
                        t_end_key,
                        mean_free_path_key,
                        stopping_power_key,
                        rtol_key,
                        atol_key,
                        density_key,
                        ]:
                
                if line.lower().startswith(key):
                    #inform user that key is used twice
                    _string_there(key, setting_dct)
                    #convert value to float and store
                    setting_dct[key] = read_out_standard_line(line)
                    #break loop if key is found
                    #break
                
            
    #deal with dose_rate unit conversion:
    if flux_unit == 'Gy/s':
        dose_rate = setting_dct['e Flux']

    elif flux_unit == 'e/A2s':
        dose_rate = e_flux_to_dose_rate(setting_dct['e Flux'],
                                        setting_dct[liquid_thickness_key])
    elif flux_unit == 'e/nm2s':
        dose_rate = e_flux_to_dose_rate(setting_dct['e Flux'],
                                        setting_dct[liquid_thickness_key])*100
    else:
        raise ValueError('Dose rate unit "{}" not recognized. It must be one of the following: {}'.format(flux_unit, ['Gy/s', 'e/A2s', 'e/nm2s']))                                        
    
    #store calculated dose rate
    setting_dct[dose_rate_key] = dose_rate                                        

            
    return setting_dct



def display_results(times, concentrations, setting_dct, reactants):
    
    
    f, ax = plt.subplots(figsize = (10,5))
    
    relevant_c = setting_dct['atol']*2
    ##filter used species
    #create output container
    used_species = []
    #loop over solver results. Range is used because the species position is relevant
    for i in range(0,concentrations.shape[0]):
        #store species if its max. concentration reaches RELEVANT_C
        
        if concentrations[i].max() >= relevant_c:
            used_species.append(i)      
    
    #define line layout
    amount_of_used_species = len(used_species)
    color = plt.cm.nipy_spectral(np.linspace(0,1,amount_of_used_species+2))[1:-1]
    symbols = cycle(['-','--', '-.', ':',])
    #title layouting
    gys = '$\\frac{\\mathrm{Gy}}{\\mathrm{s}}$'
    angstr = '$\\frac{\\mathrm{e}^-}{\\mathrm{\\AA}^2\\mathrm{s}}$'
    nms = '$\\frac{\\mathrm{e}^-}{\\mathrm{nm}^2\\mathrm{s}}$'
    liquid_thickness_str = setting_dct['liquid thickness']*1e9
    if liquid_thickness_str == int(liquid_thickness_str):
        liquid_thickness_str = int(liquid_thickness_str)
    else:
        liquid_thickness_str == round(liquid_thickness_str,1)
        
    for i, n in enumerate(used_species):
        
        species=reactants[n]
        
        
        s = next(symbols)
        ax.loglog(times, concentrations[n], ls=s, c=color[i], label=species)
    #lims = ax.get_ylim()
    
    #define title string based on e flux unit:
    flux_unit = setting_dct['Flux unit']
    dose_rate = setting_dct['dose rate']
    e_flux = setting_dct['e Flux']
    if flux_unit == 'e/A2s':
        title_str=f'{e_flux} {angstr} and {liquid_thickness_str} nm H$_2$O $\\Rightarrow$ {Decimal(dose_rate):.2e} {gys}'
    elif flux_unit == 'e/nm2s':
        title_str=f'{e_flux} {nms} and {liquid_thickness_str} nm H$_2$O $\\Rightarrow$ {Decimal(dose_rate):.2e} {gys}'
    elif flux_unit == 'Gy/s':
        title_str=f'{e_flux} {gys}'

    ax.set(ylim = (relevant_c, concentrations.max()*10),
            xlim = (1e-12,None),
            xlabel="Time / s",
            ylabel="Concentration / M",
            title = title_str)
    ax.tick_params(direction='in', which = 'both')

    #ax.grid()    
    #create amount of columns based on amount of species
    legend_columns = round(amount_of_used_species/10)
    if legend_columns == 0:
        legend_columns = 1
        
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=legend_columns)
 
    plt.show()
    plt.tight_layout()
    #save results
    
    flux_unit_save = ''.join(flux_unit.split("/"))
    reaction_set_name = setting_dct['reaction set name']
    concentrations_name = setting_dct['c0 name']
    
    
    #save figure
    f.savefig(f'{reaction_set_name}_{concentrations_name}_{e_flux}_{flux_unit_save}_Output.png', bbox_inches='tight', dpi=300)
    
    
    ##use Pandas to save the data to a human readable file
    #create DataFrame
    #save results
    if setting_dct['check_box'] != 0:

        ##use Pandas to save the data to a human readable file
        #create DataFrame
        df=pd.DataFrame(concentrations.T,columns=[r+ ' / M' for r in reactants])
        #Assign time in seconds as index
        df['time / s'] = times
        df.set_index('time / s', inplace=True)
        #create name for saving
         
        # savename = f'{reaction_set_name}_{concentrations_name}_{e_flux}_{flux_unit_save}_Output'
        
        if setting_dct['save_name'] != '':
          savepath = setting_dct['save_path']
          df.to_csv( os.path.join(savepath))
         
          setting_dct['file_ready'] =True
        else:
          
          df.to_csv(os.path.join(os.getcwd()+'/'+'Test.csv'))
          setting_dct['file_ready'] =True
         
        df['time / s'] = times
        
    else:
        df=pd.DataFrame(concentrations.T,columns=[r+ ' / M' for r in reactants])
        #Assign time in seconds as index
        df['time / s'] = times
        df.set_index('time / s', inplace=True)
        df['time / s'] = times
    return df



   
def main_GUI(reaction_set, init_concentrations, setting_dct):
    """
    Main execution program to excecute via a graphical user interface.
    
    INPUT:
        
        reaction_set: string or path-like object that refers to the reactions text file
        
        init_concentrations: string or path-like object that refers to the initial concentrations and g-values text file
        
        setting_dct: dict, as created by the GUI. It must contain the following keyword arguments:
        'e Flux', 'Flux unit', 'duration', 'liquid thickness', 'mean free path',
        'stopping power', 'rtol', 'atol','save_name','save_path'
    """
    # gui_instance.wait_message()
    #adjust dose rate based on unit
        #deal with dose_rate unit conversion:
    if setting_dct['Flux unit'] == 'Gy/s':
        dose_rate = setting_dct['e Flux']

    elif setting_dct['Flux unit'] == 'e/A2s':
        dose_rate = e_flux_to_dose_rate(setting_dct['e Flux'],
                                        setting_dct['liquid thickness'])
    elif setting_dct['Flux unit'] == 'e/nm2s':
        dose_rate = e_flux_to_dose_rate(setting_dct['e Flux'],
                                        setting_dct['liquid thickness'])*100
    setting_dct['dose rate'] = dose_rate
    
    
    #split rate constants and create a dictionary containing educts and products
    rate_constants, rmap = create_reactions(reaction_set)
    #create a set of reactants
    reactants = _make_reactant_set(rmap)
    #link initial concentrations to reactants and assign 0 if no other value is given
    dict_reactants, dict_gvals = get_initial_concentrations(reactants,
                                                            concentrations=init_concentrations)
    
    linked_reactions = _get_product_connections(reactants, rmap)
    initial_concentrations, gval_arr = _make_species_arrays(reactants,
                                                            dict_reactants,
                                                            dict_gvals)
    
    #define time interval
    time_span = 0., setting_dct['duration'] 
    
    #starting timestamp
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)
    
        
    #perform ODE solving
    solutions = solve_ivp(lambda t, c0: radiolysis(t,
                                                   c0,
                                                   gval_arr,
                                                   reactants,
                                                   linked_reactions, 
                                                   rmap,
                                                   rate_constants,
                                                   dose_rate= setting_dct['dose rate'],
                                                   density=setting_dct['density']
                                                   ),
                      time_span,
                      initial_concentrations,
                      method= 'LSODA',
                      rtol=setting_dct['rtol'],
                      atol=setting_dct['atol'])
    
    #retreive output parameters of interest
    times = solutions.t

    concentrations = solutions.y
    
    #end timestamp
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    
    print("End Time =", current_time)

    
    #create names for plotting
    setting_dct['reaction set name'] = reaction_set.split('.txt')[0].split('/')[-1]
    setting_dct['c0 name']  = init_concentrations.split('.txt')[0].split('/')[-1]
    
    #plot and store
    
    
    # df = display_results(times, concentrations, setting_dct, reactants)
    # return df
    return times, concentrations, setting_dct, reactants



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

    
    parts = re_object.findall(reactant)
    
    for p in parts:
        #remove parenthesis part:
        reactant = ''.join(reactant.split(p))
    
    return reactant



def get_charge(reactant):
    
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
  
        ####WARNING##### MUST BE ADJUSTED
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
    
    #replace tabs by spaces:
    reaction_string = ' '.join(reaction_string.split('\t'))
    #split educts and products into separate lists
    educts, products = reaction_string.split('-->')
    educts = educts.split(' + ')
    products = products.split(' + ')
    
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



def check_reaction_file(file_str):
    
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
    df.to_excel(file_str.split('.')[0]+'_checked.xlsx')
    
    return df



def check_concentration(df):
    
    # try:
        
    #     df = pd.read_csv(file)
        
        
    # except:
    #     df = pd.read_excel(file)
    
    #set time column as index
    df.set_index('time / s', inplace=True)

    #remove unit if necessary:
    species = [s.split()[0] for s in df.columns] 
    
    
    #elements as list:
    elements = [*extract_atoms(species)]
    
    
    #do fancy numpy stuff to get a quick result
    out = np.zeros((df.shape[0], len(elements)))
    
    #iterate over species
    for c,s in zip(df.columns,species):

        #get dictionary with atoms in s
        atoms = extract_atoms([s])
        
        #iterate over each atom
        for a in atoms:
            
            #get position in element list
            idx = elements.index(a)
            #add concentration times the atom for each column
            out[:,idx] += df[c] * atoms[a]
    
    #summarize in DataFrame to store column name
    df_out = pd.DataFrame(out, columns=elements)
    
    #match index with simulation
    df_out['time / s'] = df.index
    df_out.set_index('time / s', inplace=True)
    
    #show & store results:
    
    # savename = ''.join(file.split('.')[:-1])
    savename = "consistency"
    f, ax = plt.subplots(figsize=(5,5), dpi=100)
    for s in df_out:
        
        ax.loglog(df_out[s], label=s)
    
    ax.tick_params(direction='in', which = 'both')
    ax.legend(loc=0)
    ax.set(xlabel= 'time / s' , ylabel='Concentration / M')
    f.tight_layout()
    plt.show()
    
    # f.savefig(savename+'_checked.png', dpi=600, bbox_inches='tight')
    # plt.close('all')
    
    #df_out.to_excel(savename+'_checked.xlsx')
    
    # return df_out


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
    out=check_reaction_file(filename)
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


# if __name__ == '__main__':
    
#     filename = 'Water_Schneider_MV4.txt'
    
#     df_stoch = check_reaction_file(filename)
    
    
    
    #%%
    # conc_file = 'Test.csv'
    
    # df_atoms = check_concentration(conc_file)



   
class TableModel(QtCore.QAbstractTableModel):
    """
    Defining a subclass of QAbstractTableModel for showing incorrect recations as table.
    """

    def __init__(self, data):
        super(TableModel, self).__init__()
        self._data = data

    def data(self, index, role):
        if role == Qt.DisplayRole:
            value = self._data.iloc[index.row(), index.column()]
            return str(value)

    def rowCount(self, index):
        return self._data.shape[0]

    def columnCount(self, index):
        return self._data.shape[1]

    def headerData(self, section, orientation, role):
      
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return str(self._data.columns[section])

            if orientation == Qt.Vertical:
                return str(self._data.index[section])
          

class TableWindow(QtWidgets.QMainWindow):
     """
       Defining a Window for showing the table
     """
     def __init__(self, data,parent=None):
        # super(TableWindow, self).__init__(parent)
        super().__init__()
        container = QWidget()

        self.l = QVBoxLayout()
        self.textEdit = QtWidgets.QTextEdit()
        self.textEdit.setReadOnly(True)
        self.textEdit.setSizeAdjustPolicy(self.textEdit.AdjustToContents)
        self.l.addWidget(self.textEdit)
        
        layout_view = QVBoxLayout()
        layout_view.addLayout(self.l)
        
        self.table = QtWidgets.QTableView()
        layout_view.addWidget(self.table)
        container.setLayout(layout_view)
        self.model = TableModel(data)
        self.table.setModel(self.model)
        self.table.resizeColumnsToContents()
        self.table.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        
       
        
        self.table.setGeometry(QtCore.QRect(0,0,1500,2500))
    
        self.setCentralWidget(container)
        
     
    
        
class Window(QWidget):
    """
    Defining main Window for GUI
    """
    
    def __init__(self):
        QWidget.__init__(self)
        
        
        Height = 750
          
        # setting  the fixed Height of window
        self.setFixedHeight(Height)
        self.values={'e Flux':'1e8','Flux unit':'Gy/s','liquid thickness':'1e-7','duration':'1e3',
                     'mean free path':'3.8e-7','stopping power':'2.36',
                     'rtol':'1e-4','atol':'1e-30','density':'1','save_path':'','save_name':'' ,'check_box':1,'simulate_start' :False,'file_ready':False}
        
        self.Gvalue = 'not specified'
        self.Reactions = 'not specified'
       
        self.df = pd.DataFrame()
        self.setWindowTitle("Settings")
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        self.field15 = QtWidgets.QTextEdit()
        self.field15.setPlainText("If you find this tool helpful please cite our work: <br/> <html><b>B. Fritsch <i>et al.</i>, Radiolysis-Driven Evolution of Gold Nanostructures – Model Verification by Scale Bridging <i>in situ</i> Liquid-Phase Transmission Electron Microscopy and X-Ray Diffraction, <i>Advanced Science</i>, 2022</b</html>")
       
        #self.field15.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.field15.setReadOnly(True)
        self.tb = QTextBrowser()
        self.tb.setOpenExternalLinks(True)
        self.tb.setStyleSheet("""QTextBrowser {background-color: '#EEEEEE' ;}""")
        layout.addWidget( self.tb)
        
        text = self.field15.toPlainText()
        self.tb.append(text)
        # self.tb.moveCursor(QTextCursor.Start)
        self.tb.append('<a href=https://doi.org/10.1002/advs.202202803>DOI:10.1002/advs.202202803</a>')
       
        self.btnup1 = QPushButton("Load Reactions",self)
        self.btnup1.setGeometry(50,10,10,20)
        self.btnup1.clicked.connect(self.get_Reactions_File)
        layout.addWidget(self.btnup1)
        
        self.layout4 = QFormLayout()       
        self.field11 = QLineEdit(self.Reactions)  
        self.layout4.addRow('File Path',self.field11)
        layout.addLayout(self.layout4)
        
        self.btn4 = QPushButton("Check Reactions",self)
        self.btn4.setGeometry(50,10,100,20)
        self.btn4.clicked.connect(self.check_Reactions_File)
        layout.addWidget(self.btn4) 
           
        self.btnup2 = QPushButton("Load Concentrations and G-Values",self)
        self.btnup2.clicked.connect(self.get_GValue_File)
        layout.addWidget(self.btnup2)
        
        self.layout5 = QFormLayout()
        self.field12 = QLineEdit(self.Gvalue)  
        self.layout5.addRow('File Path', self.field12)
        layout.addLayout(self.layout5)
        
        
        # Create and connect the combo box to switch between basic and advanced settings
        self.pageCombo = QComboBox()
        self.pageCombo.addItems(["Basic", "Advanced"])
        self.pageCombo.activated.connect(self.switchPage)
        # Create the stacked layout"""
        self.stackedLayout = QStackedLayout()
        # Create the first page
        self.Basic = QWidget()
        self.BasicLayout = QFormLayout()
        self.field1 = QLineEdit(self.values['e Flux'])
        self.BasicLayout.addRow('Dose Rate',self.field1)
      
        self.btn1 = QRadioButton(self.values['Flux unit'])
        self.btn1.setChecked(True)
        self.btn1.name = "Gy/s"
        self.btn1.toggled.connect(self.dose_rate_unit)
        
        self.btn2 = QRadioButton("e/(Å²s)")
        self.btn2.setChecked(False)
        self.btn2.name = "e/A2s"
        self.btn2.toggled.connect(self.dose_rate_unit)
        
        self.btn3 = QRadioButton("e/(nm²s)")
        self.btn3.setChecked(False)
        self.btn3.name = "e/nm2s"
        self.btn3.toggled.connect(self.dose_rate_unit)
        
        layout2 = QHBoxLayout()
        layout2.addWidget(self.btn1)
        layout2.addWidget(self.btn2)
        layout2.addWidget(self.btn3)
        layout.addWidget( self.BasicLayout.addRow('Dose Rate Unit',layout2))
        
        self.field2 = QLineEdit(self.values['liquid thickness'])  
        self.BasicLayout.addRow("Liquid Thickness (m)", self.field2)
        self.field2.setEnabled(False)
        
        self.field16 = QLineEdit(self.values['density'])  
        self.BasicLayout.addRow("Liquid Density (kg/L)",self.field16)   
        self.field16.setEnabled(False)
        
        self.field3 = QLineEdit(self.values['duration'])    
        self.BasicLayout.addRow("Duration (s)",self.field3)
        
        self.field4 = QLineEdit(self.values['mean free path'])  
        self.BasicLayout.addRow("Mean Free Path (m)",self.field4)
        self.field4.setEnabled(False)
        
        self.field5 = QLineEdit(self.values['stopping power'])    
        self.BasicLayout.addRow("Stopping Power (MeV cm²/g)",self.field5)
        self.field5.setEnabled(False)
        
        self.field6 = QLineEdit(self.values['rtol'])    
        self.BasicLayout.addRow("Relative Tolerance",self.field6)
        self.field6.setEnabled(False)
       
        self.field7 = QLineEdit(self.values['atol'])    
        self.BasicLayout.addRow("Absolute Tolerance (M)",self.field7)   
        self.field7.setEnabled(False)
        
       
        
        self.btnup4 = QPushButton("Storage Directory",self)
        self.btnup4.clicked.connect(self.save_path)
        layout.addWidget(self.btnup4)
        
           
        self.b1 = QCheckBox("save",self)
        self.b1.setChecked(True)
        self.b1.stateChanged.connect(self.clickCheckBox)
        
        
        layout3 = QHBoxLayout()
        layout3.addWidget(self.btnup4)
        layout3.addWidget(self.b1)
      
        layout.addWidget(self.BasicLayout.addRow(layout3))

        self.field13 = QLineEdit(os.getcwd())  
        layout.addWidget(self.BasicLayout.addRow('Save Path',self.field13))
          
        self.Basic.setLayout(self.BasicLayout)
        self.stackedLayout.addWidget(self.Basic)

        # Add the combo box and the stacked layout to the top-level layout
        layout.addWidget(self.pageCombo)
        layout.addLayout(self.stackedLayout)
  
             
        Runbtn = QPushButton("&Run Simulation")
        btnBox = QDialogButtonBox()
        btnBox.addButton(Runbtn, QDialogButtonBox.AcceptRole)
        btnBox.accepted.connect(self.accept)
        layout.addWidget(btnBox)
        
        
        CheckConsistencybtn = QPushButton("Check Consistency",self)
        self.checkbtnBox  = QDialogButtonBox()
        self.checkbtnBox.setEnabled(False)
        self.checkbtnBox.addButton(CheckConsistencybtn, QDialogButtonBox.AcceptRole)
        self.checkbtnBox.accepted.connect(self.check_consistency)
        layout.addWidget(self.checkbtnBox)
           
   
    def switchPage(self):
        """
        check whether basic or advanced settings is selected 
        and make fields active or deactive according to selected settings.
        """
        self.stackedLayout.setCurrentIndex(self.pageCombo.currentIndex())
        #  if advanced settings is selected these fields: "Mean Free Path (m)","Stopping Power (MeV cm²/g)"
        # "Relative Tolerance" and "Absolute Tolerance (M)" become editable
        if self.pageCombo.currentIndex()==1:
            self.field4.setEnabled(True)
            self.field5.setEnabled(True)
            self.field6.setEnabled(True)
            self.field7.setEnabled(True)
            # self.field16.setEnabled(True)
        # else if basic settings is selected above mentioned fields are not editable
        else:
            self.field4.setEnabled(False)
            self.field5.setEnabled(False)
            self.field6.setEnabled(False)
            self.field7.setEnabled(False)
            self.field16.setEnabled(False)
        
    def get_Reactions_File(self):
        """
        get reaction file uploaded using openfiledialog
        """
        self.Reactions , _ = QFileDialog.getOpenFileName(self,"Open Reaction File","",'Text Files (*.txt *.text)')
        #  if Reactions is uploaded its directory will be shown in the textbox
        if  self.Reactions:
            self.field11.setText(self.Reactions)
            self.field11.setEnabled(False)
            
         
                                      
                
    def get_GValue_File(self):
       self.Gvalue ,_ = QFileDialog.getOpenFileName(self,"Open Gvalues File","",'Text Files (*.txt *.text)')
      
       #  if Gvalues is uploaded its directory will be shown in the textbox
       if  self.Gvalue:
        
         self.field12.setText( self.Gvalue.split('.txt')[0])
         self.field12.setEnabled(False)

    def check_Reactions_File(self):
       
        """
         check whether there are duplicate reactions and 
         check whether reaction file has error or not if this is the case 
         faulty reactions will be used as data of TableWindow class
        """
   
        d = check_reaction_file(self.Reactions)
        data = pd.DataFrame(d)
        check_duplicate = display_identical_reactions(self.Reactions)
        self.duplicate_df = pd.DataFrame(check_duplicate)
        if not self.duplicate_df.empty:
         
         self.w= TableWindow(self.duplicate_df)
         self.w.setWindowTitle("Duplicated reactions")
         self.w.setGeometry(500, 100, 500, 200)
         self.w.textEdit.setPlainText("The following duplicated reactions were found:")
         self.w.show()
            
        
            
        
        self.data_table= data[(data['Element test']==False) | (data['Charge test']==False)]
        if self.data_table.empty:
     
         dlg = QMessageBox(self)
         dlg.setWindowTitle("Info")
         dlg.setText("Reaction File is correct")
         dlg.setIcon(QMessageBox.Information)
         button = dlg.exec_()   
        
        else:  
         
            self.w1= TableWindow(self.data_table)
            self.w1.setWindowTitle("Erroneous reactions")
            self.w1.setGeometry(500, 100, 1000, 300)
            self.w1.textEdit.setPlainText("The following erroneous reactions were found:")
            self.w1.show()
         
        
       
        
    def save_path(self):
       
     """
     through this method user determine save file path.
     """
     self.file_path,_ = QFileDialog().getSaveFileName(self ,'Select a directory') 
     self.values['save_path'] =  self.file_path
     self.values['save_name'] = (self.file_path.split('/')[-1])
     # if basic settings is selected
     if self.stackedLayout.currentIndex()==0:  
       self.field13.setText(self.file_path)
       self.field13.setEnabled(False)
     # else if advanced settings is selected
     else:       
      self.field14.setText(self.file_path)
      self.field14.setEnabled(False)
       
     
       
    
    def clickCheckBox(self):
     """
     check whether save checkbox is checked or not. If it is not checked result will not be saved
     """
     # if "save" check_box is not checked result will not be saved
     if self.b1.isChecked() == False:
         self.values['check_box']=0
         self.field13.setText("Result will not be stored") 
     # else if check_box is checked result will be saved at the selected directory 
     # using "storage directory" button 
     else:
         self.values['check_box']=1
         self.field13.setText(os.getcwd()) 
         
    
    def fieldCheck(self,a):   
        """
        check whether the user has entered text as input or not. 
        If this is the case warning message will be shown.

        Parameters
        ----------
        a : user entered value 
            
        """
        try:
            return float(a)
        except ValueError:
            dlg = QMessageBox(self)
            dlg.setWindowTitle("Warning")
            # check if value for 'Dose Rate' field is text
            if a==self.field1.text(): 
                dlg.setText("Only Numbers are allowed as input for e Flux")
                dlg.setIcon(QMessageBox.Warning)
                button = dlg.exec_()   
            # check if entered value for  "Liquid Thickness (m)" filed is text
            elif a==self.field2.text():
                dlg.setText("Only Numbers are allowed as input for liquid thickness")
                dlg.setIcon(QMessageBox.Warning)
                button = dlg.exec_() 
            # check if entered value for "Liquid Density (kg/L)" Field is text
            elif a==self.field16.text():
                dlg.setText("Only Numbers are allowed as input for liquid density")
                dlg.setIcon(QMessageBox.Warning)
                button = dlg.exec_()
            # check if  entered value for "Relative Tolerance" Field is text 
            elif a== self.field6.text():
                dlg.setText("Only Numbers are allowed as input for relative tolerance  ")
                dlg.setIcon(QMessageBox.Warning)
                button = dlg.exec_() 
            # check if entered value for "Absolute Tolerance (M)" Field is text 
            elif a== self.field7.text():
                dlg.setText("Only Numbers are allowed as input for absolute tolerance  ")
                dlg.setIcon(QMessageBox.Warning)
                button = dlg.exec_() 
            # check if entered value for "Mean Free Path" Field is text
            elif a== self.field4.text():
                dlg.setText("Only Numbers are allowed as input for mean free path ")
                dlg.setIcon(QMessageBox.Warning)
                button = dlg.exec_() 
            # check if entered value for "Stopping Power" Field is text
            elif a== self.field5.text():
                dlg.setText("Only Numbers are allowed as input for stopping power ")
                dlg.setIcon(QMessageBox.Warning)
                button = dlg.exec_() 
            
  
  
    
    def accept(self):
        """
       when run button is clicked this function will be executed 
       and all the values entered by user will be stored in dictionary

        """
        # check if Gvalue or Reactions is not uploaded
        if self.Gvalue == 'not specified' or self.Reactions == 'not specified':
            if self.Gvalue == 'not specified' and self.Reactions !='not specified':
                 
                dlg = QMessageBox(self)
                dlg.setWindowTitle("Warning")
                dlg.setText('"Concentrations and G Value" file is not selected!')
                button = dlg.exec_()
               
            elif  self.Reactions == 'not specified'  and self.Gvalue != 'not specified' :
                dlg = QMessageBox(self)
                dlg.setWindowTitle("Warning")
                dlg.setText('"Reactions" file is not selected')
                button = dlg.exec_()
            else:
                dlg = QMessageBox(self)
                dlg.setWindowTitle("Warning")
                dlg.setText("No file is selected")
                button = dlg.exec_()
                      
        # else if both Gvalue and Reactions are uploaded
        else:
             #  if basic settings is selected 
             if self.pageCombo.currentIndex()==0:
               self.values['e Flux'] = self.fieldCheck(self.field1.text())
               self.values['liquid thickness'] = self.fieldCheck(self.field2.text())
               self.values['density'] = 1
               self.values['duration'] = self.field3.text()
               self.values['mean free path'] = 3.8e-7
               self.values['stopping power'] = 2.36
               self.values['rtol'] = 1e-6
               self.values['atol'] = 1e-30
               # check if "save" check_box for saving the result is not checked
               if self.b1.isChecked() == False :
                 self.values['check_box']=0
                 
             #else if advanced settings is selected 
             elif self.pageCombo.currentIndex()==1:
               
               self.values['e Flux'] = self.fieldCheck(self.field1.text())
               self.values['liquid thickness'] = self.fieldCheck(self.field2.text())
               self.values['density'] = self.fieldCheck(self.field16.text())
               self.values['duration'] = self.field3.text()
               self.values['mean free path'] = self.fieldCheck(self.field4.text())
               self.values['stopping power'] = self.fieldCheck(self.field5.text())
               self.values['rtol'] = self.fieldCheck(self.field6.text())
               self.values['atol'] = self.fieldCheck(self.field7.text())
               # check if "save" check_box for saving the result is not checked 
               if self.b1.isChecked() == False :
                 self.values['check_box']=0
             # after storing all entered values by user in dictionary 
             # simulation begins and a message will be shown and thread begins to work
             if all(value !=None for value in self.values.values()):
                 self.Warning_message() 
                 self.threadworks() 
     

    def Warning_message(self):
        """
        warning message indicating start of simulation 

        """
        QMessageBox.about(self,"Warning","Simulation will start after agreeing to this message. This may take a while, please wait.")
 
    
    
    def threadworks(self):
        """
        Create a QThread object 
        """
        self.thread = QThread()
        QCoreApplication.processEvents()
        # Create a object of Worker class
        self.worker = Worker(self.Reactions,self.Gvalue,self.values)
        # Move worker to the thread
        self.worker.moveToThread(self.thread)
        # Connect signals and slots
        self.thread.started.connect(self.worker.run)
        self.worker.result.connect(self.display)
        self.worker.finished.connect(self.thread.quit)
        self.thread.start()
        self.thread.exit()
    
    
    
    def display(self, times, concentrations,setting_dct,reactants ):
        """
 
        Parameters
        ----------
        times : 
        concentrations : string or path-like object that refers to the initial concentrations and g-values text file
        setting_dct :  dict, as created by the GUI.
        reactants : string or path-like object that refers to the reactions text file

        """
        self.df = display_results(times, concentrations,setting_dct,reactants)
        self.checkbtnBox.setEnabled(True)
    
    
    
    def check_consistency(self):
        """
        check consistency if the button check consistency is clicked

        """
        check_concentration(self.df)
             
        
        
    def dose_rate_unit(self):
        """
        check which dose rate unit is selected
        """
        radioButton = self.sender()
        #check whether dose rate radio button is checked
        if radioButton.isChecked()== True:
            
            if radioButton.text()=="e/(Å²s)" or radioButton.text()=="e/(nm²s)" :
                # based on selected dose rate unit,text fields related to "liquid thickness" and "density" will be activated
                self.field2.setEnabled(True)
                self.field16.setEnabled(True)
                # in addition if advanced settings is already selected and
                # user select one of the above mentioned dose rate radio buttons,
                # then text field related to "stopping power" will be activated
                if self.stackedLayout.currentIndex()==1:
                    self.field5.setEnabled(True)
                self.values['Flux unit'] = radioButton.name
               
            # if dose rate is Gy/s 
            #text fields related to "liquid thicknes" , "density2" and "stopping power"
            # will be deactivated
            if radioButton.text()=="Gy/s":
              
                self.field2.setEnabled(False)
                self.field16.setEnabled(False)
                if self.stackedLayout.currentIndex()==1:
                    self.field5.setEnabled(False)
                self.values['Flux unit'] = radioButton.name  




class Worker(QObject):
    """
    class for passing reaction, gvalue and values entered using GUI to 
    main_gui method and get values necessary for display method
    
    """
    result = pyqtSignal(object,object,object,object)  
    finished = pyqtSignal()
    def __init__(self,Reactions,Gvalues,values):
        super().__init__()
        self.reactions,self.gvalues,self.values = Reactions,Gvalues,values

    def run(self):
        
        QCoreApplication.processEvents()
        times, concentrations, setting_dct, reactants =  main_GUI(self.reactions,self.gvalues,self.values)
        self.result.emit(times, concentrations, setting_dct, reactants)
        self.finished.emit()




if __name__ == "__main__":
  
    app = QApplication(sys.argv)
    window = Window()
    window.setWindowFlags(QtCore.Qt.WindowCloseButtonHint | QtCore.Qt.WindowMinimizeButtonHint)
    
    app.processEvents()    
    window.show()
    app.exec_()
    sys.exit(app.exec_())

        
