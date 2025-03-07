# -*- coding: utf-8 -*-
"""
tbd
"""
import logging
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import AuRaCh as au

logging.basicConfig(level=logging.WARNING)


def reaction_throughput(simulation_df, reaction_set, time):   
    
    rate_constants, reaction_map = au.create_reactions(reaction_set)           
    #get output matrix
    pseudo_rate_map = np.zeros((time.shape[0],
                                rate_constants.shape[0]))
    
    for n, k in enumerate(rate_constants):
        educt_c = np.zeros(len(reaction_map[n]["educts"]))

        for i, t in enumerate(time):
            # iterate over all educts in reaction
            for j, educt in enumerate(reaction_map[n]["educts"]):
                educt_c[j] = simulation_df[educt].loc[t]
                # multiply rate constants * PRODUCT(educt_concentration ** prefix)
            pseudo_rate = k * np.prod(educt_c ** reaction_map[n]["educt prefix"])
            #append rate
            pseudo_rate_map[i, n] = pseudo_rate
            
    return pseudo_rate_map


def analyze_reaction_throughput(reaction_set, simulation_file, save_data=False, show=True,
                                ):
    
    #load simulation file
    if simulation_file.endswith(".xlsx"):
        simulation_df = pd.read_excel(simulation_file)
    else:
        simulation_df = pd.read_csv(simulation_file)
        
    # create a set of reactant
    try:
        time = simulation_df['time / s']
        simulation_df.set_index('time / s', inplace=True)
    except KeyError as keyerror:
        if simulation_df.index.name == 'time / s':
            time = simulation_df.index
        else:
            raise keyerror
    pseudo_rate_map = reaction_throughput(simulation_df, reaction_set, time)
    reaction_throughput_df = pd.DataFrame(pseudo_rate_map.T, columns=time)   
    #for plotting
    savename = f'{reaction_set.split(".")[0]} Reaction Throughput'
    if show:
        rate_constants, __ = au.create_reactions(reaction_set) 
        Y, X = np.meshgrid(np.arange(rate_constants.shape[0]), time)
        Y += 1
        
        f, axes = plt.subplots(1, 3, figsize=(16,9), dpi=300, sharey=True, layout='constrained')
        ax0, ax1, ax2 = axes
        xlim = (1e-9, time.max())
        
        cond = (X >= xlim[0]) & (X <= xlim[1])
        
        c0=ax0.pcolormesh(X, Y, pseudo_rate_map, 
                      cmap='viridis'
                      )
        
        ax0.set(xscale='log', xlim=xlim, title='Reaction throughput / $\\mathrm{\\frac{M}{s}}$', ylabel='Reaction')
        f.colorbar(c0, ax=ax0,)
        c1=ax1.pcolormesh(X, Y, pseudo_rate_map, 
                       norm=colors.LogNorm(vmin=1e-20,
                                          vmax=pseudo_rate_map[cond].max()
                                      ),
                      cmap='cubehelix'
                      )
        
        ax1.set(xscale='log', xlim=xlim, title='Reaction throughput / $\\mathrm{\\frac{M}{s}}$')
        f.colorbar(c1, ax=ax1,)
        c2 = ax2.pcolormesh(X, Y, pseudo_rate_map/pseudo_rate_map.max(axis=0),            
                          cmap='inferno'
                          )
        ax2.set(xscale='log', xlim=xlim, title='Relative reaction activity')
        f.colorbar(c2, ax=ax2)
        for ax in axes:
            ax.set_xlabel('$t$ / s')
            ax.invert_yaxis()
        f.suptitle(f'{reaction_set.split(".")[0]}')
        for end in ['png', 'svg']:
            f.savefig(f'{savename}.{end}', dpi=300, transparent=True)
        plt.show()
        plt.close('all')
    
    if save_data:
        file_type = simulation_file.split('.')[-1]
        if file_type == 'xlsx':
            reaction_throughput_df.to_excel(f'{savename}.xlsx')
        else:
            reaction_throughput_df.to_csv(f'{savename}.{file_type}')
        
    return reaction_throughput_df


if __name__ == '__main__':
    
    reaction_throughput_df= analyze_reaction_throughput('Reactions_H2O_sparsened.txt',
                                'Reactions_H2O_sparsened_C0_H2O_10000.0_Gys_Output.xlsx', 
                                save_data=True, show=True,
                                    )