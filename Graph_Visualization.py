import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import AuRaCh as au


"""
Created on Fri Nov 12 16:51:45 2021

@author: b.fritsch, nzarg
"""


def nodecolor(node, disabled=[], 
              second_node=None,
              show_only_disabled=False,
              disabled_color='gainsboro', 
              base_color='C0',
              color_dict={}):
    
    if second_node is None:
        second_node = node
    
    disabled_cond = any((node in disabled, second_node in disabled))
    if show_only_disabled:
        disabled_cond = any((node not in disabled, second_node not in disabled))
    
    if disabled_cond:
        return disabled_color
    
    for entry in color_dict:
        if entry in node:
            print('FOund an entry:', entry)
            return color_dict[entry]
   
    if 'Au' in node:
        return 'C1'
    elif 'Ag' in node:
        return 'silver'
    elif 'Cl' in node:
        return 'C2' 
    elif 'Br' in node:
        return 'tab:purple'
    elif 'Fe' in node:
        return 'tab:cyan'
    elif any(['Na' in node,
              'Cs' in node,
              'K' in node]):
        return'tab:pink'
    elif 'I' in node:
        return 'tab:olive'
    elif 'N' in node:
          return 'C3' 
    elif 'C' in node:
         return 'dimgray'
    elif 'R' in node:
        return 'tab:brown'
    else:
        return base_color
 

def draw_graph(reaction_set, disabled=[], figsize=(12,12), equal_node_size=False,
               seed=101, spring_iterations=55, spring_strength=0.9, color_dict={}
               ):
    """
    graph_draw_args go to nx.draw

    Parameters
    ----------
    reaction_set : TYPE
        DESCRIPTION.
    disabled : TYPE, optional
        DESCRIPTION. The default is [].
    figsize : TYPE, optional
        DESCRIPTION. The default is (12,12).


    Returns
    -------
    None.

    """
    
    sp = reaction_set
    
    rate_constants, reaction_map = au.create_reactions(sp)
    
    out = []
    out1 = []
    """
    calculate educt and product weight in each reaction and save as Dataframe
    
    """
    for key in reaction_map:
        reaction = reaction_map[key]
        for  spec in reaction['educts']:
            for k, prod in enumerate(reaction['products']):
                if rate_constants[key] == 0:
                    continue
                educt_weight = np.log10((rate_constants[key]) * (reaction['product prefix'][k]))
                educt_w_nolog= (rate_constants[key]) * (reaction['product prefix'][k])
                educt_list = [spec, prod, educt_weight, educt_w_nolog ,'R'+str(key)]
                out.append(educt_list)

        for spec in reaction['products']:
            for k, educt in enumerate(reaction['educts']):
                if rate_constants[key] == 0:
                    continue
                product_weight = np.log10((rate_constants[key]) * (reaction['educt prefix'][k]))
                product_w_nolog = (rate_constants[key]) * (reaction['educt prefix'][k])
                product_list = [educt, spec,product_weight, product_w_nolog, 'R'+str(key)]
                out1.append(product_list)
             
    for x in range(len(out)):
         out[x].insert(4,out1[x][2]/out[x][2])
         
    
    df = pd.DataFrame(out, columns = ['source', 'target', 'educt_weight',
                                      'educt_w_nolog','weight_div',
                                    'label'])

   
 
    df['Node'] = df['source']
    df['source'] = df['source'].apply(au.make_nice_label) 
    df['target'] = df['target'].apply(au.make_nice_label)
    
    educt_w_total_no_normalized = df.groupby('source')['educt_weight'].sum()
    
    ## normalize the weights to value in range 0 to 1 ##
    educt_max_weight = df['educt_weight'].max()
    educt_min_weight = df['educt_weight'].min()
    df['educt_normalized_weight'] = df['educt_weight']  / educt_max_weight
    df['educt_normalized_weight'] -= df['educt_normalized_weight'].min()*1.1
    df['educt_normalized_weight'] /= df['educt_normalized_weight'].max()
    df['educt_shifted_weight'] =df['educt_weight'] 
    if educt_min_weight < 0:
        df['educt_shifted_weight'] -=  educt_min_weight
        
        
    """
    Create Graph from Dataframe
    
    """
    G = nx.from_pandas_edgelist(df, source = 'source', target = 'target',
                                edge_attr='educt_normalized_weight', create_using=nx.DiGraph())
    
    #G1 = nx.from_pandas_edgelist(df1, source = 'source', target = 'target',edge_attr='product_weight', create_using=nx.DiGraph())
    node_num_connection = {}
    nodelist = list(G.nodes())
 
    list_degree= list(G.degree())
    list_degree.sort()
    nodes , degree = map(list, zip(*list_degree))
    
    
    for node in nodelist:
        node_num_connection.update( {node : G.degree(node)} )
    node_num_connection = {k: v for k, v in sorted(node_num_connection.items(), key=lambda item: item[1])}
     

    educt_node_total_weight = df.groupby('source')['educt_normalized_weight'].sum()
    educt_w_nolog_nonormalized = df.groupby('source')['educt_w_nolog'].prod()
    w_div =  df.groupby('source')['weight_div'].sum()
    educt_node_weight_dct = {}
  
    for node in nodelist:
         for n in df['Node']:
             if node == au.make_nice_label(n):
                 educt_node_weight_dct.update({(node,n) : [educt_node_total_weight[node] ,educt_w_nolog_nonormalized[node], educt_w_total_no_normalized[node],w_div[node]]})        
                                  
    figsize = (12,12)
    

    plt.figure(figsize=figsize, layout='constrained')
    plt.clf()
    pos = nx.spring_layout(G, seed=seed, iterations=spring_iterations, k=spring_strength)

    weight_node = np.array([w[0] for node, w in educt_node_weight_dct.items()])
    w_min = weight_node.min()
    w_max = weight_node.max()
    weight_node = (weight_node - w_min+1 ) / ( w_max -  w_min)
    node_size = np.array([nx.betweenness_centrality(G,
                    weight='educt normalized weight',
                    normalized=True)[n] for n in nodes])*5000
    if equal_node_size:
        node_size = np.full(node_size.shape, 500)        
    fac=1.5
    disabled_draw = [au.make_nice_label(node) for node in disabled]
    
    nx.draw(G, 
            node_color = [nodecolor(node, disabled_draw, show_only_disabled=False, color_dict=color_dict)  for node in nodes],
            with_labels = True, 
            edge_color = [nodecolor(u, disabled_draw, v, show_only_disabled=False, color_dict=color_dict) for u , v in G.edges()],
            pos = pos,
            nodelist = nodes,
            node_size = node_size,
            connectionstyle='arc3, rad = 0.3',
            alpha=0.8,
            width = np.array([w[list(w.keys())[0]] for __, __, w in G.edges(data=True)])*fac
            )
    
    
    disabled_str = 'all species'
    if disabled:
        disabled_str = 'disabled ' + ' '.join(disabled)
    
    for end in ['pdf', 'svg', 'png']:
        plt.savefig(f'Graph {sp} {disabled_str}.{end}', dpi=600, transparent=True) 
    plt.show()
    plt.close('all')




if __name__ == "__main__":
    
    
    reaction_set = 'Reaction_Database.txt'
    
    disabled = ['HCl', 'H+', 'OH-', 'Cl', 'Au(III)Cl4-', 'ClOH-', 'Au(I)Cl2-',
                'Cl2-', 'eh-', 'HO2-', 'H', 'Au2Cl6^2-', 'OH', 'O2-', 'HO2',
                'O2', 'H2', 'Au', 'Cl-', 'Au(III)Cl4-', 'H2O2', 'H2O']
    
    draw_graph(reaction_set, disabled=[], figsize=(12,12), equal_node_size=True,
                   seed=101, spring_iterations=55, spring_strength=0.9,
                   )