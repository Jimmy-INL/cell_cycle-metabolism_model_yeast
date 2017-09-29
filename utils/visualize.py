import pandas as pd

def interactions(interaction_objects,cobra_m,flux_distr_df,flux_distr_path = False, sim_sol_path = False ):
    """
    Create pandas dataframe holding flux values for rxns that participate in interaction
    flux_distr_df can be created from simulation solution file
    Input:
        List of interaction objects
        Cobra model
        Flux distribution dataframe, or path to simulation solution and flux distributions.
    Outputs a multi-index dataframe.
        Top-layer: CDC_node-met_enzyme pair. 
        Bottom layer: interacting reactions
    """
    
    if flux_distr_path and sim_sol_path:
        #to be created
        pass
        
    elif flux_distr_path and sim_sol_path:
        raise NameError, "either flux_distr_path: %r or sim_sol_path: %r is False. Define both or none" %(flux_distr_path,sim_sol_path)
    
    else:
        multi_index = [ [],
                        [], ]
        interaction_data = [ [] for state in flux_distr_df.columns ]
        for intx in interaction_objects:
            #get list of rxn id's and list of rxn names
            rxn_lst = [rxn for rxn in intx.get_reaction_list(cobra_m)]
            rxn_ids = [rxn.id for rxn in rxn_lst]
            rxn_names = [rxn.name for rxn in rxn_lst]
            intx_name = intx.name
            for rxn_name in rxn_names:
                multi_index[0].append( intx_name )
                multi_index[1].append( rxn_name )
            for state in flux_distr_df.columns:
                for rxn_id in rxn_ids:
                    flux = flux_distr_df[state][rxn_id]
                    interaction_data[state].append(flux)
        interactions_df = pd.DataFrame.from_records(data = interaction_data, index = flux_distr_df.columns, columns = multi_index)
    return interactions_df

def exchanges(cobra_m,flux_distr_df,flux_distr_path = False, sim_sol_path = False ):
    """
    Create pandas dataframe holding flux values for non_zero exchange reactions
    flux_distr_df can be created from simulation solution file
    Input:
        Cobra model
        Flux distribution dataframe, or path to simulation solution and flux distributions.
    Outputs a multi-index dataframe.
        Top-layer: CDC_node-met_enzyme pair. 
        Bottom layer: interacting reactions
    """
    
    
    if flux_distr_path and sim_sol_path:
        #to be created
        pass
        
    elif flux_distr_path and sim_sol_path:
        raise NameError, "either flux_distr_path: %r or sim_sol_path: %r is False. Define both or none" %(flux_distr_path,sim_sol_path)
    
    #Create a list of rxn id's for boundary rxns that cary flux for at least one state
    exchange_rxns = [rxn.id for rxn in cobra_m.exchanges if flux_distr_df.loc[rxn.id].any()]
    #Use exchange_rxns as mask to create df only containing non-zero exchanges
    exchange_df = flux_distr_df.loc[exchange_rxns]
    #Change rxn id's for full names
    exchange_df.index = [ cobra_m.reactions.get_by_id(rxn).name for rxn in exchange_rxns]
    return exchange_df.transpose()


def KEGG_trends(cobra_m,flux_distr_df,KEGG_dct,flux_distr_path = False, sim_sol_path = False):
    """
    Input: 
        KEGG_dct: A dct containing lists with KEGG pathway identifiers. Key = name of collection of paths
        cobra_m: constraint-based model in which reactions are annotated with KEGG 'PATH'
        flux_distr: flux_distr_df can be created from simulation solution file (to write)
    To write: compute amount of overlapping interactions between trends
    Output:
        Dataframe with cumulative flux for set of pathways for every state
    """
    #assert whether model is annotated with KEGG pathways
    for rxn in cobra_m.reactions:
        assert 'PATH' in rxn.annotation, "Annotation Error, rxn %r not annotated with 'PATH' " %rxn.id
    
    column_names = []
    KEGG_trend_data = []
    for trend_name,path_id_lst in KEGG_dct.iteritems():
        
        #Retrieve list of rxns that are annotated with 1 or more paths from path_id_lst
        rxn_ids = [ rxn.id for rxn in cobra_m.reactions if not set(cobra_m.reactions.get_by_id(rxn.id).annotation['PATH']).isdisjoint(path_id_lst)]
        #Sum fluxes of all path_id annotated rxns for every boolean state to create series
        trend_cum_flux_s = flux_distr_df.loc[rxn_ids].sum(axis=0)
        KEGG_trend_data.append(trend_cum_flux_s)
        column_names.append(trend_name)
    
    return pd.concat(KEGG_trend_data,axis=1,keys=column_names)