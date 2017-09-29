from IPython.display import display
#import utils.instruction_lisa
#from utils import analysis_lisa
import pandas as pd
import matplotlib.pyplot as plt
import json
pd.options.display.max_columns = 166
#from bioservices import KEGG
#import numpy as np
from scripts import openfl


def instruction(m,interaction_array,VoI):
    instruction = {}
    array = utils.instruction_lisa.create_instruction(interaction_array,VoI)
    for k,v in array.iteritems():
        met_rxns = [ m.reactions.get_by_id(rxn).name for rxn in v[1] ] #Get names of interacting met. reactions
        instruction[k] = [ v[0], met_rxns , v[2]+v[3] ] #Input v[0]: CDC-nodes, met. rxns, v[2]+v[3]: direction + nature of interaction
        
    instruction_df = pd.DataFrame.from_dict(instruction)
    display(instruction_df)
    return instruction_df

def open_solution(quinstring,path = "output_files/solutions/"):
    """returns solution dict when quinstring is provided"""
    with open(path+quinstring+".json","rb") as f:
        string = f.read()
        solution = json.loads(string)
    if type(solution) != dict:
        solution =  eval(solution)
    return solution

def interactions(m,solution,interaction_array,flux_distr_path = "input_files/flux_distributions/",iterations = 14):
    """"""
    interacting_rxns = [rxn for field in interaction_array.values() for rxn in field[1] ] #List of rxns that can export/import over system boundary
    print len(interacting_rxns), "nr of interacting reactions"
    interact_dct = {m.reactions.get_by_id(rxn).name:{} for rxn in interacting_rxns} #initialize exchange dct with abbr names
    for i in range(iterations):
        ternary = solution['ternary'][i]
        with open(flux_distr_path+ternary+".json","rb") as f:
            string = f.read()
            flux_distr = json.loads(string)
        assert flux_distr != False #Solution not allowed due to death (sometimes) or non-feasible (always)
        for rxn in interacting_rxns:
            name = m.reactions.get_by_id(rxn).name 
            interact_dct[name][i] = flux_distr[rxn]

    #Create and display df
    interact_df = pd.DataFrame.from_dict(interact_dct)
    display(interact_df)
    return interact_df

def bool_ternary(solution):
    """Show boolean states and accompanying ternarystrings, return pandas dataframe"""
    if len( solution['ternary'] ) !=  len( solution['bool']['Clb1_2'] ):
        print "length ternary / bool_ternare are different"
        print solution['bool']
        print solution['ternary']
    bool_dct = solution['bool']
    bool_dct['ternary'] = solution['ternary']
    bool_df = pd.DataFrame.from_dict(bool_dct)
    display(bool_df)
    return bool_df
   
def exchange(m,solution,flux_distr_path = "input_files/flux_distributions/",iterations = 14):
    """For every state display the non-zero exchange reactions"""
    exchange_rxn = [rxn for rxn in m.exchanges] #List of rxns that can export/import over system boundary
    exchange_dct = {rxn.name:{} for rxn in exchange_rxn} #initialize exchange dct
    for i in range(iterations):
        ternary = solution['ternary'][i]
        with open(flux_distr_path+ternary+".json","rb") as f:
            string = f.read()
            flux_distr = json.loads(string)
        assert flux_distr != False #Solution not allowed due to death (sometimes) or non-feasible (always)
        for rxn in exchange_rxn:
            exchange_dct[rxn.name][i] = flux_distr[rxn.id]

    #Delete exchange reactionts that are 0 for every boolean state
    del_keys = []
    for key in exchange_dct:
        save = False
        for i in exchange_dct[key]:
            if abs(exchange_dct[key][i]) > 10**-5:
                save = True
        if save == False:
            del_keys.append(key)
    for key in del_keys:
        del exchange_dct[key]   
    exchange_df = pd.DataFrame.from_dict(exchange_dct)
    display(exchange_df)
    return exchange_df

def biom_met_exchange(m,phase_objx,solution,flux_distr_path = "input_files/flux_distributions/",iterations=14):
    #Take all sink metabolites from the biomass function to visualize exchange
    biom_met_exchange_dct = {met.name:{} for met,value in m.reactions.r_4041.metabolites.iteritems() if value < 0}
    
    for i in range(iterations):
        ternary = solution['ternary'][i]
        binary = solution['binary_objx'][i]
        obj_rxn = phase_objx[binary]
        mets_produced = {met.name:coeff for met,coeff in obj_rxn.metabolites.iteritems()}
        with open(flux_distr_path+ternary+".json","rb") as f:
            string = f.read()
            flux_distr = json.loads(string)        
        assert flux_distr != False, "Visualizing solutions containing Infeasible sol is not possible"
        
        obj_value = flux_distr[obj_rxn.id]   
        for met in biom_met_exchange_dct:
            biom_met_exchange_dct[met][i] = 0 #Initialize value for every biom_met at 0
            if met in mets_produced:
                coeff = mets_produced[met]
                met_quantity = obj_value * coeff #Multiply flux at iteration with coeff in obj. function
                biom_met_exchange_dct[met][i] += met_quantity #not produced mets stay at 0
    
    met_exchange_df = pd.DataFrame.from_dict(biom_met_exchange_dct)
    display(met_exchange_df)
    return met_exchange_df

""" KEGG pathways that together comprise the three macro molecular classes of lipids, proteins and nucleotides  """
lipid_metabolism = ["Fatty acid biosynthesis","Fatty acid elongation","Fatty acid degradation",
                    "Synthesis and degradation of ketone bodies","Cutin, suberine and wax biosynthesis",
                   "Steroid biosynthesis","Primary bile acid biosynthesis","Secondary bile acid biosynthesis",
                   "Steroid hormone biosynthesis","Glycerolipid metabolism","Glycerophospholipid metabolism",
                   "Ether lipid metabolism","Sphingolipid metabolism","Arachidonic acid metabolism",
                   "Linoleic acid metabolism","alpha-Linolenic acid metabolism","Biosynthesis of unsaturated fatty acids",
                   "Lipid biosynthesis proteins"]
amino_acid_metabolism = ["Alanine, aspartate and glutamate metabolism","Glycine, serine and threonine metabolism",
                 "Cysteine and methionine metabolism","Valine, leucine and isoleucine degradation",
                 "Valine, leucine and isoleucine biosynthesis","Lysine biosynthesis","Lysine degradation",
                 "Arginine biosynthesis","Arginine and proline metabolism","Histidine metabolism",
                 "Tyrosine metabolism","Phenylalanine metabolism","Tryptophan metabolism",
                    "Phenylalanine, tyrosine and tryptophan biosynthesis","Amino acid related enzymes"]
nucleotide_metabolism = ["Purine metabolism","Pyrimidine metabolism"]
glycolysis = ["Glycolysis / Gluconeogenesis"]
Citrate_cycle = ["Citrate cycle (TCA cycle)"]

important_paths = [lipid_metabolism,amino_acid_metabolism,nucleotide_metabolism,glycolysis,Citrate_cycle]

def macromol_trends(m,ternarystrings,flux_folder,flname=False):
    important_path_rxn = {"lipid_metabolism":[],"amino_acid_metabolism":[],"nucleotide_metabolism":[],"glycolysis":[],"TCA":[]}
    
    #Iterate over all reactions in model to store in the three categories
    for rxn in m.reactions:
        try:
            if type(rxn.annotation['PATH'])==list:
                for path in rxn.annotation['PATH']:
                    if path in important_paths[0]:
                        important_path_rxn["lipid_metabolism"].append(rxn)
                    elif path in important_paths[1]:
                        important_path_rxn["amino_acid_metabolism"].append(rxn)
                    elif path in important_paths[2]:
                        important_path_rxn["nucleotide_metabolism"].append(rxn)
                    elif path in important_paths[3]:
                        important_path_rxn["glycolysis"].append(rxn)
                    elif path in important_paths[4]:
                        important_path_rxn["TCA"].append(rxn)
                    else:
                        pass
            else:
                if path in important_paths[0]:
                    important_path_rxn["lipid_metabolism"].append(rxn)
                elif path in important_paths[1]:
                    important_path_rxn["amino_acid_metabolism"].append(rxn)
                elif path in important_paths[2]:
                    important_path_rxn["nucleotide_metabolism"].append(rxn)
                elif path in important_paths[3]:
                    important_path_rxn["glycolysis"].append(rxn)
                elif path in important_paths[4]:
                    important_path_rxn["TCA"].append(rxn)
                else:
                    pass
        except KeyError: #Reaction has no annotation key 'PATH'
            pass

    #Print how many reactions in which category
    for k,v in important_path_rxn.iteritems():
        print k,len(v), "annotated reactions"

    #Analyze and print overlap between categories
    overlap1 = set(important_path_rxn["lipid_metabolism"]) & set(important_path_rxn["amino_acid_metabolism"])
    overlap2 = set(important_path_rxn["lipid_metabolism"]) & set(important_path_rxn["nucleotide_metabolism"])
    overlap3 = set(important_path_rxn["nucleotide_metabolism"]) & set(important_path_rxn["amino_acid_metabolism"])
    overlap4 = set(important_path_rxn["glycolysis"]) & set(important_path_rxn["TCA"])
    print "overlap lipid and amino_acid metabolism ",len(overlap1)
    print "overlap lipid and nucleotide metabolism ",len(overlap2)
    print "overlap nucleotide and amino_acid metabolism ",len(overlap3)
    print "overlap glycolysis and TCA ",len(overlap4)
    
    #Calculate flux through paths
    paths = {path:None for pathset in important_paths for path in pathset} #Get all KEGG pathways in map
    cumflux_dct,path_rxns = analysis_lisa.init_cumfluxpath_dict(m,paths) #Prepare dict to safe flux through paths
    iteration = 0
    fba_dct = {}
    first = True
    for ternarystring in ternarystrings:
        with open(flux_folder+ternarystring+".json","rb") as f:
            string = f.read()
            fluxdist = json.loads(string)
        if first == True:
            fba_dct = {rxn_id:{} for rxn_id in fluxdist.keys()}
            first = False
        for rxn_id in fba_dct:
            fba_dct[rxn_id][iteration] = fluxdist[rxn_id]
        iteration += 1
    cumflux_dct = analysis_lisa.calc_cumflux(cumflux_dct,path_rxns,fba_dct) #Fill cumflux dict with flux
    
    #Analyse flux through paths
    print "iterations is hardcoded at 14"
    lipid_flux = [0 for x in xrange(14)]
    amino_flux = [0 for x in xrange(14)]
    nucleotide_flux = [0 for x in xrange(14)]
    glycolysis_flux = [0 for x in xrange(14)]
    TCA_flux = [0 for x in xrange(14)]
    

    for path in cumflux_dct:
        for iteration in cumflux_dct[path]: #Iterate over boolean states in cum. flux dictionary & sort in 3 categories
            path_flux = cumflux_dct[path][iteration]
            if path in lipid_metabolism:
                lipid_flux[iteration] += path_flux
            elif path in amino_acid_metabolism:
                amino_flux[iteration] += path_flux
            elif path in nucleotide_metabolism:
                nucleotide_flux[iteration] += path_flux
            elif path in glycolysis:
                glycolysis_flux[iteration] += path_flux
            elif path in Citrate_cycle:
                TCA_flux[iteration] += path_flux
            else:
                print path
                raise "above printed path does not occur in any category"
    #Plot results
    plt.plot(range(14),lipid_flux,label="lipid flux")
    plt.plot(range(14),amino_flux,label="amino acid flux")
    plt.plot(range(14),nucleotide_flux,label="nucleotide flux")
    plt.plot(range(14),glycolysis_flux,label="glycolysis flux")
    plt.plot(range(14),TCA_flux,label="TCA cycle flux")
    
    plt.xticks(range(14), ("0: Cell Size",
                             "1: Start",
                             "2: G1",
                             "3: G1",
                             "4: G1",
                             "5: S",
                             "6: G2",
                             "7: M",
                             "8: M",
                             "9: M",
                             "10: M",
                             "11: M",
                             "12: G1-stat",
                             "13: G1-stat"),
                               rotation=70
                       )
    plt.xlabel("cumulative flux")
    plt.legend(bbox_to_anchor=(0.0, 1.03, 5., 20.5022), loc=3,ncol=3,borderaxespad=0.0,title="Cum. flux for glycolysis,TCA & lipid-, amino acid- and nucleotide pathways", )
    if flname:
        plt.savefig(flname,bbox_inches='tight')
    plt.show()
    
def query_kegg(keggID):
    k = KEGG()
    res = k.get(keggID)
    d = k.parse(res)
    return d    
    
def annotate_genes(genes,pathways):
    ignored_pathways = [ 'Metabolic pathways','Carbon metabolism','Biosynthesis of amino acids','Biosynthesis of secondary metabolites',
                        'Microbial metabolism in diverse environments','Biosynthesis of antibiotics' ]
    gene_not_annotated_list  = []
    gene_annotation = {}
    for gene in genes:
        keggID = "sce:"+str(gene)
        d = query_kegg(keggID)
        if type(d) == dict:
            if 'PATHWAY' in d:
                if pathways:
                    path = [str(x) for x in d['PATHWAY'].values() if x in pathways ]
                else:
                    path = [str(x) for x in d['PATHWAY'].values() if x not in ignored_pathways ]
                if type(path) == list:
                    gene_annotation[gene] = path
                elif type(path) == str:
                    raise "why is this a string?", path
                    gene.annotation['PATH'] = [path]
                    gene_annotation[gene] = [path]
            else:
                print str(gene),'has no linked pathways'
                gene_annotation[gene] = []
        else:
            "No KEGG annotation found with this gene"
            gene_not_annotated_list.append(gene.id)
            
    return gene_annotation,gene_not_annotated_list
    
def proteome(m,proteome_df,gene_annotation,VoI_lst,solution_path,flux_distr_path):
    """
    Genes are annotated with their paths so result can be filtered on KEGG pathways
    Function first recreates VoI by loading fba_df with the flux distributions.
    Finally solutions are run through the proteomics analysis to obtain:
        - Score for every enzyme
        - Flux peaks in every phase (absolute)
    
    Returns a dataframe with Gene names, enzyme conc. peaks from data, enzyme/flux peaks from VoI's, names of rxn of enzyme
    ?????Remove pathways????
    Filtering on paths can be performed by creating a list of genes with certain paths from the gene_annotation dict
    -----------------------------------------------------------------------------------------------------------------------
    Input: 
        - VoI_lst : List of VoI's to compare
        - solution_path : folder path to get solution so sequence of flux distr can be recreated
        - flux_distr_path : folder path to get flux distributions
        - proteome_df : dataframe containing data on concentration levels over the phases
    Output:
        - Multi-column-indexed Dataframe
        
    Every time we create one or more columns of data, they have to be add to the df_data list, and to the arrays list of lists
    """
    
    #Initialize arrays to contain (multi-index) headers and df_data to contain the list of data to become the df columns
    arrays = [ [],
               [], ]
    df_data = []
    
    #Fill first columns: Genes and their pathways
    genes = []
    pathways = []
    for gene,path_lst in gene_annotation.iteritems():
        genes.append(gene) #Will form index
        pathways.append(path_lst)
        
    df_data.append(pathways)
    arrays[0].append("KEGG paths")
    arrays[1].append("KEGG paths") 
    
    #Add proteomics data to dataframe
    G1_data = []
    S_data = []
    G2_M_data = []
    for gene in genes:
        G1_data.append(proteome_df.loc[gene]['G1'])
        S_data.append(proteome_df.loc[gene]['S'])
        G2_M_data.append(proteome_df.loc[gene]['G2_M'])
    for i in range(3):
        arrays[0].append("Data")
    df_data.append(G1_data)
    arrays[1].append("G1") #2-levelled column-index
    df_data.append(S_data)
    arrays[1].append("S") #2-levelled column-index
    df_data.append(G2_M_data)
    arrays[1].append("G2_M") #2-levelled column-index
    
    fba_dct = {}
    for VoI in VoI_lst:
        
        #Recreate flux distr dataframe from solution
        solution = openfl.openjson(solution_path,VoI)
        bool_data = solution["bool"]
        phase_dct,exotic_CDC = analysis_lisa.def_CDC_phase(bool_data)
        tern_lst = solution["ternary"]
        iterations = len(tern_lst)
        for i,tern in enumerate(tern_lst):
            fba_dct[i] = openfl.openjson(flux_distr_path,tern) #Load flux distribution for every boolean state
        fba_df = pd.DataFrame.from_dict(fba_dct).transpose()
        fba_dct = fba_df.to_dict()          
        
        #Analyze flux distr. dataframe on proteomics vs. flux peaks to obtain score + absolute peaks for every protein
        
        # 1. Create lists for dataframe + highest level header
        G1_m_peaks = []
        S_m_peaks = []
        G2_M_m_peaks = []
        scores = []
        for i in range(4):
            arrays[0].append(VoI+" model peaks") #Still have to insert the actual VoI string
        
        # 2. For every protein in 'genes' calculate cum. flux for all their rxns
        prot_cumflux = {}
        for gene in genes:
            rxns = [rxn.id for rxn in m.genes.get_by_id(gene).reactions]
            prot_cumflux[gene] = []
            for state in range(iterations):
                cumflux = 0
                for rxn in rxns: 
                    cumflux += fba_dct[rxn][state]
                prot_cumflux[gene].append(cumflux)            
            
            # 3. Analyze flux peaks against proteomics data and return flux peaks + score for that protein
            G1_m,S_m,G2_M_m,score = analysis_lisa.flux_per_prot(phase_dct,prot_cumflux,gene,proteome_df)
            G1_m_peaks.append(G1_m)
            S_m_peaks.append(S_m)
            G2_M_m_peaks.append(G2_M_m)
            scores.append(score)
        #Append protein data to the dataframe data and add low level headers
        df_data.append(G1_m_peaks)
        arrays[1].append("G1 peaks") #2-levelled column-index
        df_data.append(S_m_peaks)
        arrays[1].append("S peaks") #2-levelled column-index
        df_data.append(G2_M_m_peaks)
        arrays[1].append("G2_M peaks") #2-levelled column-index
        df_data.append(scores)
        arrays[1].append("score") #2-levelled column-index
        
    df = pd.DataFrame(np.column_stack(df_data),index = genes, columns = arrays)
    return df
    
                       
        
        
    
    
    


"""exception clause depracated"""    
#for rxn in exchange_rxn:
#            try:
#            except TypeError:
                #            exchange_dct[rxn.name][i] = 'death'


