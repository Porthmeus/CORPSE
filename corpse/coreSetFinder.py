# Porthmeus
# 18.08.21

import warnings
import pandas as pd
import numpy as np
import re

class coreSetFinder:
    def __init__(self):
        self.name="coreSetFinder"

    def getCoreSet(
            array,
            global_lower = 0,
            global_upper = None,
            local = None,
            subset = None):
        '''calculate the genes/rxns which are considered expressed in a sample
        @array = pandas.DataFrame - contains the expression/activity values for each gene/rxn (row) and sample (column) - should contain only expression values and gene names as index and sample names as columnnames.
        @global_lower = float[0-100] - defines the global lower threshold as percentile from the whole data set, genes/rxn which have an expression/activity lower than this threshold are considered inactive.
        @global_upper = float[0-100, > global_lower] - defines the global upper threshold for genes/rxns expression/activity as percentile of the whole data set, genes/rxns which have an expression/activity higher than this value are considered always active. If = None, this threshold will be not employed (genes are considered active either depending on the local and/or on the global_lower threshold). If the global_upper threshold is employed, the local threshold can not = None  -> will be automatically set to 50.
        @local = float[0-100, global_lower < local < global_upper] - defines a local expression/activity threshold as percentile which is calculated individually for each gene/rxn across the data set. Genes which have an higher expression will be considered active, if expression > global_lower, conversely genes/rxn with expression/activity < local will be considered inactive, if expression/activity < global_upper.
        @subset = list - defines either the indeces or index names of the array rows to subset the data set to only those genes/rxns in the list (to subset only metabolic active genes for example)

        Value: a pandas.DataFrame of the same dimension than the input array (or the subset of it) with containing 0 if the gene/rxn in that sample is inactive or 1 if its active. Additionally a string containing a summary of thresholds applied.
        '''
        
        # sanity check
        if global_lower == None or global_lower < 0 or global_lower > 100:
            raise ValueError("Global lower threshold must be an integer between 0 and 100")

        if global_upper != None:
            if global_upper < global_lower:
                warnings.warn("Global lower threshold is higher than global upper threshold - will use only global lower threshold and discard local and global upper thresholds")
                global_upper = None
                local = None
            elif global_upper > 100:
                warnings.warn("Global upper threshold must be <= 100 - will not use global upper threshold")
                global_upper = None

        if local != None:
            if global_lower > local:
                warnings.warn("Local threshold is lower than the global lower threshold - will use only the global lower threshold")
                local = None
            elif global_upper != None and global_upper < local:
                warnings.warn("Global upper threshold is larger than the local threshold - will use only the local and global lower threshold")
                global_upper = None
            elif local > 100:
                warnings.warn("Local threshold must be <= 100 - will not use local threshold")
                local = None
                
            
        
        # subset data if necessary
        if subset != None:
            if type(subset[0]) == str:
                subset = [x for x in subset if x in array.index]
                array = array.loc[subset]
            elif type(subset[0]) == int:
                subset = [x for x in subset if x < array.shape[0]]
                array = array.iloc[subset]
            elif type(subset[0]) == bool:
                subset = subset[0:array.shape[0]]
                array = array[subset]

        # test the lower threshold
        glt = np.percentile(np.array(array), global_lower)
        resDF = array > glt   
        
        # if there is a local treshold, test for it
        if local != None:
            local_string = "L"+re.sub("\.0$","",str(local))
            for gene in array.index:
                localt = np.percentile(array.loc[gene], local)
                colVals = [True if all([x,y]) else False for x,y in zip(array.loc[gene] > localt,resDF.loc[gene])]
                resDF.loc[gene,:] = colVals
        else:
            local_string = ""

        # if there is a global upper threshold, test for it
        if global_upper != None:
            gu_string = "GU" + re.sub("\.0$","",str(global_upper))
            gut = np.percentile(np.array(array), global_upper)
            upperArray = array > gut
            resDF[upperArray] = True
        else:
            gu_string = ""
        
        # create the string
        out_string = re.sub("\|+$","","|".join(["GL" +  re.sub("\.0$","",str(global_lower)) ,local_string, gu_string]))
        
        return resDF.astype(int), out_string
