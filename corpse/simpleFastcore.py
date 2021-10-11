# Porthmeus
# 08.10.21

# this is a wrapper in order to make running fastcore much more easier than it is currently implemented in troppo or matlab

import cobra as cb
import cobamp as ca
from troppo.methods.reconstruction.fastcore import FASTcore, FastcoreProperties
import warnings


class simpleFastcore():
    ''' This class is an attempt to simplify the usage of fastcore - simply supply the model which should be used and a set of core reactions as a list of reaction IDs and call self.run() - and get the ready to use core model. All steps in between can be executed as well, see self.fastcc() and self.fastcore(). Note: these are all just wrappers for functions implemented in cobrapy, cobamp and troppo'''
    def __init__(self, model, core_set, solver = None, max_boundaries = 1000, flux_threshold = 1E-4):

        self.model = model.copy()
        self.solver = solver
        self.max_boundaries = max_boundaries
        self.flux_threshold = flux_threshold
        
        if all([True if type(x) == int else False for x in core_set]):
            self.core_idx = core_set
            self.core_set = [self.model.reactions[i].id for i in core_set]
        elif all([True if type(x) == str else False for x in core_set]):
            self.core_set = core_set
            self.core_idx = [i for i,x in enumerate(self.model.reactions) if x.id in self.core_set]
        else:
            raise ValueError("core_set needs to be either index of name_ID of reactions in the model")

        self.check_solver()
        self.check_boundaries()
        self.status = ["initiated"]
        
    
    def check_solver(self):
        if self.solver == None:
            self.solver = str(self.model.solver.interface.Model).split(".")[1].replace("_interface","").upper()
        print("Will use " + self.solver + " for linear problems")
    
    def check_boundaries(self):
        # check the solver which is selected and adjust infinity boundaries if necessary
        if self.solver == None:
            self.check_solver()
            self.check_boundaries()

        print("Will adjust model boundaries which are set to inf/-inf and set it to " +str(self.max_boundaries))
        inf = float("inf")
        for rxn in self.model.reactions:
            if rxn.upper_bound == inf:
                rxn.upper_bound = self.max_boundaries
            elif rxn.upper_bound == -1*inf:
                rxn.upper_bound = -1*self.max_boundaries
            if rxn.lower_bound == inf:
                rxn.lower_bound = self.max_boundaries
            elif rxn.lower_bound == -1*inf:
                rxn.lower_bound = -1*self.max_boundaries

        self.status = ["boundaries checked"]
    
    def fastcc(self):
        # this function wil return the consistent model
        # make a consistent model
        print("Creating consistent model")
        init_rxn = len(self.model.reactions)
        self.model = cb.flux_analysis.fastcc(model = self.model, flux_threshold = self.flux_threshold)
        print("Reduced the initial model from {init} to {final} consistent model".format(init = str(init_rxn) , final = str(len(self.model.reactions)))) 
        # get the cobamp model
        #self.ca_model = ca.wrappers.COBRAModelObjectReader(self.model)    
        # map the core set to the ids of the cobamp model
        self.core_idx = [i for i,x in enumerate(self.model.reactions) if x.id in self.core_set]
        if len(self.core_idx) == 0:
            raise ValueError("No core set left in the consistent model, check your input!")
        
        self.status.append("fastcc")

        return(self.model)
        

    def fastcore(self):
        # this function will return the specific model, make sure to run all preparation steps beforhand
        # initiate the fastcore model 
        S = cb.util.create_stoichiometric_matrix(self.model)
        lb = [x.lower_bound for x in self.model.reactions]
        ub = [x.upper_bound for x in self.model.reactions]
        #lb,ub = self.ca_model.get_model_bounds(False, True)
        fastcoresolver = FASTcore(S, lb, ub,
                FastcoreProperties(core= self.core_idx,
                    solver = self.solver,
                    flux_threshold = self.flux_threshold)
                )

        # run fastcore
        specific_idx = fastcoresolver.fastcore()
        self.specific_idx_caMod = specific_idx
        print("Fastcore is done, will adjust the model")

        # get the specific model
        #specific_IDs = [x for i,x in enumerate(self.ca_model.r_ids) for i in specific_idx]
        #self.specific_IDs = specific_IDs
        rm_rxns = [x for i,x in enumerate(self.model.reactions) if i not in specific_idx]
        self.model.remove_reactions(rm_rxns)
        
        if len(specific_idx) != len(self.model.reactions):
            warnings.warn("Something odd happened and the reactions extracted by fastcore and the ones which remained in the model have not the same length. Please check the results carefully!")

        self.status.append("fastcore")

        return(self.model)

    def get_model(self):

        return(self.model.copy())

    def run(self):

        self.fastcc()
        self.fastcore()
        return(self.model.copy())
        
