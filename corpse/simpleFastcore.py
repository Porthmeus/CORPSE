# Porthmeus
# 08.10.21

# this is a wrapper in order to make running fastcore much more easier than it is currently implemented in troppo or matlab

import cobra as cb
import cobamp as ca
from troppo.methods.reconstruction.fastcore import FASTcore, FastcoreProperties
import warnings
import time
from contextlib import contextmanager,redirect_stderr,redirect_stdout
from os import devnull


class simpleFastcore():
    ''' This class is an attempt to simplify the usage of fastcore - simply supply the model which should be used and a set of core reactions as a list of reaction IDs and call self.run() - and get the ready to use core model. All steps in between can be executed individually as well, see self.FVA_consistency(), self.fastcc(), self.fastcc_repeat() and self.fastcore(). Note: these are all just wrappers for functions implemented in cobrapy, cobamp and troppo'''
    def __init__(self, model, core_set = [],  max_boundaries = 1000, zero_cutoff = None):

        self.model = model.copy()
        self.model_ori = model.copy()
        self.status = []
        self.solver = None
        self.max_boundaries = max_boundaries
        if zero_cutoff == None:
            self.zero_cutoff = self.model.tolerance
        else:
            self.zero_cutoff = zero_cutoff

        
        self.set_core_set(core_set = core_set)
        self.core_set_ori = core_set
        self.check_solver()
        self.check_boundaries()
        self.fastcc_runs = 0
        

    @contextmanager
    def suppress_stdout_stderr(self):
        """A context manager that redirects stdout and stderr to devnull"""
        # this is just needed to catch the chatter of the troppo functions
        with open(devnull, 'w') as fnull:
            with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
                yield (err, out)
        
    
    def check_solver(self):
        ''' Gets the solver which is set for the model'''
        if self.solver == None:
            self.solver = str(self.model.solver.interface.Model).split(".")[1].replace("_interface","").upper()
        print("# Will use " + self.solver + " as linear problem solver")
        self.status.append("Solver_checked")
    
    def check_boundaries(self):
        ''' Checks the model reactions and sets the boundaries to max_boundaries (see __init__()) if the value is currently set to Inf/-Inf'''
        # check the solver which is selected and adjust infinity boundaries if necessary
        if self.solver == None:
            self.check_solver()
            self.check_boundaries()

        print("# Will adjust model boundaries which are set to inf/-inf and set it to {boundary}/-{boundary}".format(boundary = str(self.max_boundaries)))
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

        self.status.append("boundaries_checked")

    def check_core_rxns(self):
        ''' Checks if the core set is not empty after extracting the consistent model'''
        self.core_idx = [i for i,x in enumerate(self.model.reactions) if x.id in self.core_set]
        if len(self.core_idx) == 0:
            raise ValueError("No core set left in the model provided, check your input!")
    
    def reset_model(self):
        ''' resets the model to the original input model in case one wants to use a different consistency algorithm'''
        # reset the model to the initial state including the checks
        self.model = self.model_ori.copy()
        self.status =[]
        self.check_solver()
        self.check_boundaries()
        self.set_core_set(self.core_set_ori)

        
    def set_core_set(self, core_set):
        ''' Set the core reactions in the model 
        @ core_set - list of reaction IDs or integers corresponding to the index of the reactions in the model'''

        if all([True if type(x) == int else False for x in core_set]):
            self.core_idx = core_set
            self.core_set = [self.model.reactions[i].id for i in core_set]
        elif all([True if type(x) == str else False for x in core_set]):
            self.core_set = core_set
            self.core_idx = [i for i,x in enumerate(self.model.reactions) if x.id in self.core_set]
        else:
            raise ValueError("core_set needs to be either index of name_ID of reactions in the model")
        
    
    def fastcc(self):
        ''' Creates a consistent model using the fastcc algorithm '''

        print("# Creating consistent model using fastcc")
        # supress a warning which is usually thrown by cobrapy
        warnings.filterwarnings("ignore", ".*need to pass in a list.*")

        # get the number of reactions to report afterwards
        tictic = time.perf_counter()
        init_rxn = len(self.model.reactions) 

        # reduce the model
        self.model = cb.flux_analysis.fastcc(model = self.model, zero_cutoff = self.zero_cutoff)
        cons_rxn = len(self.model.reactions)
        toc = time.perf_counter()
        
        
        # internally count how often fastcc was run
        self.fastcc_runs = self.fastcc_runs +1

        print("# Reduced the initial model from {init} to {final} reactions ({frac}%) in the consistent model in {tictoc}s in the first fastcc run".format(init = str(init_rxn),
            final = str(cons_rxn),
            frac = str(round(1-cons_rxn/init_rxn,3)*100),
            tictoc = str(round(toc-tictic,3)))) 

        # add the status to the status list
        if "consistent" not in self.status: 
            self.status.append("consistent")

    def FVA_consistency(self):
        ''' Creates a consistent model using FVA - this is the default as it seems to produce most consistent results'''
        # use FVA to get a consistent model
        
        print("# Creating consistent model using FVA")
        init_rxn = len(self.model.reactions)
        tic = time.perf_counter()
        blocked_rxns = cb.flux_analysis.find_blocked_reactions(self.model, zero_cutoff = self.zero_cutoff)
        self.model.remove_reactions(blocked_rxns)
        toc = time.perf_counter()
        cons_rxn = len(self.model.reactions)

        print("# Reduced the initial model from {init} to {final} reactions ({frac}%) in the consistent model in {tictoc}s in the first fastcc run".format(init = str(init_rxn),
            final = str(cons_rxn),
            frac = str(round(1-cons_rxn/init_rxn,3)*100),
            tictoc = str(round(toc-tic,3)),
            no = self.fastcc_runs)) 

        # add the status to the status list
        if "consistent" not in self.status: 
            self.status.append("consistent")

    def fastcc_repeat(self):
        '''This is a last resort attempt to create consistent model of the simpleFastcore object. It will run fastcc over and over until no inconsistent reactions can be found anymore. Use this with care, as this usually removes reactions which are actually not inconsistent and produces a different model for the context specific extraction process as basis''' 
        # this function wil return the consistent model
        
        # supress a warning which is usually thrown by cobrapy
        warnings.filterwarnings("ignore", ".*need to pass in a list.*")

        
        # get the number of reactions to initialize a loop and record some statistics
        tictic = time.perf_counter()
        init_rxn = len(self.model.reactions) 
        cons_rxn = init_rxn +1
        current_rxns = init_rxn
        print("# Creating consistent model using repeating fastcc")
        i = self.fastcc_runs
        print("Iteration\tNoInputRxns\tNoOutputRxn\tFracReduced\tTimeNeeded")

        while cons_rxn != current_rxns:
            current_rxns = cons_rxn
            i = i+1
            tic = time.perf_counter()
            self.model = cb.flux_analysis.fastcc(model = self.model, zero_cutoff = self.zero_cutoff)
            toc = time.perf_counter()

            cons_rxn = len(self.model.reactions)
            print(str(i)+ "\t" +
                    str(current_rxns) + "\t" +
                    str(cons_rxn) + "\t" + 
                    str(round(1-cons_rxn/current_rxns,3)) + "\t" +
                    str(round(toc-tic,2)))

        print("# Reduced the initial model from {init} to {final} reactions ({frac}%) in the consistent model in {tictoc}s with {no} fastcc runs".format(init = str(init_rxn),
            final = str(cons_rxn),
            frac = str(round(1-cons_rxn/init_rxn,3)*100),
            tictoc = str(round(toc-tictic,3)),
            no = i)) 

        # add the status to the status list
        if "consistent" not in self.status: 
            self.status.append("consistent")

        

    def fastcore(self):
        ''' Runs fastcore on the model and core_set stored in the simpleFastcore object'''

        # this function will return the specific model, make sure to run all preparation steps beforhand

        # map the core set to the ids of the model
        self.check_core_rxns()

        tic = time.perf_counter()
        # initiate the fastcore model extractor 
        S = cb.util.create_stoichiometric_matrix(self.model)
        lb = [x.lower_bound for x in self.model.reactions]
        ub = [x.upper_bound for x in self.model.reactions]
        with self.suppress_stdout_stderr(): # remove the noise from the troppo fastcore implementation
            fastcoresolver = FASTcore(S, lb, ub, # this is everything which is needed from the model
                    FastcoreProperties(core= self.core_idx,
                        solver = self.solver,
                        flux_threshold = self.zero_cutoff)
                    )

            # run fastcore
            specific_idx = fastcoresolver.fastcore()

        self.specific_idx_caMod = specific_idx
        toc = time.perf_counter()
        print("# Fastcore was done in {tictoc}s, will adjust the model".format(tictoc = str(round(toc-tic,3))))

        # get the specific model
        rm_rxns = [x for i,x in enumerate(self.model.reactions) if i not in specific_idx]
        self.model.remove_reactions(rm_rxns)
        
        if len(specific_idx) != len(self.model.reactions):
            warnings.warn("Something odd happened and the reactions extracted by fastcore and the ones which remained in the model have not the same length. Please check the results carefully!")

        self.status.append("context_specific")

    def get_model(self):
        ''' Returns a copy of the model'''

        return(self.model.copy())

    def run(self):
        ''' Runs the complete fastcore pipeline - makes the model consistent, and extracts the context specific model with fastcore - returns a copy of the context specific model'''
        tic = time.perf_counter()
        
        try:
            # try to run fastcore with FVA - that seems more stable
            self.FVA_consistency()
            self.fastcore()
        except ValueError as VE:
            raise VE
        except Exception as exception:
            # if that does not work try to run fastcore with the fastcc consistent model
            print(exception)
            warnings.warn("WARNING: Tried to run fastcore with FVA which failed, will try to run it with fastcc")
            try:
                # try the FVA approach
                self.reset_model()
                self.fastcc()
                self.fastcore()
            except ValueError as VE:
                raise VE
            except Exception as exception:
                print(exception)
                # if that does not work either, try the repeated fastcc
                warnings.warn("WARNING: Tried to run fastcore with fastcc and FVA which failed, will try to run it with repeated fastcc")
                self.reset_model()
                self.fastcc_repeat()
                self.fastcore()
                warnings.warn("WARNING: Using repeated fastcc might remove consistent reactions from the initial model - check your results carefully!")

        toc = time.perf_counter()
        print("# Total runtime: " + str(round(toc-tic,3)) + "s")
        return(self.model.copy())
        
