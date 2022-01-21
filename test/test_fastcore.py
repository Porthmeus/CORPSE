# Porthmeus
# 8.10.21

import time
import cobra.test
import random
import pandas as pd

%run corpse/simpleFastcore.py

# simple test
eco = cobra.test.create_test_model("textbook") # load the model

fast_mod = simpleFastcore(model = eco, core_set = core_eco)
fast_mod.run()
fast_mod.fastcc()
core_eco = random.sample([x.id for x in fast_mod.model.reactions], 5) # sample 5 reactions
fast_mod.set_core_set(core_eco)

fast_mod.reset_model()
fast_mod.FVA_consistency()
fast_mod.fastcore()

fast_mod.reset_model()
fast_mod.fastcc_repeat()
fast_mod.fastcore()

core_mod = fast_mod.get_model()


eco.solver = "glpk"
fast_mod = simpleFastcore(model = eco, core_set = core_eco)
core_mod = fast_mod.run()


fast_mod = simpleFastcore(model = eco, core_set = core_eco)
fast_mod.fastcc()
fast_mod.fastcore()

# testing 
mod = cb.io.read_sbml_model("/home/taube/Work/miTarget/Recon2/COLORMORE22/colormore22/colormore22.xml")
core_set = pd.read_csv("/home/taube/Work/miTarget/eMed_FUTURE/Pipeline_MM/temp/SPLITGL25|L50_colormore22-F02234_L1_S1_L001.csv", index_col = 0)

diet = pd.read_csv("/home/taube/Work/miTarget/eMed_FUTURE/Pipeline_MM/resources/diets/colormore22_MatjesAbsorption.csv", index_col = 0)

dietRxn = [mod.reactions.get_by_id(x) for x in diet.index if x in [y.id for y in mod.reactions]]
for rxn in dietRxn:
    if rxn in diet.index:
        rxn.lower_bound = -1*diet.loc[rxn.id,:]
    else:
        rxn.lower_bound = 0

core_set = [x for x in core_set.index if core_set.loc[x][0] == 1]

mod_cplex = simpleFastcore(model = mod, core_set = core_set)
mod1 =mod_cplex.run()

mod.solver = "glpk"
mod_glpk = simpleFastcore(mod, core_set = core_set)
mod2 = mod_glpk.run()
mod_glpk.reset_model()
mod3 = mod_glpk.run()


mod_glpk.reset_model()
mod_glpk.FVA_consistency()
mod_glpk.fastcore()
mod4 = mod_glpk.get_model()
mod_glpk.reset_model()
mod_glpk.FVA_consistency()
mod_glpk.fastcore()
mod5 = mod_glpk.get_model()

eco = cobra.test.create_test_model("textbook")
core_eco = random.sample([x.id for x in eco.reactions], 5)


fast_eco = simpleFastcore(eco, core_eco)
fast_eco.run()

eco.solver = "glpk"
fast_eco = simpleFastcore(eco, core_eco)
fast_eco.run()


