# CORPSE - CObRaPy tool SEt
CORPSE is an extension to the CobraPy environment to make some extra functionality - I initiated this because I missed functions for reconstruction of specific models.

The toolset includes functions for:
+ mapping gene expression values to reactions
+ get essential genes/rxn depending on the global lower and upper and local thresholds
+ calculate distances on FVA results for different samples/simulations of the same model 
+ extract tissue specific models using fastcore

## Installation

Simply use pip to install the package:

    pip install git+https://github.com/Porthmeus/CORPSE.git

Alternatively clone/download the package and install it manually

    git clone https://github.com/Porthmeus/CORPSE.git
    cd CORPSE
    python setup.py install


## Usage

### Fastcore

Here a simple example of how to use fastcore. First we load the E. coli model from cobrapy as a toy model and we take 5 random reactions as a core set.

    import cobra
    from corpse import simpleFastcore
    import random
    
    eco = cobra.io.load_model("textbook") # load the model
    core_eco = random.sample([x.id for x in eco.reactions], 5) # sample 5 reactions
    
To run fastcore instantiate a fastcore object and call run():
    
    fast_mod = simpleFastcore(model = eco, core_set = core_eco)
    core_mod = fast_mod.run()

This will run fastcc and fastcore on the model and return the tissue specific model. One can run also only fastcc or fastcore on the model provided:
    
    # run only fastcc
    fast_mod.fastcc()
    consist_mod = fast_mod.get_model()

    # run only fastcore
    fast_mod.fastcore()
    core_mod = fast_mod.get_model()

To change the solver for the problems, simply change the solver of the initial model object provided:

    eco.solver = "cplex"
    fast_mod = simpleFastcore(model = eco, core_set = core_eco)

### Omics Mapper

To define the active reactions, usually one uses omics data to map these to the reactions. `corpse` offers to map transcriptomic or proteomic data (or any other measure which can be associated to genes/proteins) to the model in questions. Given you have transcriptomic data where the rows contain genes and the columns samples one can use the following code to map.

```
import corpse
import cobra
import pandas

mod = cobra.io.read_sbml_model("path/to/model.xml")
df = pandas.read_csv("path/to/csv", index_col = 0)
mapper = corpse.omicsMapper()
RAS_df = mapper.mapExpressionToReaction(model = mod, dataframe = df)
```
Set `protein = True` if you want to use the gene names instead of the IDs of the model.

## TODO:
Write a good documentation with examples how to use the library
