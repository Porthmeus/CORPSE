# CORPSE - CObRaPy tool SEt
CORPSE is an extension to the CobraPy environment to make some extra functionality - I initiated this because I missed functions for reconstruction of specific models.

The toolset incudes functions for:
+ mapping gene expression values to reactions
+ get essential genes/rxn depending on the global lower and upper and local thresholds
+ calculate distances on FVA results for different samples/simulations of the same model 

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
    import cobra.test
    import corpse
    import random
    
    eco = cobra.test.create_test_model("textbook") # load the model
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

## TODO:
Write a good documentation with examples how to use the library
