# Porthmeus
# 25.08.21

import troppo
import time

model = cb.io.read_sbml_model("/home/taube/Work/miTarget/Recon2/COLORMORE22/colormore22/colormore22.xml")
dataframe = pd.read_csv("/home/taube/Work/miTarget/eMed_FUTURE/Pipeline_MM/results/data/TPM_emed_future.csv", index_col =0)
ix = troppo.omics.id_converter.idConverter(dataframe.index, old = "ensembl_gene_id", new = "hgnc_id")
ix2 = [ix[x] if x in list(ix.keys()) else x for x in dataframe.index ]
dataframe.index = ix2

protein = False
orIsSum = True
column = dataframe.columns[0]

gpr = model.reactions[15].gene_reaction_rule
gpr2 = model.reactions[0].gene_reaction_rule

tic = time.time()
mapExpressionToReaction(model, dataframe.iloc[:,1:20])
toc = time.time()
print(toc - tic)

import multiprocessing
from joblib import Parallel, delayed

def myFunc(lst):
    newlst = []
    for i in lst:
        if i%2 == 0:
            newlst.append(i*i)
    return(newlst)

def myFunc2(i, add):
    if i%2==0:
        return(i*i+add)

d = 9
joblib.Parallel(n_jobs=num_cores)(joblib.delayed(myFunc2)(i = i, add =d) for i in range(11))


colnames = list(dataframe.columns[1:10])
num_cores = 3
results = joblib.Parallel(n_jobs = num_cores)(joblib.delayed(mapSampleToModel)(model = model,
    dataframe = dataframe,
    column = sample,
    protein = protein,
    orIsSum = orIsSum) for sample in colnames)


e = {"a" : 1, "b" : 2, "c" :3 ,"d":4, "e":5}
gpr = "(a and b) or (c and d) or (a and d) or (c and b)"
gpr = "( a and b ) or ( c and d ) or ( a and d ) or ( c and b )"
gpr = "e and (d and (b or a))"
