# Porthmeus
# 29.09.21

# calculate the sample distance based on the differences in the FVA

import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as spc

class FVAjuggler:
    def __init__(self):
        self.name = "FVAjuggler"

    def filterFVA(self,min_mat, max_mat, method ="any"):
        # filter the FVA results by variance

        # remove numerical stuff
        min_mat[np.absolute(min_mat) < 1E-6] = 0
        max_mat[np.absolute(max_mat) < 1E-6] = 0

        min_var = np.var(min_mat, 1)
        max_var = np.var(max_mat, 1)

        
        t_max = 0 #np.quantile(max_var,thrsld)
        t_min = 0 #np.quantile(min_var,thrsld)
        print(t_max)
        print(t_min)
        if method == "any":
            sel = np.any(pd.DataFrame({"min":min_var > t_min, "max":max_var>t_max}),1)
        elif method == "all":
            sel = np.all(pd.DataFrame({"min":min_var > t_min, "max":max_var>t_max}),1)
        else:
            raise ValueError("Method must be 'all' or 'any'")

        min_mat = min_mat.loc[sel,:]
        max_mat = max_mat.loc[sel,:]
        return(min_mat, max_mat)

    def clusterFVA(self, min_mat, max_mat, thrld = 0.75):
        # clusters the data by correlation and returns the most representative data point for each cluster
        # thrld - defines the minimum correlation within one cluster, should range between 0 and 1
        
        if thrld < 0 or thrld > 1:
            raise ValueError("thrld must range between 0 and 1")

        thrld = 1-thrld
        
        # filter the matrices for non variance
        a,b = filterFVA(min_mat,max_mat)

        # calc correlation matrix
        min_max = pd.concat([a,b],1)
        corMat = min_max.transpose().corr()
        corMatAb = np.absolute(corMat)

        # get linkage information
        linkage = spc.linkage(min_max, method = "complete", metric="correlation")
        idx = spc.fcluster(linkage, 0.25, "distance")
        
        # get the center of each cluster by means of the most 
        centers = list()
        for i in range(min(idx), max(idx)+1):
            sel = idx==i
            mat = corMatAb.loc[sel,sel]
            sums = np.sum(mat,1)
            n_i = np.where(sums == np.amax(sums))
            name = mat.index[n_i[0][0]]
            centers.append(name)
        
        # return the matrices for the cluster centers and the cluster association
        cluster = {corMatAb.index[i]:x for i,x in enumerate(idx)}
        return(min_mat.loc[centers,:], max_mat.loc[centers,:], cluster)



    def calcSampleTaub(self, min1, min2, max1, max2):
        # get the distance for each reaction in two FVAs

        idx = range(np.size(min1))
        d = np.max(pd.DataFrame({"x1" : min1-max1,
            "x2" : min1 - max2,
            "x3" : min2 - max1,
            "x4" : min2 - max2},
            index = idx), 1)
        m = np.mean(pd.DataFrame({"x1" : max1-min1,
            "x2" : max2 - min2},
            index = idx),1)

        nonzero = list(m != 0)
        d.loc[nonzero] = d.loc[nonzero]/m.loc[nonzero]
        d.loc[[not x for x in nonzero]] = -1
        d = d+2
        d = np.log2(d)
        return(d)

    def calcSampleMoors(self, min1, min2, max1, max2):
        # Karlis idea to implement the distances
        
        idx = range(np.size(min1))
        d =np.sum( pd.DataFrame({"min" : np.absolute(min1-min2),
            "max": np.absolute(max1-max2)},
            index = idx), 1)
        m = np.sum(pd.DataFrame({"x1" : max1-min1,
            "x2" : max2 - min2},
            index = idx),1)
        nonzero = list(m != 0)
        d.loc[nonzero] = d.loc[nonzero]/m.loc[nonzero]
        return(d)

    def calcSampleJacc(self, min1, min2, max1, max2):
        # calculate jaccard distance

        idx = range(np.size(min1))
        d = np.max(pd.DataFrame({"x1" : min1-max1,
            "x2" : min1 - max2,
            "x3" : min2 - max1,
            "x4" : min2 - max2},
            index = idx), 1)
        d[d>0] = -0
        d = d*-1
        m = np.sum(pd.DataFrame({"x1" : max1-min1,
            "x2" : max2 - min2},
            index = idx),1)
        nonzero = list(m != 0)
        d.loc[nonzero] = d.loc[nonzero]/m.loc[nonzero]
        return(d)
        
            

    def calcFVAdistPerSamplePair(self, min_mat, max_mat, dist_method = "Moors", filt_method = "any", cluster = True):
        # calculate a distance matrix for all sample pairs and reaction, but do the calculation for each sample across all reactions first - this is should be much faster computationally 
        # dist_method is either "Moors" or "Taub" or "Jacc"
        # filt_method is either "any" or "all"
        
        # pre filter or cluster the data
        if cluster:
            min_mat,max_mat,cluster = clusterFVA(min_mat, max_mat)
        else:
            min_mat,max_mat = filterFVA(min_mat,max_mat, method = filt_method)

        samples = min_mat.shape[1]
        rxns = min_mat.shape[0]
        d3 = np.zeros((samples,samples,rxns))
        print(d3.shape)
        
        for i in range(samples-1):
            for j in range(i+1,samples):
                min1 = np.array(min_mat)[:,i]
                max1 = np.array(max_mat)[:,i]
                min2 = np.array(min_mat)[:,j]
                max2 = np.array(max_mat)[:,j]
                if dist_method == "Moors":
                    d = calcSampleMoors(min1,min2,max1,max2)
                elif dist_method == "Taub":
                    d = calcSampleTaub(min1,min2,max1,max2)
                elif dist_method == "Jacc":
                    d = calcSampleJacc(min1,min2,max1,max2)
                else:
                    raise ValueError("method must be one of: 'Moors', 'Taub', 'Jacc'")               
                d3[i,j,:] = d3[j,i,:] = d
        return(d3, cluster)

    def calcFVAdistPerRxn(self, min_mat, max_mat, dist_method = "Moors", filt_method = "any", cluster = True):
        # calculate a distance matrix for all sample pairs and reaction, but do the calculation for each reaction across all sample pairs first - this is should be slower in computation, but has the advantage to filter results already early
        # dist_method is either "Moors" or "Taub" or "Jacc"
        # filt_method is either "any" or "all"

        # pre filter or cluster the data
        if cluster:
            min_mat,max_mat,cluster = clusterFVA(min_mat, max_mat)
        else:
            min_mat,max_mat = filterFVA(min_mat,max_mat, method = filt_method)

        samples = min_mat.shape[1]
        rxns = min_mat.shape[0]
        d3 = np.zeros((samples,samples,rxns))
        print(d3.shape)
        
        for h in range(rxns):
            for i in range(samples-1):
                for j in range(i+1,samples):
                    min1 = np.array(min_mat.iloc[h,i])
                    max1 = np.array(max_mat.iloc[h,i])
                    min2 = np.array(min_mat.iloc[h,j])
                    max2 = np.array(max_mat.iloc[h,j])
                    if dist_method == "Moors":
                        d = calcSampleMoors(min1,min2,max1,max2)
                    elif dist_method == "Taub":
                        d = calcSampleTaub(min1,min2,max1,max2)
                    elif dist_method == "Jacc":
                        d = calcSampleJacc(min1,min2,max1,max2)
                    else:
                        raise ValueError("method must be one of: 'Moors', 'Taub', 'Jacc'")               
                    d3[i,j,:] = d3[j,i,h] = d
        return(d3, cluster)


## testing stuff
#min_mat = pd.read_csv("minFVA_MatjesAbsorption.colormore22-zeroBiomass.csv", index_col=0)
#max_mat = pd.read_csv("maxFVA_MatjesAbsorption.colormore22-zeroBiomass.csv", index_col=0)
#
#min_ar = np.array(min_mat)
#max_ar = np.array(max_mat)
#
#min1 = min_mat.iloc[:,23]
#min2 = min_mat.iloc[:,1]
#max1 = max_mat.iloc[:,23]
#max2 = max_mat.iloc[:,1]
#
#%timeit calcFVAdistPerRxn(min_mat.iloc[0:20,0:20], max_mat.iloc[0:20,0:20])
#%timeit calcFVAdistPerSamplePair(min_mat.iloc[0:20,0:20], max_mat.iloc[0:20,0:20], dist_method = "Taub")
#
