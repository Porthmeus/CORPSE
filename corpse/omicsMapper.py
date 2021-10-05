# Porthmeus
# 25.08.21

import pandas as pd
import cobra as cb
import numpy as np
import warnings
import multiprocessing
import joblib

class omicsMapper:
    def __init__(self):
        self.name = "omicsMapper"

    def parseData(self, model, dataframe, column = 0, protein = False):
        '''Takes a pandas.DataFrame with gene/protein expression data and a cobra.Model and extracts only those data in the DataFrame where the gene/protein ids of the model matches the index of the data frame and returns the values in a dic with index:value pairs.
    Keyword arguments:
        @ model - a cobra.Model object representing the metabolic model
        @ dataframe - a pandas.DataFrame containing the expression data, the index should correspond to the gene/protein names in the metabolic model
        @ column - index of the column containing the expression data - either an integer corresponding to the column or the column name
        @ protein - whether to map the expression to the gene.ids of the model or to the gene-product - basically a switch between model.gene.id (False, default) and model.gene.name (True).

    Value:
        A dictionary with gene:value pairs.
        '''
       
        # get all genes in the model
        if protein:
            genes = np.unique([x.name for x in model.genes if x.name in dataframe.index])
            missing = np.unique([x.name for x in model.genes if x.name not in dataframe.index])
        else:
            genes = np.unique([x.id for x in model.genes if x.id in dataframe.index])
            missing = np.unique([x.id for x in model.genes if x.id not in dataframe.index])
        
        
        # get values from the data frame
        if type(column) == str:
            df = dataframe.loc[genes,column]
        elif type(column) == int:
            df = dataframe.loc[genes].iloc[:,column]
        else:
            warnings.warn("column argument was neither integer nor string - try to convert to integer")
            try:
                column = int(column)
                df = dataframe.loc[genes].iloc[:,column]
            except:
                warnings.warn("could not convert columns to integer and extract column from the data frame - will fall back to column = 0")
                df = dataframe.loc[genes].iloc[:,0]
        
        # convert the results into dictionary and add all the missing values with exrpression = 0
        df = dict(df)
        for x in missing:
            df[x] = 0
        
        return(df)


    def mapGPR(self, gpr, expression, orIsSum = True):
        '''The function evaluates a GPR exrpession and maps the expression to it.
    Keyword arguments:
        @ gpr - the gene reaction rule to evaluate
        @ expression - a dictionary containing gene:expression value pairs
        @ osIsSum - how should the OR operator of the GPR be handled - if True, the pairs of an OR operator are summed, if False, the maximum value of the two is used
    Value:
        Returns a single two values - the first representing the evaluation of the rule with the given expression data, the second value is the length of the GPR in characters (its needed for recursive calls in the function)
    Note: The GPR is evaluated from left to right, meaning if no parenthesis is set, the operators are evaluated in the order of appearance. This means that: "gene1 and gene2 or gene3" is the same as "(gene1 and gene2) or gene3", but is different from "gene1 and (gene2 or gene3)".
        '''
        
        # sanity check - if string is empty return 0
        if gpr == "":
            return(0,len(gpr))

        n = 0
        pairs = [None,None] # contains the values of which should be compared by an operator
        word = "" # the gene name in GPR
        op = "" # contains the operation on the pairs
        # go through the complete GPR string
        while n < len(gpr):
            # debug
            #print("BEGIN: "+word + " - " +str(pairs) + " - " +op + " - " +str(n)+ " - " + gpr[n])

            # the idea is following: go through the string if one encounters either
            # a space, "(", or ")" - for each of the cases evaluate whether two
            # matching pairs + operator has been found and calculate the result for
            # each pair
            # evaluate the space as a separator
            if gpr[n] == " ":
                # make sure its no just an additional space before a parenthesis or
                # in the beginning of the string
                # if the space was only an additional space - remove it
                if word in [""," ","(",")"]:
                    word = ""
                else:
                    # check whether the last word was a gene name or an operator
                    # if the word was an operator check which operator it was and
                    # what kind of evaluating function has to be used
                    if word.lower() == "and":
                        op = "min"
                        word = ""
                    elif word.lower() == "or":
                        if orIsSum:
                            op = "sum"
                            word = ""
                        else:
                            op = "max"
                            word = ""

                    else: 
                        # if it is a gene name - check which one of the pairs it is
                        if pairs[0] == None:
                            pairs[0] = expression[word]
                            word = ""
                        # if the word was the second of an evaluation pair,
                        # evaluate the pair according to the operator and set it to
                        # the first value of the next pair
                        elif pairs[1] == None:
                            pairs[1] = expression[word]
                            word = ""
                            pairs[0] = self.compareBinary(pairs, function = op)
                            pairs[1] = None
                            op = ""
                        else:
                            warnings.warn("Both pairs are not None and the word {word} was found".format(word = word))
            # if one encounters a closed bracket - return the current value of the
            # first pair - which should be the stored result of all evaluations of
            # operators in these parenthesis
            elif gpr[n] == ")":
                # if there is no word saved, the expression was evaluated already
                if word == "" and op == "":
                    val = pairs[0]
                # otherwise evaluate it and return the result
                else:
                    pairs[1] = expression[word]
                    val = self.compareBinary(pairs, op)
                #print("returned:" + str(val) +", "+str(n))
                return(val,n)
            
            # if one encounters an opening bracket - recursively evaluate the
            # internal of the bracket and return the result
            elif gpr[n] == "(":
                val,nn = self.mapGPR(gpr = gpr[n+1:], expression = expression, orIsSum= orIsSum)
                n = n+nn+1 # adjust the counter to after the bracket
                # write the result in the correct position of the pair
                if pairs[0] == None:
                    pairs[0] = val
                    word = ""
                elif pairs[1] == None:
                    pairs[1] = val
                    word = ""
                    pairs[0] = self.compareBinary(pairs, function = op)
                    pairs[1] = None
                    op =""
                else:
                    raise ValueError("found '(' in the GPR but pairs list is already populated")
            else:
                word = word + gpr[n]

            # debug
            #print("END: "+word + " - " +str(pairs) + " - " +op + " - " +str(n)+ " - " + gpr[n])

            # set iterator
            n = n+1 
        
        # if only one word is given you end up here with no value in pairs
        #print(str(pairs) + " - " + word + " - " +op)
        if pairs[0] == None and pairs[1] == None:
            val = expression[word]
        # add the final word if not done yet and evaluat the final pair
        elif pairs[0] != None and pairs[1] == None and word != "":
            pairs[1] = expression[word]
            val = self.compareBinary(pairs, function = op)
        else:
            val = pairs[0]

        # debug
        #print("returned:" + str(val) +", "+str(n))
        return(val,n)


    def compareBinary(self, pairs, function):
        ''' compares two values either by sum, min or max - helper function for evalGPR '''
        if function == "sum":
            return(sum(pairs))
        elif function == "min":
            return(min(pairs))
        elif function == "max":
            return(max(pairs))
        elif function == "":
            print(pairs)
            raise ValueError("dont know")
            return(pairs[0])
        else:
            raise ValueError("Argument '{f}' for function is not defined".format(f = function))

    def mapRxn(self, rxn, expression, protein =False, orIsSum = True):
        '''same as mapGPR only takes reaction object from a cobra.Model to map the expression values'''
        # check whether protein or gene reaction rule should be evaluated
        if protein:
            gpr = rxn.gene_name_reaction_rule
        else:
            gpr = rxn.gene_reaction_rule
        val,n = self.mapGPR(gpr = gpr, expression = expression, orIsSum= orIsSum)
        # return the value
        return(val)

    def mapSampleToModel(self, model, dataframe, column = 0, protein = False, orIsSum = True):

        dfx = self.parseData(model, dataframe = dataframe, column = column, protein = protein)
        vals = [self.mapRxn(rxn, expression = dfx, protein= protein, orIsSum= orIsSum) for rxn in model.reactions]
        return(vals)

    def mapExpressionToReaction(self,
            model,
            dataframe,
            column = None,
            protein = False,
            orIsSum = True,
            num_cores = multiprocessing.cpu_count()-1):
        '''Maps expression values to a reactions of the a cobra.Model object. It takes the expression values as pandas.DataFrame and will map all columns to the model reactions
    Keyword arguments:
        @ model - a cobra.Model object representing the metabolic model
        @ dataframe - a pandas.DataFrame containing the expression data, the index should correspond to the gene/protein names in the metabolic model
        @ column - a list of indeces for the columns of the dataframe containing the expression data - either as integer corresponding to the column or string of the column name
        @ protein - whether to map the expression to the gene.ids of the model or to the gene-product - basically a switch between model.gene.id (False, default) and model.gene.name (True).
        @ osIsSum - how should the OR operator of the GPR be handled - if True, the pairs of an OR operator are summed, if False, the maximum value of the two is used
        @ num_cores - how many cores should be used - defaults to all available cores -1
    Value:
        A pandas.DataFrame containing reactions in rows and samples as columns containing the rxn activity as expression value.
    Note: The GPR is evaluated from left to right, meaning if no parenthesis is set, the operators are evaluated in the order of appearance. This means that: "gene1 and gene2 or gene3" == "(gene1 and gene2) or gene3" != "gene1 and (gene2 or gene3)".
        '''
        

        # get all reactions of the model
        rxns = [x.id for x in model.reactions]

        # sanity check for the columns
        if column == None:
            column = dataframe.columns
        # if column is not a list yet, try to cast it to one
        if type(column) != list:
            column = list(column)
        if any(dataframe.shape[1] < x for x in column if type(x) == int):
            raise ValueError("Subscript out of bounds - found index in columns, which is larger than the dimension of the expression dataframe")
        if any([x not in dataframe.columns for x in column if type(x) == str]):
            raise ValueError("Subscript out of bounds - found index in columns, which does not match any of the columns in the expression dataframe")

        # remove some memory footprint by reducing the expression data to only the relevant genes
        if protein:
            genes = [gene.name for gene in model.genes]
        else:
            genes = [gene.id for gene in model.genes]

        dataframe= dataframe.loc[[x for x in dataframe.index if x in genes]]

        # get the sample names and create an empty dataframe
        colnames = [x if type(x) == str else dataframe.columns[x] for x in column]
        rxn_exp = pd.DataFrame(np.zeros((len(rxns),len(colnames))), index = rxns, columns = colnames,dtype = float)

        # iterate through the samples and get the relevant expression values - prefer "threads" is speeding up the process enormously
        results = joblib.Parallel(n_jobs = num_cores, verbose = 10, prefer ="threads")(joblib.delayed(self.mapSampleToModel)(model = model,
            dataframe = dataframe,
            column = sample,
            protein = protein,
            orIsSum = orIsSum) for sample in colnames)
        
        
        results = pd.DataFrame(np.array(results).transpose(),
                index = [rxn.id for rxn in model.reactions],
                columns = colnames)
        return(results)
