#!/usr/bin/env python3
# coding:utf-8

# default libraries
import sys
if sys.version_info < (3, 6):  # minimum is python3.6
    raise AssertionError("Minimum python version to run PPanGGOLiN is 3.6. Your current python version is " +
                         ".".join(map(str, sys.version_info)))
import argparse
import logging
import pkg_resources
import tempfile
import os
import pdb

# installed libraries
import numpy as np
from scipy.stats import ttest_rel, ttest_ind, wilcoxon, ranksums, iqr, fisher_exact, chi2_contingency
from statsmodels.stats.multitest import multipletests

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats import checkPangenomeInfo, writePangenome, ErasePangenome

def compareSubparser(subparser): 
    parser = subparser.add_parser("compare", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group(title="Required arguments", description="One of the following arguments is required :")
    required.add_argument('-p', '--pangenome', required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-i', "--input", required=False, type=str, help="The list of all genomes including the two conditions")
    required.add_argument('-o', "--output", required=False, type=str, help="The file returns the dictionnary with p-val")

    return parser

def performComparisons(pangenome, dataset1, dataset2):
    """Perform comparations based on the occurences of genes in a pangenome between two list of genomes corresponding to 2 conditions.
    Reads a pangenome object and return a dictionnary where gene families are the keys and the values are lists of 3 elements (p-value, oddsratio, V-cramer).
    
    :param pangeome: a pangenome  
    :type pangeome: :class:`ppanggolin.Pangenome`
    :param dataset1: The name of different strains for condition1
    :type dataset1: list[str]
    :param dataset2: The name of different strains for condition2
    :type dataset2: list[str]
    :return: a dictionnary of family genes as key and list of 3 elements  p-value, oddsratio, v-cramer, p-value corrected
    :type: dict[ str , list[float, float, float, float] ]
    """

    results={}
    uncorrected_pval = []
    for f in pangenome.geneFamilies:
        dataset1_fampresence = 0
        dataset2_fampresence = 0
        dataset1_famabsence = 0
        dataset2_famabsence = 0
        sub_org_list = [org.name for org in pangenome.organisms if ((org.name in dataset1) or (org.name in dataset2))]
        for org in sub_org_list:
            if org in dataset1:
                dataset1_fampresence+=1
            else:
                dataset1_famabsence+=1
            if org in dataset2:
                dataset2_fampresence+=1
            else:
                dataset2_famabsence+=1
            #table of contingency  
            contingency_table = np.array([[dataset1_fampresence, dataset2_fampresence], [dataset1_famabsence, dataset2_famabsence]])
            
        def cramerV(ct):
            X2 = chi2_contingency(ct, correction=False)[0]
            n = np.sum(ct)
            minDim = min(ct.shape)-1
            return(np.sqrt((X2/n) / minDim))

        pvalue, oddsratio, V, pval_corr = (float("nan"), float("nan"), float("nan"), float("nan"))

        try:
            oddsratio, pvalue = fisher_exact(contingency_table)
            V = cramerV(contingency_table)
        except:
            #print(fam+" "+str(contingency_table))
            pass
        results[f.name] = [pvalue, oddsratio, V, pval_corr]  
        uncorrected_pval.append(pvalue)
    #comprehensive lists
    # all_corrected_pvals = multipletests([results[f.name][0] for f in pangenome.geneFamilies], alpha=0.05, method='hs', is_sorted=False, returnsorted=False)
    # r=[]
    # for f in pangenome.geneFamilies:
    #     r.append(results[f.name][0])
    # all_corrected_pvals = multipletests(r, alpha=0.05, method='hs', is_sorted=False, returnsorted=False)

    all_corrected_pvals = multipletests(uncorrected_pval, alpha=0.05, method='hs', is_sorted=False, returnsorted=False)[1]
    for index, f in enumerate(pangenome.geneFamilies):
        #pdb.set_trace()
        results[f.name][-1] = all_corrected_pvals[index]
        #print(index)
    return(results)
    
def launch(args):
    """ launch the comparison"""
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    checkPangenomeInfo(pangenome, needFamilies=True, needAnnotations=True, disable_bar=args.disable_prog_bar)
    #intersec_name = set(args.dataset1).intersection(args.dataset2)
    #if len(intersec_name) != 0:
        #raise Exception(f"You provided same genomes, must be different argument to compare common part : '{intersec_name}'")
    #if args.condition1== args.condition2:
        #raise Exception (f"You provided same conditions, must be different arguments common part : '{args.condition1}'")
    #if os.stat(args.file).st_size == 0:
        #raise Exception(f"Your provided an empty file")
    extract_condition(file=args.input)
    (data1, data2)= extract_condition(file)
    performComparisons(pangenome, data1, data2)
    output_file(file2= args.output, performComparisons(pangenome, dataset1=args.input, dataset2=args.input))
    #writePangenome(pangenome, pangenome.file, args.force, disable_bar=args.disable_prog_bar)

def extract_condition(file):
     """ extract condition names and list of genomes for these two condition from file 
         :param file: contains genome_names, condition1 and condition2 with 1 if the genome is present for condition1 or condition2 and 0 if the genome is absent
         :type os.File: tsv file """
    list_condition1=[]
    list_condition2=[]
    with open(file,"r") as tsvfile :
        for i, line in enumerate(tsvfile) :
            elements = line.split('\t')
            genome_name = elements[0]
            if i == 0:
                condition1 = elements[1]
                condition2 = elements[2]  
            else:
                if (elements[1] == 1):
                    list_condition1.append(genome_name)
                if (elements[2] == 1):
                    list_condition2.append(genome_name)
                #prévoir un cas on a 1 pour les 2 conditions
                                          
        return(list_conition1, list_condition2, condition1, condition2)

def ouput_file(file2, res): 
    with open (file2, "w") as tsvfile:
        for fam_name in res:
            file2.write(fam_name+"\t"+res[fam_name][0])
