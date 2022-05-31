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
    required.add_argument('-dn1', "--condition1", required=False, type=str)
    required.add_argument('-dn2', "--condition2", required=False, type=str)
    required.add_argument('-fi', "--file", required=False, type=str, help="The list of all_genomes_input_ppanggolin_antibio.list")

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
    :rtype: dict[ str , list[float, float, float, float] ]
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
    print(performComparisons(pangenome, dataset1=args.file, dataset2=args.file))
    file(file=args.file)
    #writePangenome(pangenome, pangenome.file, args.force, disable_bar=args.disable_prog_bar)




def file(file):
    list_susceptible=[]
    list_resistant=[]
    with open(file,"r") as tsvfile :
        for line in tsvfile :
            elements = line.split('\t')
            genome_name = elements[0]
            genome_type = elements[1]
            if ("Resistant" in genome_type): 
                list_resistant.append(genome_name)
            elif ("Susceptible" in genome_type): 
                list_susceptible.append(genome_name)
            else: 
                pass
        return(list_susceptible, list_resistant)
            