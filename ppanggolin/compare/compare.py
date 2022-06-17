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
import tqdm
import operator
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
    required.add_argument('-in', "--input_file_compare", required=False, type=str, help="The list of all genomes including the two conditions")
    required.add_argument('-out', "--output_file_compare", required=False, type=str, help="The file returns the dictionnary with p-val")
    return parser

def performComparisons(pangenome, dataset1, dataset2):
    """Perform comparations based on the occurences of genes in a pangenome between two list of genomes corresponding to 2 conditions.
    Reads a pangenome object and return a dictionnary where gene families are the keys and the values are lists of 4 elements (p-value, oddsratio, V-cramer,p-value-corrected).
    
    :param pangeome: a pangenome  
    :type pangeome: :class:`ppanggolin.Pangenome`i
    :param dataset1: The name of different strains for condition1
    :type dataset1: list[str]
    :param dataset2: The name of different strains for condition2
    :type dataset2: list[str]
    :return: a dictionnary of family genes as key and list of 4 elements  p-value, oddsratio, v-cramer, p-value corrected
    :type: dict[ str , list[float, float, float, float]]
    """
      
    results={}
    uncorrected_pval = []
    for f in tqdm.tqdm(pangenome.geneFamilies, unit = "gene family"):
        dataset1_fampresence = 0
        dataset2_fampresence = 0
        dataset1_famabsence = 0
        dataset2_famabsence = 0
        org_names = set([org.name for org in f.organisms])
        for org in dataset1:
            if org in org_names:
                dataset1_fampresence+=1
            else:
                dataset1_famabsence+=1
        for org in dataset2:
            if org in org_names:
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
            pass
        results[f.name] = [0, f.namedPartition, pvalue, oddsratio, V, pval_corr]  
        uncorrected_pval.append(pvalue)
    #comprehensive lists
    # all_corrected_pvals = multipletests([results[f.name][0] for f in pangenome.geneFamilies], alpha=0.05, method='hs', is_sorted=False, returnsorted=False)
    # r=[]
    # for f in pangenome.geneFamilies:
    #     r.append(results[f.name][0])
    # all_corrected_pvals = multipletests(r, alpha=0.05, method='hs', is_sorted=False, returnsorted=False)
    all_corrected_pvals = multipletests(uncorrected_pval, alpha=0.05, method='hs', is_sorted=False, returnsorted=False)[1]
    for index, f in enumerate(pangenome.geneFamilies):
        results[f.name][5] = all_corrected_pvals[index]
        results[f.name][0] = index

    logging.getLogger().debug("end of performCompararison step")     
    return(results)
    
    
def launch(args):
    """ launch the comparison"""

    pangenome = Pangenome()
    logging.getLogger().debug("start step")     
    if os.stat(args.input_file_compare).st_size == 0:
        raise Exception(f"Your provided an empty file")
    pangenome.addFile(args.pangenome)
    checkPangenomeInfo(pangenome, needFamilies=True, needAnnotations=True, disable_bar=args.disable_prog_bar)
    (data1, data2, c1, c2) = extract_condition(args.input_file_compare)
    intersec_name = set(data1).intersection(set(data2))
    if len(intersec_name) != 0:
        raise Exception(f"You provided same genomes, must be different argument to compare common part : '{intersec_name}'")
    if c1== c2:
        raise Exception (f"You provided same conditions in file, must be different arguments common part : '{args.condition1}'")
    rs = performComparisons(pangenome = pangenome, dataset1 = data1, dataset2 = data2)
    output_file(args.output_file_compare, rs)
    #writePangenome(pangenome, pangenome.file, args.force, disable_bar=args.disable_prog_bar)
    logging.getLogger().debug("end of launch step")     

def extract_condition(file):
    """ extract condition names and list of genomes for these two condition from matrice file 
         :param file: contains genome_names, condition1 and condition2 with 1 if the genome is present for condition1 or condition2 and 0 if the genome is absent
         :type os.File: tsv file 
    """

    list_genomes_presence_cond1 = []
    list_genomes_presence_cond2 = []
    with open(file,"r") as tsvfile :
        for i, line in enumerate(tsvfile) :
            elements = [e.strip() for e in line.split('\t')]
            genomes_name = elements[0]
            if i == 0:
                condition1 = elements[1]
                condition2 = elements[2]  
            else :
                if (elements[1] == "1"):
                    list_genomes_presence_cond1.append(genomes_name)
                if (elements[1] == "0"):
                    list_genomes_presence_cond2.append(genomes_name)
                if (elements[1] == "1" and elements[2] == "1") :
                    raise Exception(f"Genome cannot be present at same time for two conditions : '{genomes_name}") 
        logging.getLogger().debug("end of extraction step")     
        return(list_genomes_presence_cond1, list_genomes_presence_cond2, condition1, condition2)

def output_file(file2, res):
     """ extract condition names and list of genomes for these two condition from matrice file
          :param file: contains genome_names, condition1 and condition2 with 1 if the genome is present for condition1 or condition2 and 0 if the genome is absent
          :type os.File: tsv file
     """

     with open(file2+".tsv", "w") as tsvfile:
         tsvfile.write("Id\tGene_family_name\tpartition\tp_value\toddsratio\tV\tpval_cor\n")
         for fam_name, stat_values in sorted(res.items(), key = lambda x : x[1][2]):
             index = str(stat_values[0])              
             partition = str(stat_values[1])              
             p_value_str = str(stat_values[2])
             oddsratio_str = str(stat_values[3])
             V_str = str(stat_values[4])
             p_value_corrected_str = str(stat_values[5])
             tsvfile.write(index+"\t"+fam_name+"\t"+partition+"\t"+p_value_str+"\t"+oddsratio_str+"\t"+V_str+"\t"+p_value_corrected_str+"\n")

