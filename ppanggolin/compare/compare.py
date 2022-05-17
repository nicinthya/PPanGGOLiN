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
import numpy
from scipy.stats import ttest_rel, ttest_ind, wilcoxon, ranksums, iqr, fisher_exact, chi2_contingency
from statsmodels.stats.multitest import multipletests

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats import checkPangenomeInfo, writePangenome, ErasePangenome

def compareSubparser(subparser): 
    parser = subparser.add_parser("compare", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group(title="Required arguments", description="One of the following arguments is required :")
    required.add_argument('-p', '--pangenome', required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-dn1', "--datasetname1", required=False, nargs=1, type=str)
    required.add_argument('-dn2', "--datasetname2", required=False, nargs=1, type=str)
    required.add_argument('-d1', "--dataset1", required=False, nargs='+')
    required.add_argument('-d2', "--dataset2", required=False, nargs='+')

    return parser

def performComparisons(pangenome, dataset1, dataset2):
    
    for f in pangenome.geneFamilies:
        dataset1_fampresence = 0
        dataset2_fampresence = 0
        dataset1_famabsence = 0
        dataset2_famabsence = 0
        sub_org_list = [org for org in pangenome.organisms if ((org in dataset1) or (org in dataset2))]
        for org in sub_org_list:
            if org in args.dataset1:
                dataset1_fampresence+=1
            else:
                dataset1_famabsence+=1
            if org in args.dataset2:
                dataset2_fampresence+=1
            else:
                dataset2_famabsence+=1
            #[f for f in pangenome.geneFamilies][0].organisms
            results[f.name]=np.array([[dataset1_fampresence, dataset2_fampresence], [dataset1_famabsence, dataset2_famabsence]])
    #[f.name for f in pangenome.geneFamilies]
    print(results)

#[f for f in pangenome.geneFamilies if f.name == "GCF_001398295.1_7396_3_21_genomic_CDS_0944"][0].organisms

def launch(args):
    """ launch the comparison"""
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    checkPangenomeInfo(pangenome, needFamilies=True, needAnnotations=True, disable_bar=True)
    intersec_name = set(args.dataset1).intersection(args.dataset2)
    intersec_condition = set(args.datasetname1).intersection(args.datasetname2)

    if len(intersec_name) != 0:
        raise Exception(f"You provided same genomes, must be different argument to compare common part : '{intersec_name}'")
    if args.datasetname1== args.datasetname2:
        raise Exception (f"You provided same list for conditions, must be different arguments common part : '{intersec_condition}'")
    if org in f.organisms :
        performComparisons(pangenome, dataset1=args.dataset1, dataset2=args.dataset2)
        writePangenome(pangenome, pangenome.file, args.force, disable_bar=args.disable_prog_bar)


