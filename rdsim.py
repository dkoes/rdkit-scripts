#!/usr/bin/env python

import sys,string,argparse
from rdkit.Chem import AllChem as Chem
from optparse import OptionParser
import os, gzip
from rdkit import DataStructs

parser = argparse.ArgumentParser(description='Compute similarity to query ligand')
parser.add_argument('query',metavar='query ligand')
parser.add_argument('test',metavar='test ligands')

args = parser.parse_args()

query = open(args.query).readline().split()[0].replace('$','')
qmol = Chem.MolFromSmiles(query)
qfp = Chem.RDKFingerprint(qmol)

for line in open(args.test):
    try:
        smi = line.split()[0].replace('$','')
        mol = Chem.MolFromSmiles(smi)
        fp = Chem.RDKFingerprint(mol)
        print smi,DataStructs.FingerprintSimilarity(qfp,fp)
    except:
        pass #ignore rdkit errors
