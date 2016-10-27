#!/usr/bin/env python

import argparse, sys, rdkit, collections
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolAlign
from rdkit.Chem.rdShape import Align
from rdkit.Chem.rdMolAlign import AlignMol
from rdkit.Chem import rdMolAlign, rdShapeHelpers

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Molecular shape comparison')
    parser.add_argument('ref',help='Reference sdf file for comparison')
    parser.add_argument('test',help="SDF file to compare with")
    parser.add_argument('-c','--collate',action='store_true',help='Only output a single value per unique mol name')
    parser.add_argument('-t','--tanimoto',action='store_true',help='Output Tanimoto distance')
    args = parser.parse_args()

    #load reference file molecules
    refmols = [mol for mol in Chem.SDMolSupplier(args.ref)]

    collated = collections.defaultdict(list)
    #for each test mol compare to all refmols
    for mol in Chem.SDMolSupplier(args.test):
        try:
            vals = []
            for r in refmols:
                o3a = rdMolAlign.GetO3A(r, mol)
                o3a.Align()                
                if args.tanimoto:
                    score = 1.0 - rdShapeHelpers.ShapeTanimotoDist(r, mol)
                else:
                    score = o3a.Score()
                vals.append(score)
            tc = max(vals)
            if args.collate:
                collated[mol.GetProp("_Name")].append(tc)
            else:
                print mol.GetProp("_Name"), tc
        except:
            pass
        
    if args.collate:
        namevals = [ (name, max(vals)) for (name, vals) in collated.iteritems()]
        namevals.sort(key=lambda (n,v): v)
        
        for (n,v) in namevals:
            print n,v
