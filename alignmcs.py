#!/usr/bin/env python3

import argparse, sys, rdkit, collections
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolAlign
from rdkit.Chem.rdMolAlign import AlignMol
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdMolTransforms

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Align query molecules to a reference molecule using their maximum common substructure. May produce multiple alignments per an input.')
    parser.add_argument('ref',help='Reference SDF file for comparison')
    parser.add_argument('query',help="SDF file to compare with")
    parser.add_argument('out',help='Output SDF')
    parser.add_argument('-s','--strict',action='store_true',help='Use strict atom/bond matching')
    args = parser.parse_args()

    #load reference file molecules
    refmols = [mol for mol in Chem.SDMolSupplier(args.ref)]

    if len(refmols) > 1:
        print("Only using first reference molecule")
    refmol = refmols[0]
    
    out = Chem.SDWriter(args.out)
    for mol in Chem.SDMolSupplier(args.query):
        if args.strict:
            mcs = rdFMCS.FindMCS([refmol,mol])
        else:
            mcs = rdFMCS.FindMCS([refmol,mol],atomCompare=rdFMCS.AtomCompare.CompareAny,bondCompare=rdFMCS.BondCompare.CompareAny)
        
        submol = Chem.MolFromSmarts(mcs.smartsString)

        for refmatch in refmol.GetSubstructMatches(submol):
            for qmatch in mol.GetSubstructMatches(submol):
                AlignMol(mol, refmol, atomMap=list(zip(qmatch,refmatch)))
                out.write(mol)
