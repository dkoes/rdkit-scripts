#!/usr/bin/env python

'''A script for enumerate all possible stereo-isomers around chiral centers.
Stereo bonds are not enumerate currently since they are more problematic
to deal with in RDKit.'''

import rdkit
from rdkit.Chem import AllChem as Chem
import itertools,sys

def make_iso_smiles(mol):
    '''enumerate chiral centers and return list of smiles'''
    #for some reason rdkit can't detect chiral centers without a structure
    Chem.EmbedMultipleConfs(mol)
    Chem.rdmolops.AssignAtomChiralTagsFromStructure(mol)
    chirals = [] # chiral atoms
    Chem.rdmolops.AssignStereochemistry(mol) #make sure this is done
    for a in mol.GetAtoms():
        if a.GetChiralTag() == Chem.rdchem.CHI_TETRAHEDRAL_CCW or a.GetChiralTag() == Chem.rdchem.CHI_TETRAHEDRAL_CW:
            chirals.append(a)
    cmask = itertools.product(*[[0,1]]*len(chirals))
    if len(chirals) == 0:
        return [Chem.MolToSmiles(mol)]
    else:
        ret = []
        for m in cmask:
            for i in xrange(len(chirals)):
                if m[i]:
                    chirals[i].SetChiralTag(Chem.rdchem.CHI_TETRAHEDRAL_CCW)
                else:
                    chirals[i].SetChiralTag(Chem.rdchem.CHI_TETRAHEDRAL_CW)
            ret.append(Chem.MolToSmiles(mol,True))
        return ret


if len(sys.argv) < 2:
    print "Need smiles file(s)"
    sys.exit(1)
    
for f in sys.argv[1:]:
    for line in open(f):
        vals = line.split()
        smi = vals[0]
        name = ''
        if len(vals) > 1:
            name = ' '.join(vals[1:])
        mol = Chem.MolFromSmiles(smi)
        smis = make_iso_smiles(mol)
        for smi in smis:
            print smi,name
