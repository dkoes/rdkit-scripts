#!/usr/bin/python

#given a smarts rxn file (backwards), a core scaffold smarts file (with name containing connecting atoms)
# an sdf file with the results of clustering scaffolds and an input sdf file and an
#output prefix, output the extracted reactant conformations aligned to their scaffold

import sys,gzip,argparse
from rdkit.Chem import AllChem

def subMol(mol, match):
	#not sure why this functionality isn't implemented natively
	#but get the interconnected bonds for the match
	atoms = set(match)
	bonds = set()
	for a in atoms:
		atom = mol.GetAtomWithIdx(a)
		for b in atom.GetBonds():
			if b.GetOtherAtomIdx(a) in atoms:
				bonds.add(b.GetIdx())
	return AllChem.PathToSubmol(mol,list(bonds))

#return index of closest scaffold in scaffold to mol
def closestScaffold(scaffolds, pattern, core, mol):
	ret = -1
	match = mol.GetSubstructMatch(pattern) 
	if match:
		sub = subMol(mol, match)
		cmatch = sub.GetSubstructMatch(core)
		if cmatch:
			min = float('inf')			
			for (i,(s,smatch)) in enumerate(scaffolds):
				r = AllChem.GetBestRMS(s, sub, maps=[zip(cmatch,smatch)])
				if r < min:
					min = r
					ret = i
					mmatch =  mol.GetSubstructMatch(core)
					AllChem.GetBestRMS(s,mol,maps=[zip(mmatch,smatch)])
	return ret

#MAIN

parser = argparse.ArgumentParser()
parser.add_argument('-r','--rxn', help="Reaction file")
parser.add_argument('-c','--core',help="Core scaffold with connecting atoms in name")
parser.add_argument('-i','--input',help="Input conformers")
parser.add_argument('-s','--scaffolds',help="Scaffold conformers")
parser.add_argument('-o','--output',help="Output prefix")

args = parser.parse_args()

rxnf = open(args.rxn)
rxnsm = rxnf.readline().split()[0] #ignore any name
rxn = AllChem.ReactionFromSmarts(rxnsm)
rxn.Initialize()

if rxn.GetNumReactantTemplates() != 1:
	print "Need backwards reaction"
	sys.exit(-1)
	

coref = open(args.core)
corel = coref.readline()
coreconnects = corel.split()[1:]
core = AllChem.MolFromSmarts(corel.split()[0])

inscaffolds = AllChem.SDMolSupplier(args.scaffolds,False)
if inscaffolds is None:
	print "Could not open ",args.scaffolds
	sys.exit(-1)

inmols = AllChem.SDMolSupplier(args.input)
if inmols is None:
	print "Could not open ",args.input
	sys.exit(-1)

smart = AllChem.MolToSmarts(rxn.GetReactantTemplate(0))
pattern = AllChem.MolFromSmarts(smart)

#read in scaffolds
scaffolds = list()
for mol in inscaffolds:
	#compute match of core
	 cmatch = mol.GetSubstructMatch(core)
	 scaffolds.append((mol,cmatch))
	 
#setup output file, one for each reactant product
outputs = list()
for i in xrange(rxn.GetNumProductTemplates()):
	outputs.append(list())
	for j in xrange(len(scaffolds)):
		sdwriter = AllChem.SDWriter("%s_%d_%d.sdf" % (args.output,i,j))
		outputs[i].append(sdwriter)

for mol in inmols:
	#for each mol, decompose it into its reactants
	mol = AllChem.AddHs(mol)
	#figure out which scaffold conformation is closest
	c = closestScaffold(scaffolds, pattern, core, mol)
	prods = rxn.RunReactants([mol])
	if c >= 0:
		for p in prods: #there may be multiple possible products
			for (i,react) in enumerate(p):
				react = AllChem.RemoveHs(react)
				react.SetProp('_Name',AllChem.MolToSmiles(react))
				outputs[i][c].write(react)
		
		 