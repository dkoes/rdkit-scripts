#!/usr/bin/python

#given a smarts rxn file, a core scaffold smarts file and an sdf file, extract 
#the matching scaffolds from of the sdf file

import sys,gzip
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

if len(sys.argv) < 5:
	print "Need reaction file, core scaffold file, sdf input file and sdf output"
	sys.exit(1)
	
rxnf = open(sys.argv[1])
rxnsm = rxnf.readline()
rxn = AllChem.ReactionFromSmarts(rxnsm)
rxn.Initialize()

if rxn.GetNumProductTemplates() == 1:
	product = rxn.GetProductTemplate(0)
	reactants = list()
	for i in xrange(rxn.GetNumReactantTemplates()):
		reactants.append(rxn.GetReactantTemplate(i))
elif rxn.GetNumReactantTemplates() == 1:
	product = rxn.GetReactantTemplate(0)
	reactants = list()
	for i in xrange(rxn.GetNumProductTemplates()):
		reactants.append(rxn.GetProductTemplate(i))
else:
	print "Can have only one product"
	sys.exit(1)

coref = open(sys.argv[2])
core = AllChem.MolFromSmarts(coref.readline())

inmols = AllChem.SDMolSupplier(sys.argv[3])
if inmols is None:
	print "Could not open ",sys.argv[3]
	sys.exit(-1)

outf = gzip.open(sys.argv[4],'w')
sdwriter = AllChem.SDWriter(outf)
if sdwriter is None:
	print "Could not open ",sys.argv[4]
	sys.exit(-1)

smart = AllChem.MolToSmarts(core)
pattern = AllChem.MolFromSmarts(smart)
#read through input
for mol in inmols:
	if mol is not None:
		try:
			mol = AllChem.AddHs(mol)
			match = mol.GetSubstructMatch(pattern) #just one? why not, we're only sampling
			if match:
				sub = subMol(mol, match)
				cmatch = sub.GetSubstructMatch(core)
				if cmatch:
					sub = AllChem.RemoveHs(subMol(sub,cmatch))
					sdwriter.write(sub)
		except (KeyboardInterrupt, SystemExit):
			raise			
		except Exception as e:
			print "Exception occurred",mol.GetProp('_Name'),e
	else:
		print "ERROR"

sdwriter.close()
outf.close()

