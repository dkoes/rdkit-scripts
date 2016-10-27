#!/usr/bin/python

#given a smarts rxn file, a core scaffold smarts file (with name containing connecting atoms)
# and an sdf file, extract the matching scaffolds from of the sdf file and cluster
#them greedily to identify a set of scaffold conformations

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

#compute the distances between the matching connecting atoms 
#and return true if all distance are small enough
def checkConnect(center, cmatch, mol, match, connectIndices, connect):
	cconf = center.GetConformer(0)
	mconf = mol.GetConformer(0)
	for i in connectIndices:
		cidx = cmatch[i]
		midx = match[i]
		cpt = cconf.GetAtomPosition(cidx)
		mpt = mconf.GetAtomPosition(midx)
		dist = cpt.Distance(mpt)
		if dist > connect:
			return False
	return True

#find all the scaffolds in mols that are within rmsd of center and where the connecting
#atoms are within connect, return new cluster and new mols
def createCluster(center,cmatch, mols, pattern, core, rmsd, connectIndices, connect):
	cluster = list()
	newmols = list()
	for (mol,match) in mols:
		r = AllChem.GetBestRMS(mol,center,maps=[zip(cmatch,match)])
		if r < rmsd and checkConnect(center, cmatch, mol,match,connectIndices, connect):
			cluster.append((mol,match,r))
		else:
			newmols.append((mol,match))
	cluster.sort(key = lambda (m,mtch,r): r )
	return (cluster, newmols)

#find the mol in mols that has the maximum minimum distance between the first
#mol in each cluster
#ACTUALLY, for the tight tolerances we need, this really doesn't make a difference
#and just slows things down, so just pick the first available conformer
def computeNext(clusters,mols):
	if len(mols) > 0:
		return mols[0]
	else:
		return (None,None)
	max = 0
	best = (None,None)
	for (mol,match) in mols:
		min = float('inf')
		for cl in clusters:
			cmol = cl[0][0]
			cmatch = cl[0][1]
			r = AllChem.GetBestRMS(cmol,mol,maps=[zip(match,cmatch)])
			if r < min:
				min = r
		if min > max:
			max = min
			best = (mol,match)
	return best

#MAIN
if len(sys.argv) < 5:
	print "Need reaction file, core scaffold file, sdf input file and sdf output"
	sys.exit(1)
	
parser = argparse.ArgumentParser()
parser.add_argument('-r','--rxn', help="Reaction file")
parser.add_argument('-c','--core',help="Core scaffold with connecting atoms in name")
parser.add_argument('-i','--input',help="Input conformers")
parser.add_argument('-o','--output',help="Clustered core scaffold output")
parser.add_argument("--rmsd",type=float,default=0.5,help="Maximum RMSD for cluster membership")
parser.add_argument("--connect",type=float,default=0.1,help="Maximum allowed deviation of connecting atoms for cluster membership")
parser.add_argument("--sample",type=int,default=1,help="Amount to sample conformations")
args = parser.parse_args()

rxnf = open(args.rxn)
rxnsm = rxnf.readline().split()[0] #ignore any name
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

coref = open(args.core)
corel = coref.readline()
coreconnects = corel.split()[1:]
core = AllChem.MolFromSmarts(corel.split()[0])

inmols = AllChem.SDMolSupplier(args.input)
if inmols is None:
	print "Could not open ",args.input
	sys.exit(-1)

sdwriter = AllChem.SDWriter(args.output)
if sdwriter is None:
	print "Could not open ",args.output
	sys.exit(-1)

smart = AllChem.MolToSmarts(product)
pattern = AllChem.MolFromSmarts(smart)

#figure out the indices of connected atoms in the smart core pattern
connectIndices = list()
for c in coreconnects:
	cm = AllChem.MolFromSmarts(c)
	a = cm.GetAtoms()[0]
	if a.HasProp('molAtomMapNumber'):
		mapnum = a.GetProp('molAtomMapNumber')
		for sma in core.GetAtoms(): 
			if sma.HasProp('molAtomMapNumber') and sma.GetProp('molAtomMapNumber') == mapnum:
				connectIndices.append(sma.GetIdx())

#read all core scaffold molecules into memory
mols = list()
cnt = 0
for mol in inmols:
	if cnt % args.sample == 0 and mol is not None:
		try:
			mol = AllChem.AddHs(mol)
			match = mol.GetSubstructMatch(pattern) #just one? why not, we're only sampling
			if match:
				sub = subMol(mol, match)
				cmatch = sub.GetSubstructMatch(core)
				if cmatch:
					sub = subMol(sub,cmatch)
					mols.append((sub,sub.GetSubstructMatch(core)))
		except (KeyboardInterrupt, SystemExit):
			raise		
		except Exception as e:
			print "Exception occurred",mol.GetProp('_Name'),e
	cnt += 1

if len(mols) == 0:
	print "No molecules!"
	sys.exit(-1)
print "Done reading"	
clusters = list() #these are just defined by a list of all the scffolds assigned to the cluster
(center, cmatch) = mols[0]

while len(mols) > 0:
	(cluster, mols) = createCluster(center,cmatch,mols, pattern, core, args.rmsd, connectIndices, args.connect)
	clusters.append(cluster)
	(center, cmatch) = computeNext(clusters,mols)

print len(clusters)
for cl in clusters:
	cmol = cl[0][0]
	cmol.SetProp("ClusterSize",str(len(cl)))
	AllChem.GetBestRMS(clusters[0][0][0],cmol) #align to very first 
	sdwriter.write(cmol)
sdwriter.close()

