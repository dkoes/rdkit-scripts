#!/usr/bin/python

#given a rxn smarts (either forward or backwards, but there is assumed to be
#only one product) identify a minimal core scaffold that has every connecting atom
#from each reactant and only those scaffold atoms that are necessary to minimally
#connect these atoms
#this assume the input smarts is fully annoated with atom mappings for the heavy atoms
#output the smarts (with atom numbers) of the core scaffold

import sys
from rdkit.Chem import AllChem
import igraph

Debug = False

if len(sys.argv) < 1:
	print "Need reaction file"
	sys.exit(1)
if len(sys.argv) > 2:
	Debug = True
	
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

#record what atom mapping number belong to each reactant
reactantMapNumbers = dict()
for (i, react) in enumerate(reactants):
	for a in react.GetAtoms():
		if a.HasProp('molAtomMapNumber'):
			mapnum = a.GetProp('molAtomMapNumber')
			reactantMapNumbers[mapnum] = i

#now create an igraph from the product
molgraph = igraph.Graph()
molgraph.add_vertices(product.GetNumAtoms())

for a in product.GetAtoms():
	s = a.GetIdx()
	v = molgraph.vs[s]
	v['react'] = -1 #reactant this atom belongs to, -1 if none
	v['heavy'] = a.GetAtomicNum() != 1 #note this things [#1,#6] is H only, but that should be okay - if H can go there, no need to worry about scaffold connectivity 
	if a.HasProp('molAtomMapNumber'):
		mapnum = a.GetProp('molAtomMapNumber')
		if mapnum in reactantMapNumbers:
			v['react'] = reactantMapNumbers[mapnum]
	for na in a.GetNeighbors():
		t = na.GetIdx()
		molgraph.add_edge(s,t)

#now identify the connecting heavy atoms for reactants
connecting = list()
for v in molgraph.vs:
	if v['react'] >= 0 and v['heavy']:
		r = v['react']
		for nv in v.neighbors():
			if nv['react'] != r and nv['heavy']:
				connecting.append(v.index)
				break

if Debug:
	for i in connecting:
		print product.GetAtomWithIdx(i).GetSmarts(),molgraph.vs[i]['react']
	
#compute the shortest paths between all these connecting atoms
#the atoms along these paths make the minimal scaffold
scaffold = set()
for i in connecting:
	paths = molgraph.get_all_shortest_paths(i,connecting,mode=igraph.ALL)
	for p in paths:
		for v in p:
			scaffold.add(v)

if Debug:
	print "\nScaffold"
	for i in scaffold:
		print product.GetAtomWithIdx(i).GetSmarts(),molgraph.vs[i]['react']

#need to identify the bonds between the scaffold atoms
bonds = set()
for i in scaffold:
	a = product.GetAtomWithIdx(i)
	for b in a.GetBonds():
		end = b.GetOtherAtomIdx(i)
		if end in scaffold:
			bonds.add(b.GetIdx())

if Debug:
	print "Bonds",list(bonds)
	
subMol = AllChem.PathToSubmol(product,list(bonds))
print AllChem.MolToSmarts(subMol),
for i in connecting:
	print product.GetAtomWithIdx(i).GetSmarts(),
print ""