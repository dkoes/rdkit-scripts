#!/usr/bin/env python

'''Given a pdb file and an sdf file of the same ligand, find the correspondance
between the atoms and output a pdb with the same atom names as the pdb input
but with the coordinates of the sdf file. Atoms missing in the sdf 
will not have their coordinates changed, but note that hydrogens are stripped'''

import sys
from rdkit.Chem import AllChem

if len(sys.argv) != 3:
    print "Need sdf and pdb files"
    sys.exit(1)
    
sdf = sys.argv[1]
pdb = sys.argv[2]

#figure out starting atom number
startnum = 1
for line in open(pdb):
    if line.startswith('ATOM') or line.startswith('HETATM'):
	startnum = int(line[6:11])
	break

pmol = AllChem.MolFromPDBFile(pdb)
smol = AllChem.SDMolSupplier("newsnappy.sdf").next()
pmol = AllChem.AssignBondOrdersFromTemplate(smol,pmol)

if smol.HasSubstructMatch(pmol):
    m = smol.GetSubstructMatch(pmol)
    pconf = pmol.GetConformer()
    sconf = smol.GetConformer()
    for (pi, si) in enumerate(m):
	pconf.SetAtomPosition(pi, sconf.GetAtomPosition(si))
    pdbout = AllChem.MolToPDBBlock(pmol, flavor=2).split('\n')
    num = startnum
    for line in pdbout:
	if line.startswith('ATOM') or line.startswith('HETATM'):
	    print '%s%5d%s' % (line[:6],num,line[11:])
	    num += 1
else:
    print "Could not matchup pdb and sdf"
    sys.exit(1)
