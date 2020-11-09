#!/usr/bin/env python3
from rdkit.Chem import AllChem as Chem
import numpy as np
import sys, gzip, re

rname = sys.argv[1]
lname = sys.argv[2]
rec = Chem.MolFromPDBFile(rname)  

if lname.endswith('.gz'):
    lig = next(Chem.ForwardSDMolSupplier(gzip.open(lname),sanitize=False))
else:    
    lig = next(Chem.SDMolSupplier(sys.argv[2],sanitize=False))

complex = Chem.CombineMols(rec,lig)
Chem.SanitizeMol(complex)

ff = Chem.UFFGetMoleculeForceField(complex,ignoreInterfragInteractions=False)
print("Before:",ff.CalcEnergy())

for p in range(rec.GetNumAtoms()):
    ff.AddFixedPoint(p)

ff.Minimize()

print("After:",ff.CalcEnergy())

cpos = complex.GetConformer().GetPositions()
conf = lig.GetConformer()
for (i,xyz) in enumerate(cpos[-lig.GetNumAtoms():]):
    conf.SetAtomPosition(i,xyz)

m = m = re.search(r'(.*)_(\d+).sdf(.gz)?',lname)
if m:
    outname = '%s_rec_uff_%s.sdf.gz'%(m.group(1),m.group(2))
else:
    outname = lname.replace('.gz','').replace('.sdf','')+'_rec_uff.sdf.gz'
outfile = gzip.open(outname,'wt')
out = Chem.SDWriter(outfile)
out.write(lig)
out.close()
outfile.close()




