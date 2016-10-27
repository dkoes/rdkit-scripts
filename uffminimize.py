#!/usr/bin/python

import sys,string
from rdkit.Chem import AllChem as Chem
from optparse import OptionParser

#minimize sdf structures with the UFF forcefiled and rdkit

parser = OptionParser(usage="Usage: %prog [options] <input>.sdf <output>.sdf")
parser.add_option("-v","--verbose", dest="verbose",action="store_true",default=False,
                  help="verbose output")


    
(options, args) = parser.parse_args()
input = args[0]
output = args[1]

inmols = Chem.SDMolSupplier(input)
if inmols is None:
    print "Could not open ".input
    sys.exit(-1)
    
sdwriter = Chem.SDWriter(output)
if sdwriter is None:
    print "Could not open ".output
    sys.exit(-1)
        
for mol in inmols:
    if mol is not None:
        try:
            Chem.UFFOptimizeMolecule(mol)
            if options.verbose:
                e = Chem.UFFGetMoleculeForceField(mol).CalcEnergy()
                print mol.GetProp('_Name'),e
            sdwriter.write(mol)
        except (KeyboardInterrupt, SystemExit):
            raise            
        except:
            print "Exception occurred",mol.GetProp('_Name')
    else:
        print "ERROR"