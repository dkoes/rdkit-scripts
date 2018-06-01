#!/usr/bin/python

import sys,string,math
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolTransforms

import argparse, traceback
import os, gzip

'''Given a smiles file, exhaustively enumerate 3D conformers in output sdf.  
Conformers are first sampled using rdkit's distance geometry based method
combined with energy minimization.  All conformers that, with identical torsions,
have greater than an rmsd cutoff difference are kept (i.e., there are conformational
differences such as ring pucker unrelated to torsions).

For each one of these conformers, we exhaustively sample all possible torsions
at a specified degree increment (beware exponential growth!).

'''

#convert smiles to sdf
def getRMS(mol, c1,c2):
    rms = Chem.GetBestRMS(mol,mol,c1,c2)
    return rms

def getDihedralMatches(mol):
    '''return list of atom indices of dihedrals'''
    #this is rdkit's "strict" pattern
    pattern = r"*~[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]=[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-!@[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])]~*"
    qmol = Chem.MolFromSmarts(pattern)
    matches = mol.GetSubstructMatches(qmol);
    #these are all sets of 4 atoms, uniquify by middle two
    uniqmatches = []
    seen = set()
    for (a,b,c,d) in matches:
        if (b,c) not in seen:
            seen.add((b,c))
            uniqmatches.append((a,b,c,d))
    return uniqmatches

def genConformer_r(mol, conf, i, matches, degree, sdwriter):
    '''recursively enumerate all angles for matches dihedrals.  i is where is
    which dihedral we are enumerating by degree to output conformers to out'''
    print i,matches
    if i >= len(matches): #base case, torsions should be set in conf
        sdwriter.write(mol,conf)
        return 1
    else:
        incr = math.pi*degree / 180.0
        total = 0
        deg = 0
        while deg < 360.0:
            rad = math.pi*deg / 180.0
            rdMolTransforms.SetDihedralRad(mol.GetConformer(conf),*matches[i],value=rad)
            total += genConformer_r(mol, conf, i+1, matches, degree, sdwriter)
            deg += args.degree
        return total

parser = argparse.ArgumentParser(description="Usage: %prog [options] <input>.smi <output>.sdf")
parser.add_argument("--sample", help="number of conformers to sample to get non-torsional differences (default 100)", default=100, type=int, metavar="sample")
parser.add_argument("--seed", help="random seed (default 062609)", default="062609", type=int, metavar="s")
parser.add_argument("--rms_threshold", help="cutoff for considering sampled conformers the same (default 0.25)", default="0.25", type=float, metavar="R")
parser.add_argument("-v","--verbose",action="store_true",default=False, help="verbose output")
parser.add_argument("-d","--degree", type=float,help="Amount, in degrees, to enumerate torsions by (default 5.0)",default=5.0)             
parser.add_argument("--etkdg", dest="etkdg",action="store_true",default=False,
                  help="use new ETKDG knowledge-based method instead of distance geometry") 
parser.add_argument("--max_torsions",type=int,help="Skip any molecules with more than this many torsions (default 10)",default=10)
parser.add_argument("input",help="Input smi file")
parser.add_argument("output",help="Output sdf file")
   
args = parser.parse_args()

smifile = open(args.input)

if args.etkdg and not Chem.ETKDG:
    print "ETKDB does not appear to be implemented.  Please upgrade RDKit."
    sys.exit(1)
     
split = os.path.splitext(args.output)
if split[1] == '.gz':
    outf=gzip.open(args.output,'w+')
else:
    outf = open(args.output,'w+')
     
sdwriter = Chem.SDWriter(outf)    
if sdwriter is None:
    print "Could not open ".output
    sys.exit(-1)
    
for line in smifile:
    toks = line.split()
    smi = toks[0]
    name = string.join(toks[1:])    
    
    pieces = smi.split('.')
    if len(pieces) > 1:
        smi = max(pieces, key=len) #take largest component by length
        print "Taking largest component: %s\t%s" % (smi,name)
        
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        if args.verbose:
            print smi
        try:
            Chem.SanitizeMol(mol)
            mol = Chem.AddHs(mol)
            mol.SetProp("_Name",name);
            
            rotmatches = getDihedralMatches(mol)
            if(len(rotmatches) > args.max_torsions):
                print "Too many torsions (%d). Skipping %s" %(len(rotmatches),line.rstrip())
                continue
            if args.etkdg:
                cids = Chem.EmbedMultipleConfs(mol, args.sample, Chem.ETKDG(),randomSeed=args.seed)
            else:
                cids = Chem.EmbedMultipleConfs(mol, args.sample,randomSeed=args.seed)
            if args.verbose:
                print len(cids),"conformers sampled"
                
            #energy minimize all to get more realistic results
            cenergy = []            
            for conf in cids:
                #not passing confID only minimizes the first conformer
                converged = not Chem.UFFOptimizeMolecule(mol,confId=conf)
                cenergy.append(Chem.UFFGetMoleculeForceField(mol,confId=conf).CalcEnergy())
            
            #reduce to unique set
            mol = Chem.RemoveHs(mol)
            sortedcids = sorted(cids,key = lambda cid: cenergy[cid])
            selectedcids = []
            for conf in sortedcids:
                #set torsions to zero
                for m in rotmatches:
                    rdMolTransforms.SetDihedralRad(mol.GetConformer(conf),*m,value=0)
                #check rmsd
                for seenconf in selectedcids:
                    rms = getRMS(mol,seenconf,conf) 
                    if rms < args.rms_threshold:
                        break
                else: #loop completed normally - no break, included empty
                    selectedcids.append(conf)
                    
            #now exhaustively drive torsions of selected conformers
            if args.verbose:
                print len(selectedcids),"unique (ignoring torsions) starting conformers"
                
            total = 0
            for conf in selectedcids:
                total += genConformer_r(mol, conf, 0, rotmatches, args.degree, sdwriter)
            if args.verbose:
                print "%d total conformations generated"%total
                
        except (KeyboardInterrupt, SystemExit):
            raise                
        except Exception as e:
            print traceback.print_exc()
    else:
        print "ERROR:",smi

sdwriter.close()
outf.close()
