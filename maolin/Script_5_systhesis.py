#!/usr/bin/env python
from maolin import *

if len(sys.argv) != 3:
    OEThrow.Usage("%s <acid_amine> <acid_amine_sys>" % sys.argv[0])
acid_amine_file=sys.argv[1]
outputfile=sys.argv[2]
smirks='[*:1][C;X3;!R:2](=[O:5])[N;X3;!R]([H])[C;X4]([H])([H])[H].[*:3][NH:4][C;!R](=[O])[C;X4]([H])([H])[H]>>[*:1][C;X3:2](=[O:5])[NH:4][*:3]'

ifs=OeReadfile(acid_amine_file)
ofs=OeWritefile(outputfile)
i=1
listtemp=[]
for mol in ifs.GetOEGraphMols():
    if i%2 !=0:
        listtemp.append(OEGraphMol(mol))
    if i%2 ==0:
        listtemp.append(OEGraphMol(mol))
        libgen = OELibraryGen(smirks)
        libgen.SetStartingMaterial(listtemp[0], 0)
        libgen.SetStartingMaterial(listtemp[1], 1)
        for molsys in libgen.GetProducts():
            OEWriteMolecule(ofs, molsys)
            break
        listtemp=[]
    i+=1
    

