# -*- coding: utf-8 -*-
import argparse
from ChemClass import *

##########################################
## Options and defaults
##########################################
def getOptions():
    parser = argparse.ArgumentParser(description='python *.py [option]"')
    requiredgroup = parser.add_argument_group('required arguments')
    requiredgroup.add_argument('--react', dest='reacts', help='sdf file for reacts', default='', required=True)
    requiredgroup.add_argument('--product', dest='products', help='sdf file for products', default='', required=True)
    parser.add_argument('--nmx',dest='nmx',help='numer of top max similar reacts to record for each type(acid and amines)', type=int, default=1)
    parser.add_argument('--out',dest='out',help='output rmsd file', default='output.o')
    parser.add_argument('--sdfout',dest='sdfout',help='output sdf file', default='')


    args = parser.parse_args()

    return args

##########################################
## Master function
##########################################           
def main():
    options = getOptions()
    print "read two sdf files..."
    reactparse = ChemParse()
    reactparse.sdf_reader(options.reacts)
    productparse = ChemParse()
    productparse.sdf_reader(options.products)
    
    print "generate rmsd matrix..."
    psource = []
# ######################################################
#     #small sample test
#     for i in range(100):
#         psource.append(productparse.source.next())
#     rsource = []
#     for i in range(100):
#         rsource.append(reactparse.source.next())
#     rmsdcomp = RMSDCompare(psource, rsource)
# ######################################################
       
    rmsdcomp = RMSDCompare(productparse.source, reactparse.source)
    rmsdcomp.getRMSDmatrix()
    print "\nwrite to output..."
    rmsdcomp.writeMaxtrix(options.out,options.nmx)
    if options.sdfout:
        rmsdcomp.writesdf(options.sdfout,options.nmx)
    
    
        
    
if __name__ == "__main__":
    main()