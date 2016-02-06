# -*- coding: utf-8 -*-
import argparse
from ChemClass import *
import os

##########################################
## Options and defaults
##########################################
def getOptions():
    parser = argparse.ArgumentParser(description='python *.py [option]"')
    requiredgroup = parser.add_argument_group('required arguments')
    requiredgroup.add_argument('--react', dest='reacts', help='sdf file for reacts', default='', required=True)
    requiredgroup.add_argument('--product', dest='products', help='sdf file for products', default='', required=True)
    requiredgroup.add_argument('--job',dest='job',help='temple file to submit job', default='',required=True)
    parser.add_argument('--nmx',dest='nmx',help='numer of top max similar reacts to record for each type(acid and amines)', type=int, default=1)
    parser.add_argument('--split',dest='split',help='numer of comounds per input products.sdf ', type=int, default=500)



    args = parser.parse_args()

    return args

##########################################
## Master function
##########################################           
def main():
    options = getOptions()
    with open(options.job) as fjob:
        jobline = fjob.read()
    
    productparse = ChemParse()
    productparse.sdf_reader(options.products)
    allcompound = len(productparse.source)

    
    allcompound = len(productparse.source)
    productparse.source.reset()
    
    filecount = allcompound/options.split + 1
    
    num = 1
    for f in range(filecount):
        currentpath = os.getcwd()
        newfolder = "product_%d" % num
        try:
            os.mkdir(newfolder) 
        except:
            pass
        print newfolder
        newpath = os.path.join(currentpath,newfolder)
        newsdf = os.path.join(newpath,newfolder+'.sdf')
        
        sdfw = Chem.SDWriter(newsdf)
        for count in range(options.split):
            try:
                m = productparse.source.next()
                sdfw.write(m)
            except:
                break

        sdfw.close()
        
        newpbs = os.path.join(newpath,newfolder+'.pbs')
        pbsw = open(newpbs,'w')
        pbsw.write("#PBS -N rmsdRun%d\n" % num)
        pbsw.write(jobline+'\n')
        react = os.path.join(currentpath,'4_acid_amine_torsionfilter.sdf')
        out = os.path.join(newpath,'output.o')
        pbsw.write('python mainserial.py --react=%s --product=%s --out=%s\n' % (react, newsdf, out))
        pbsw.close()
        num += 1
    

        
    
    
    
        
    
if __name__ == "__main__":
    main()