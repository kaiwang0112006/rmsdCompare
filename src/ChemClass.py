# -*- coding: utf-8 -*-
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys,rdMolAlign
from rdkit.SimDivFilters import rdSimDivPickers
from rdkit.Chem.Fingerprints import FingerprintMols
import numpy
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina,Murtagh
from rdkit.ML.Cluster import ClusterUtils
from progressorsClass import *
import multiprocessing
import random
import string

#read sdf file and do fingerprint calculation AND clustering
class ChemParse(object):
    def __init__(self):
        self.source = ''
        self.fps = ''

        
    def sdf_reader(self,input):
        self.source = Chem.SDMolSupplier(input)
        
    def get_fps(self,ftype,radius=1):
        if self.source == '' or self.source == None:
            pass
        else:
            if ftype == 'mo':
                fps=[AllChem.GetMorganFingerprint(m,radius,useFeatures=True) for m in self.source]
            elif ftype == 'tp':
                fps = [FingerprintMols.FingerprintMol(m) for m in self.source]
            elif ftype == 'mc':
                fps = [MACCSkeys.GenMACCSKeys(m) for m in self.source]
        self.fps = fps

        
    def clusterOutput(self, output, cdict):
        sdfout = Chem.SDWriter(output)
        
        for index, m in enumerate(self.source):
            classid = cdict[index][0]
            isCentroid = cdict[index][1]
            m.SetProp("class",str(classid))
            m.SetProp("isCentroid",isCentroid)
            sdfout.write(m)
        sdfout.close()
        
class RMSDCompare(object):
    def __init__(self,products,reacts):
        self.products = products
        self.reacts = reacts
        self.rmatrix = {}
    
    #the matrix will be like:{product:{'acid':{'react1':0.2,'react2':0.3...},'amines':{'react1':0.2,'react2':0.3...}}}
    #hard code, cannot be reused in other case
    def getRMSDmatrix(self):
        total = len(self.products) * len(self.reacts)
        pb = progressbar(total, "*")
        for indexpdt, mp in enumerate(self.products):
            pdtname = mp.GetProp('_Name')
            self.rmatrix[pdtname] = {'acid':{},'amines':{}}
            for indexreact, mr in enumerate(self.reacts):
                reactname = mr.GetProp('_Name')
                try:
                    rmsd = AllChem.GetBestRMS(mp, mr)/(mp.GetNumAtoms()+mr.GetNumAtoms())
                except:
                    rmsd = -1

                if 'acid' in reactname:
                    self.rmatrix[pdtname]['acid'][reactname] = rmsd
                elif 'amines' in reactname:
                    self.rmatrix[pdtname]['amines'][reactname] = rmsd
                pb.progress(indexpdt*indexreact)
                
    def getRMSDmatrixparallel(self,np):
        print ('Creating %d process' % np)

        results = multiprocessing.Queue()
        
        process = []  
        jobs = []
        for indexpdt, mp in enumerate(self.products):
            for indexreact, mr in enumerate(self.reacts):
                jobs.append((mp,mr,mp.GetProp('_Name'),mr.GetProp('_Name')))
        for i in range(np):  
            process.append(multiprocessing.Process(target=self.getRMSDmatrixTask, args=(jobs[i::np],results)))  

        # Start processes one by one  

        for p in process:  
            p.start()  

        # Wait for all processed to finish  
        for p in process:  
            p.join()  
             
        for indexpdt, mp in enumerate(self.products):
            for indexreact, mr in enumerate(self.reacts):
                r = results.get()
                pdtname,reactname,rmsd = r
                if not pdtname in self.rmatrix:
                    self.rmatrix[pdtname] = {'acid':{},'amines':{}}
                if 'acid' in reactname:
                    self.rmatrix[pdtname]['acid'][reactname] = rmsd
                elif 'amines' in reactname:
                    self.rmatrix[pdtname]['amines'][reactname] = rmsd

        print "multiprocessing done!"
        
    def getRMSDmatrixTask(self,joblist,result_queue):
        #r = []
        for i,compoundPair in enumerate(joblist):
            mp = compoundPair[0]
            mr = compoundPair[1]
            pdtname = compoundPair[2]
            reactname = compoundPair[3]
            try:
                rmsd = AllChem.GetBestRMS(mp, mr)/(mp.GetNumAtoms()+mr.GetNumAtoms())
            except:
                rmsd = 0
            #r.append([pdtname,reactname,rmsd])
            #print [pdtname,reactname,rmsd]
            result_queue.put([pdtname,reactname,rmsd])
        

        
                       
    def writeMaxtrix(self,output,nmx):
        fout = open(output,'w')
        for pdtname in self.rmatrix:
            acidNvalues = self.getNcompound(self.rmatrix[pdtname]['acid'],nmx)
            aminesNvalues = self.getNcompound(self.rmatrix[pdtname]['amines'],nmx)
            line = pdtname +"\t"
            for v in acidNvalues:
                line += str(v) + '\t'
            for v in aminesNvalues:
                line += str(v) + '\t'
                
            fout.write(line[:-1]+'\n')
        fout.close()
            
            
    def getNcompound(self,values,nmx):
        # sort dict value to like: [('33', 56), ('a', 31), ('bc', 5), ('asd', 4), ('c', 3), ('d', 0)]
        sortvalues = sorted(values.iteritems(), key=lambda d:d[1], reverse = True )
        if nmx > len(sortvalues):
            nmx = len(sortvalues)
        return [sortvalues[i] for i in range(nmx)]
        
class Fingerprint_Cluster(object):
    def __init__( self, fps):
        self.fplist = fps
        self.dist = []
        self.cdict = {} 
        self.clustdict = {}

    #generate the distance matrix
    def distance_matrix(self):
        self.dist = []
        nfps = len(self.fplist)
        for i in range(1,nfps):
            sims = DataStructs.BulkTanimotoSimilarity(self.fplist[i],self.fplist[:i])
            self.dist.extend([1-x for x in sims])
        
    #generate cluster dict as {1:[1,2,3],2:[4,5,6]...}
    def cluster_dict(self,algorithm, cutoff=0.5, method='Wards', ncluster=1):
        if algorithm == 'Butina':
            self.ClusterFps_Butina(self.dist,len(self.fplist),cutoff)
        elif algorithm == 'Murtagh':
            self.ClusterFps_Murtagh(self.dist,len(self.fplist),method,ncluster)
        
    def ClusterFps_Butina(self, dists, nfps,cutoff):
        self.cdict = {}
        cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
        for index, eachcs in enumerate(cs): 
            self.clustdict[index+1] = eachcs
            for eachid in eachcs:
                self.cdict[eachid] = [index+1]
                if eachid == eachcs[0]:
                    self.cdict[eachid].append("true")
                else:
                    self.cdict[eachid].append("flase")
            
    def ClusterFps_Murtagh(self, dists, nfps, method, ncluster):
        self.cdict = {}
        cs = None
        if method == 'Wards':
            cs= Murtagh.ClusterData(dists,len(self.fplist),Murtagh.WARDS,isDistData=1)
        elif method == 'SLINK':
            cs= Murtagh.ClusterData(dists,len(self.fplist),Murtagh.SLINK,isDistData=1)
        elif method == 'CLINK':
            cs= Murtagh.ClusterData(dists,len(self.fplist),Murtagh.CLINK,isDistData=1)
        elif method == 'UPGMA':
            cs= Murtagh.ClusterData(dists,len(self.fplist),Murtagh.UPGMA,isDistData=1)
            
        splitClusts=ClusterUtils.SplitIntoNClusters(cs[0],ncluster)
        #centroids = [ClusterUtils.FindClusterCentroidFromDists(x,dists) for x in splitClusts]
        for index, cluster in enumerate(splitClusts):
            children = cluster.GetPoints() 
            pts = [x.GetData() for x in children]  
            self.clustdict[index+1] = pts 
            for pt in pts:
                self.cdict[pt] = [index + 1] 
                if pt == pts[0]:
                    self.cdict[pt].append("true")
                else:
                    self.cdict[pt].append("flase")
                    

                    