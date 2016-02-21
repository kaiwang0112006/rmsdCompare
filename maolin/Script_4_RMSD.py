#!/usr/bin/env python
#############################################################################
# Copyright (C) 2014 OpenEye Scientific Software, Inc.
#############################################################################
# Align two compounds based on smarts match
#############################################################################
from __future__ import print_function
import sys
from openeye.oechem import *
import math

def SmartsAlign(acidmol, aminemol, ssacid, ssamine):
    unique = True
    for match1 in ssacid.Match(acidmol,unique):
        for match2 in ssamine.Match(aminemol,unique):
            match = OEMatch()
            for mp1, mp2 in zip(match1.GetAtoms(), match2.GetAtoms()):
                match.AddPair(mp1.target, mp2.target)
            overlay = False 
            return  OERMSD(acidmol, aminemol, match, overlay)
def CoordstoSpherical(x,y,z):
    r=math.sqrt(math.pow(x,2)+math.pow(y,2)+math.pow(z,2))
    sita=math.acos(z/r)
    if x!=0:
        faiabs=abs(math.atan(y/x))
    if x==0 and y>0:
        fai=math.pi*0.5
    if x==0 and y<0:
        fai=math.pi*1.5  
    if x>0 and y>0:
        fai=faiabs
    if x<0 and y>0:
        fai=math.pi-faiabs
    if x<0 and y<0:
        fai=math.pi+faiabs
    if x>0 and y<0:
        fai=2*math.pi-faiabs
    if y==0:
        fai=0.0 
    return (r,sita,fai)

def SphericaltoCoords(r,sita,fai):
    x=r*math.sin(sita)*math.cos(fai)
    y=r*math.sin(sita)*math.sin(fai)
    z=r*math.cos(sita)
    return (x,y,z)

def ProjectionCalculate(vector, normalvector):
    '''The origin must be on the normal plane
       x0-x=aA, y0-y=aB, z0-z=aC, Ax+By+Cz=0
       where x0,y0,z0 is from vector, A,B,C is from normalvector'''
    x0=vector[0]
    y0=vector[1]
    z0=vector[2]
    A=normalvector[0]
    B=normalvector[1]
    C=normalvector[2]
    a=(x0*A+y0*B+z0*C)/(math.pow(A,2)+math.pow(B,2)+math.pow(C,2))
    x=x0-a*A
    y=y0-a*B
    z=z0-a*C
    return (x,y,z)

def CrossProduct(vector1,vector2):
    x=vector1[1]*vector2[2]-vector1[2]*vector2[1]
    y=vector1[2]*vector2[0]-vector1[0]*vector2[2]
    z=vector1[0]*vector2[1]-vector1[1]*vector2[0]
    return (x,y,z)

def AngleCalculation(vector1,vector2):
    a=[vector1[0],vector1[1],vector1[2]]
    b=[vector2[0],vector2[1],vector2[2]]
    ab=a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
    am=math.sqrt(math.pow(a[0],2)+math.pow(a[1],2)+math.pow(a[2],2))
    bm=math.sqrt(math.pow(b[0],2)+math.pow(b[1],2)+math.pow(b[2],2))
    sita=math.acos((ab)/(am*bm))
    return sita

def Torsionnew(acidmol, aminemol, ssacid, ssamine):
    unique = True
    acidcoords=OEFloatArray(acidmol.GetMaxAtomIdx()*3)
    acidmol.GetCoords(acidcoords)
    aminecoords=OEFloatArray(aminemol.GetMaxAtomIdx()*3)
    aminemol.GetCoords(aminecoords)
    acididxlist=[]
    amineidxlist=[]
    for match1 in ssacid.Match(acidmol,unique):
        for match2 in ssamine.Match(aminemol,unique):
            for acidatom in match1.GetTargetAtoms():
                acididxlist.append(acidatom.GetIdx())
                for amineatom in match2.GetTargetAtoms():
                    amineidxlist.append(amineatom.GetIdx())
    for match1 in ssacid.Match(acidmol,unique):
        for match2 in ssamine.Match(aminemol,unique):
            for acidatom in match1.GetTargetAtoms():
                acididx=acidatom.GetIdx()
                acidsyb=OEGetAtomicSymbol(acidatom.GetAtomicNum())
                if acidsyb=='C':
                    acidcoordsC=[acidcoords[acididx*3],acidcoords[acididx*3+1],acidcoords[acididx*3+2]]
                    for nbor in acidatom.GetAtoms():
                        if nbor.GetIdx() not in acididxlist:
                            acidcoordsR=[acidcoords[nbor.GetIdx()*3],acidcoords[nbor.GetIdx()*3+1],acidcoords[nbor.GetIdx()*3+2]]
                    for amineatom in match2.GetTargetAtoms():
                        amineidx=amineatom.GetIdx()
                        aminesyb=OEGetAtomicSymbol(amineatom.GetAtomicNum())
                        if aminesyb=='N':
                            for nbor in amineatom.GetAtoms():
                                if nbor.GetIdx() not in amineidxlist:
                                    if OEGetAtomicSymbol(nbor.GetAtomicNum())!='H':
                                        aminecoordsR=[aminecoords[nbor.GetIdx()*3],
                                                      aminecoords[nbor.GetIdx()*3+1],
                                                      aminecoords[nbor.GetIdx()*3+2]]
                if acidsyb=='N':
                    acidcoordsN=[acidcoords[acididx*3],acidcoords[acididx*3+1],acidcoords[acididx*3+2]]
                if acidsyb=='O':
                    acidcoordsO=[acidcoords[acididx*3],acidcoords[acididx*3+1],acidcoords[acididx*3+2]]
    distanceRR=CalculateDistance(aminecoordsR,acidcoordsR)
    distanceRC=CalculateDistance(aminecoordsR,acidcoordsC)
    distanceRN=CalculateDistance(aminecoordsR,acidcoordsN)
    #print ("%f,%f,%f"%(distanceRR,distanceRC,distanceRN))
    #print ("Acid C: %f,%f,%f\n"%(acidcoordsC[0],acidcoordsC[1],acidcoordsC[2]))
    #print ("Acid N: %f,%f,%f\n"%(acidcoordsN[0],acidcoordsN[1],acidcoordsN[2]))
    #print ("Acid R: %f,%f,%f\n"%(acidcoordsR[0],acidcoordsR[1],acidcoordsR[2]))
    #print ("Amine R: %f,%f,%f\n"%(aminecoordsR[0],aminecoordsR[1],aminecoordsR[2]))

    if distanceRR>distanceRC>distanceRN:
        return 1
    else:
        return 0

def CalculateDistance(alist,blist):
    distance=math.sqrt((alist[0]-blist[0])**2+(alist[1]-blist[1])**2+(alist[2]-blist[2])**2)
    return distance

def TorsionAngleForCON(acidmol, aminemol, ssacid, ssamine):
    unique = True
    acidcoords=OEFloatArray(acidmol.GetMaxAtomIdx()*3)
    acidmol.GetCoords(acidcoords)
    aminecoords=OEFloatArray(aminemol.GetMaxAtomIdx()*3)
    aminemol.GetCoords(aminecoords)
    acididxlist=[]
    amineidxlist=[]
    for match1 in ssacid.Match(acidmol,unique):
        for match2 in ssamine.Match(aminemol,unique):
            for acidatom in match1.GetTargetAtoms():
                acididxlist.append(acidatom.GetIdx())
                for amineatom in match2.GetTargetAtoms():
                    amineidxlist.append(amineatom.GetIdx())
    for match1 in ssacid.Match(acidmol,unique):
        for match2 in ssamine.Match(aminemol,unique):
            for acidatom in match1.GetTargetAtoms():
                acididx=acidatom.GetIdx()
                acidsyb=OEGetAtomicSymbol(acidatom.GetAtomicNum())
                if acidsyb=='C':
                    acidcoordsC=[acidcoords[acididx*3],acidcoords[acididx*3+1],acidcoords[acididx*3+2]]
                    for nbor in acidatom.GetAtoms():
                        if nbor.GetIdx() not in acididxlist:
                            acidcoordsR=[acidcoords[nbor.GetIdx()*3],acidcoords[nbor.GetIdx()*3+1],acidcoords[nbor.GetIdx()*3+2]]
                    for amineatom in match2.GetTargetAtoms():
                        amineidx=amineatom.GetIdx()
                        aminesyb=OEGetAtomicSymbol(amineatom.GetAtomicNum())
                        if aminesyb=='C':
                            Dvalues=[acidcoords[acididx*3]-aminecoords[amineidx*3],
                                     acidcoords[acididx*3+1]-aminecoords[amineidx*3+1],
                                     acidcoords[acididx*3+2]-aminecoords[amineidx*3+2]]
                    for amineatom in match2.GetTargetAtoms():
                        amineidx=amineatom.GetIdx()
                        aminesyb=OEGetAtomicSymbol(amineatom.GetAtomicNum())
                        if aminesyb=='N':
                            aminecoordsnewN=[aminecoords[amineidx*3]+Dvalues[0],
                                             aminecoords[amineidx*3+1]+Dvalues[1],
                                             aminecoords[amineidx*3+2]+Dvalues[2]]
                            for nbor in amineatom.GetAtoms():
                                if nbor.GetIdx() not in amineidxlist:
                                    if OEGetAtomicSymbol(nbor.GetAtomicNum())!='H':
                                        aminecoordsR=[aminecoords[nbor.GetIdx()*3],
                                                      aminecoords[nbor.GetIdx()*3+1],
                                                      aminecoords[nbor.GetIdx()*3+2]]
                                        aminecoordsnewR=[aminecoords[nbor.GetIdx()*3]+Dvalues[0],
                                                         aminecoords[nbor.GetIdx()*3+1]+Dvalues[1],
                                                         aminecoords[nbor.GetIdx()*3+2]+Dvalues[2]]
                        if aminesyb=='O':
                            aminecoordsnewO=[aminecoords[amineidx*3]+Dvalues[0],
                                             aminecoords[amineidx*3+1]+Dvalues[1],
                                             aminecoords[amineidx*3+2]+Dvalues[2]]
                if acidsyb=='N':
                    acidcoordsN=[acidcoords[acididx*3],acidcoords[acididx*3+1],acidcoords[acididx*3+2]]
                if acidsyb=='O':
                    acidcoordsO=[acidcoords[acididx*3],acidcoords[acididx*3+1],acidcoords[acididx*3+2]]
    #Step 1: C of amine move to C of acid
    acidcoords0R=[acidcoordsR[0]-acidcoordsC[0],acidcoordsR[1]-acidcoordsC[1],acidcoordsR[2]-acidcoordsC[2]]
    acidcoords0N=[acidcoordsN[0]-acidcoordsC[0],acidcoordsN[1]-acidcoordsC[1],acidcoordsN[2]-acidcoordsC[2]]
    acidcoords0O=[acidcoordsO[0]-acidcoordsC[0],acidcoordsO[1]-acidcoordsC[1],acidcoordsO[2]-acidcoordsC[2]]
    aminecoords0R=[aminecoordsnewR[0]-acidcoordsC[0],aminecoordsnewR[1]-acidcoordsC[1],aminecoordsnewR[2]-acidcoordsC[2]]
    aminecoords0N=[aminecoordsnewN[0]-acidcoordsC[0],aminecoordsnewN[1]-acidcoordsC[1],aminecoordsnewN[2]-acidcoordsC[2]]
    aminecoords0O=[aminecoordsnewO[0]-acidcoordsC[0],aminecoordsnewO[1]-acidcoordsC[1],aminecoordsnewO[2]-acidcoordsC[2]]
    list1=[acidcoords0R,acidcoords0N,acidcoords0O,aminecoords0R,aminecoords0N,aminecoords0O]
    acidsphericalR=[]
    acidsphericalN=[]
    acidsphericalO=[]
    aminesphericals1R=[]
    aminesphericals1N=[]
    aminesphericals1O=[]
    list2=[acidsphericalR,acidsphericalN,acidsphericalO,aminesphericals1R,aminesphericals1N,aminesphericals1O]
    for i in range(len(list1)):
        listtemp=list(CoordstoSpherical(list1[i][0],list1[i][1],list1[i][2]))
        #r=math.sqrt(math.pow(list1[i][0],2)+math.pow(list1[i][1],2)+math.pow(list1[i][2],2))
        #sita=math.acos(list1[i][2]/r)
        #fai=math.atan(list1[i][1]/list1[i][0])
        list2[i].append(listtemp[0])
        list2[i].append(listtemp[1])
        list2[i].append(listtemp[2])
    #end of step 1
    #Step 2: N of amine move to N of acid
    TorsionDvalueN=[0,acidsphericalN[1]-aminesphericals1N[1],acidsphericalN[2]-aminesphericals1N[2]]
    aminesphericals2N=[aminesphericals1N[0]+TorsionDvalueN[0],aminesphericals1N[1]+TorsionDvalueN[1],aminesphericals1N[2]+TorsionDvalueN[2]]
    aminesphericals2R=[aminesphericals1R[0]+TorsionDvalueN[0],aminesphericals1R[1]+TorsionDvalueN[1],aminesphericals1R[2]+TorsionDvalueN[2]]
    aminesphericals2O=[aminesphericals1O[0]+TorsionDvalueN[0],aminesphericals1O[1]+TorsionDvalueN[1],aminesphericals1O[2]+TorsionDvalueN[2]]
    #end of step 2
    #Step 3: SphericaltoCoords
    aminescoords02R=list(SphericaltoCoords(aminesphericals2R[0],aminesphericals2R[1],aminesphericals2R[2]))
    aminecoords02N=list(SphericaltoCoords(aminesphericals2N[0],aminesphericals2N[1],aminesphericals2N[2]))
    acidcoords0Ntest=list(SphericaltoCoords(acidsphericalN[0],acidsphericalN[1],acidsphericalN[2]))
    aminescoords02O=list(SphericaltoCoords(aminesphericals2O[0],aminesphericals2O[1],aminesphericals2O[2]))
    #end of step 3
    #Step 4: Projection on the normal plane of C-N
    normalvector=acidcoords0N
    acidprojectionR=list(ProjectionCalculate(acidcoords0R, normalvector))
    acidprojectionO=list(ProjectionCalculate(acidcoords0O, normalvector))
    aminescoords02Rorigin=[aminescoords02R[0]-acidcoords0N[0],
                           aminescoords02R[1]-acidcoords0N[1],
                           aminescoords02R[2]-acidcoords0N[2]]
    amineprojectionR=list(ProjectionCalculate(aminescoords02Rorigin, normalvector))
    amineprojectionO=list(ProjectionCalculate(aminescoords02O, normalvector))
    #end of step 4
    #Step 5: Cross product
    angledirectionO=list(CrossProduct(amineprojectionO,acidprojectionO))
    angledirectionR=list(CrossProduct(amineprojectionR,acidprojectionR))
    angleO=AngleCalculation(amineprojectionO,acidprojectionO)
    angleR=AngleCalculation(amineprojectionR,acidprojectionR)
    direction=angledirectionO[0]*angledirectionR[0]+angledirectionO[1]*angledirectionR[1]+angledirectionO[2]*angledirectionR[2]
    if direction>0:
        if abs(angleO-angleR)<=math.pi:
            torsionangleradian=abs(angleO-angleR)
        if abs(angleO-angleR)>math.pi:
            torsionangleradian=2*math.pi-abs(angleO-angleR)
    if direction<0:
        if abs(angleO+angleR)<=math.pi:
            torsionangleradian=abs(angleO+angleR)
        if abs(angleO+angleR)>math.pi:
            torsionangleradian=2*math.pi-abs(angleO+angleR)
    torsionangle=torsionangleradian/math.pi*180
    return (torsionangle,angleO,angleR)

def main(argv=[__name__]):
    if len(argv) != 3:
        OEThrow.Usage("%s <acid> <amine>" % argv[0])
    smartsacid = "[NX3;H1;!R;$(N[CX4H3])][CX3;!R](=[OX1])"
    smartsamine = "[NX3;H1;!R][CX3;$(C[CX4H3])](=[OX1])"
    outputname="/home/xh276/ucc-fileserver/1/12112015/RMSD/acid_amine.sdf"
    outputname_torsionfilter="/home/xh276/ucc-fileserver/1/12112015/RMSD/acid_amine_torsionfilter.sdf"
    cutoff =float(0.8)

    ofs = oemolostream()
    if not ofs.open(outputname):
        OEThrow.Fatal("Unable to create %s" % outputname)
    ofstorsion = oemolostream()
    if not ofstorsion.open(outputname_torsionfilter):
        OEThrow.Fatal("Unable to create %s" % outputname_torsionfilter)    
    acid = oemolistream()
    if not acid.open(argv[1]):
        OEThrow.Fatal("Unable to open %s for reading" % argv[1])

    amine = oemolistream()
    if not amine.open(argv[2]):
        OEThrow.Fatal("Unable to open %s for reading" % argv[2])

    ssacid = OESubSearch()
    if not ssacid.Init(smartsacid):
        OEThrow.Fatal("Unable to parse SMARTS: %s" % smarts)
    ssamine = OESubSearch()
    if not ssamine.Init(smartsamine):
        OEThrow.Fatal("Unable to parse SMARTS: %s" % smarts)
    acidnum=0
    aminenum=0
    for acidmol in acid.GetOEGraphMols():
        if not acidmol.GetDimension() == 3:
            OEThrow.Warning("%s doesn't have 3D coordinates" % acidmol.GetTitle())
            continue
        OEPrepareSearch(acidmol, ssacid)
        if not ssacid.SingleMatch(acidmol):
            OEThrow.Warning("SMARTS fails to match fitmol %s" % acidmol.GetTitle())
            continue

        acidnum+=1
        print('Progressing Acid %d'% acidnum)
        for aminemol in amine.GetOEGraphMols():
            if not aminemol.GetDimension() == 3:
                OEThrow.Warning("%s doesn't have 3D coordinates" % aminemol.GetTitle())
                continue
            OEPrepareSearch(aminemol, ssamine)
            if not ssamine.SingleMatch(aminemol):
                OEThrow.Warning("SMARTS fails to match fitmol %s" % aminemol.GetTitle())
                continue
            aminenum+=1
            valueRMSD=SmartsAlign(acidmol, aminemol, ssacid, ssamine)
            if valueRMSD < 0.5:
                if Torsionnew(acidmol, aminemol, ssacid, ssamine):
                    OESetSDData(acidmol, "RMSD", str(valueRMSD))
                    OESetSDData(aminemol, "RMSD", str(valueRMSD))
                    OEWriteMolecule(ofs, acidmol)
                    OEWriteMolecule(ofs, aminemol)
                    torsionangle=TorsionAngleForCON(acidmol, aminemol, ssacid, ssamine)
                    if (torsionangle[2]/math.pi*180)>=60:
                        OESetSDData(acidmol, "Torsion", str(torsionangle[0]))
                        OESetSDData(aminemol, "Torsion", str(torsionangle[0]))
                        OESetSDData(acidmol, "TorsionO", str(torsionangle[1]/math.pi*180))
                        OESetSDData(aminemol, "TorsionO", str(torsionangle[1]/math.pi*180))
                        OESetSDData(acidmol, "TorsionR", str(torsionangle[2]/math.pi*180))
                        OESetSDData(aminemol, "TorsionR", str(torsionangle[2]/math.pi*180))
                        OEWriteMolecule(ofstorsion, acidmol)
                        OEWriteMolecule(ofstorsion, aminemol)
                
        aminenum=0
        amine = oemolistream()
        amine.open(argv[2])

if __name__ == "__main__":
    sys.exit(main(sys.argv))
