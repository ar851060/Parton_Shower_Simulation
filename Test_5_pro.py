#try to simulate the parton shower

import math
import numpy as np
import matplotlib.pyplot as plt


massi = 0
massf = 0
massr = 0

enini = 90
pxini = math.sqrt(1/2.0)*enini
pyini = math.sqrt(1/2.0)*enini
xini = 0.0 # math.sqrt(1/2.0)*enini
yini = 0.0 #pyini
v = 1
#plt.scatter(0,0)
plt.scatter(xini,yini)
#plt.plot((0,xini),(0,yini))

PDF = lambda x : 1/(x)
CDF = lambda x : np.log(x)
invCDF = lambda y : np.exp(y)

zmin = pow(10,-4)
zmax = 1.

thetamin = pow(10, -4)
thetamax = math.pi/2

N = 1
enthr = 1

yzmin = CDF(zmin)
yzmax = CDF(zmax)

ythetamax = CDF(thetamax)
ythetamin = CDF(thetamin)

unsplten = []
unsplten.append(enini)
unspltpx = []
unspltpx.append(pxini)
unspltpy = []
unspltpy.append(pyini)
unspltx = []
unspltx.append(xini)
unsplty = []
unsplty.append(yini)

unsplten.append(enini)
unspltpx.append(-pxini)
unspltpy.append(-pyini)
unspltx.append(xini)
unsplty.append(yini)

total_x = 0
total_y = 0

while(len(unsplten)!=0):
    # should this particle be emmitting? 
    emmission_prob = (0.35 * math.log(enini/enthr))
    rndm_emmits = np.random.uniform(0,1,1)

    #print "emmission prob: %lf and rndm_emmits: %lf" % (emmission_prob, rndm_emmits)

    if (rndm_emmits > emmission_prob): # this gluon didn't emmit so we go to the next unsplitted guy
       #plt.plot((0,pxini/enini),(0,pyini/enini))
       total_x += pxini
       total_y += pyini
       #print "Total_px = %lf Total_py = %lf" % (total_x, total_y)
       unsplten.pop(0)
       unspltpx.pop(0)
       unspltpy.pop(0)
       unspltx.pop(0)
       unsplty.pop(0)
       if len(unsplten) != 0:
           enini = unsplten[0]
           pxini = unspltpx[0]
           pyini = unspltpy[0]
           xini = unspltx[0]
           yini = unsplty[0]
       continue 
      
    
    YZ = np.random.uniform(yzmin,yzmax,N)
    YTheta = np.random.uniform(ythetamin,ythetamax,N)

    sing = np.random.uniform(-1,1,1)
    Z = invCDF(YZ)
    Theta = invCDF(YTheta) * sing / abs(sing)

    enrad = enini*Z
    enfin = enini*(1-Z)

    #print "Enrad: %lf Theta: %lf, Enrad*Theta: %lf" % (enrad,Theta,math.fabs(enrad*Theta))

    while ( math.fabs(enrad*(Theta)) < enthr): #we already decided that this guy emmits, so we should be in the right range now
        YZ = np.random.uniform(yzmin,yzmax,N)
        YTheta = np.random.uniform(ythetamin,ythetamax,N)

        sing = np.random.uniform(-1,1,1)
        Z = invCDF(YZ)
        Theta = invCDF(YTheta)*sing/abs(sing)

        enrad = enini*Z
        enfin = enini*(1-Z)
        continue
    #massless particle approximation
    prad = enrad 

    #The tagged system is defined such that the direction of the incoming parton is x'
    pxrad_tag = prad*math.cos(Theta)
    pyrad_tag = prad*math.sin(Theta)

    pyfin_tag = -pyrad_tag #conservation of momentum
    pxfin_tag = enini - pxrad_tag

    reco_energy_rad = math.sqrt(pxrad_tag*pxrad_tag + pyrad_tag*pyrad_tag)
    reco_energy_fin = math.sqrt(pxfin_tag*pxfin_tag + pyfin_tag*pyfin_tag)
    print ("eini = %lf, Z = %lf, Theta= %lf, pt = %lf, reco_erad = %lf, reco_efin = %lf" % (enini, Z, Theta, (abs(Theta))*enini*Z, reco_energy_rad, reco_energy_fin))
#    alpha = math.atan(pyfin_tag/pxfin_tag)#(pyini-prad*math.sin(Theta))/pfin

    rot_angle =  math.atan2(pyini,pxini)
    xrad = pxrad_tag * math.cos(rot_angle) - pyrad_tag * math.sin(rot_angle) + xini
    yrad = pyrad_tag * math.cos(rot_angle) + pxrad_tag * math.sin(rot_angle) + yini
    xfin = pxfin_tag * math.cos(rot_angle) - pyfin_tag * math.sin(rot_angle) + xini
    yfin = pyfin_tag * math.cos(rot_angle) + pxfin_tag * math.sin(rot_angle) + yini

    pxrad = pxrad_tag * math.cos(rot_angle) - pyrad_tag * math.sin(rot_angle)
    pyrad = pyrad_tag * math.cos(rot_angle) + pxrad_tag * math.sin(rot_angle)
    pxfin = pxfin_tag * math.cos(rot_angle) - pyfin_tag * math.sin(rot_angle)
    pyfin = pyfin_tag * math.cos(rot_angle) + pxfin_tag * math.sin(rot_angle)

    plt.scatter(xrad,yrad)
    plt.scatter(xfin,yfin)
    plt.plot((xini,xrad),(yini,yrad))
    plt.plot((xini,xfin),(yini,yfin))

    #print "initial : xini =  %lf, yini = %lf, enini = %lf pxini = %lf pyini = %lf " % (xini, yini, enini, pxini, pyini)
    #print "rad : xrad =  %lf, yrad = %lf, enrad = %lf pxrad = %lf pyrad = %lf " % (xrad, yrad, enrad, pxrad, pyrad)
    #print "fin : xfin =  %lf, yfin = %lf, enfin = %lf pxfin = %lf pyfin = %lf " % (xfin, yfin, enfin, pxfin, pyfin)
    enini = enfin
    pxini = pxfin
    pyini = pyfin
    xini = xfin
    yini = yfin

    unsplten.append(enrad)
    unspltpx.append(pxrad)
    unspltpy.append(pyrad)
    unspltx.append(xrad)
    unsplty.append(yrad)
    continue

plt.show()
