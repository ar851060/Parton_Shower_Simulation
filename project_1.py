from random import random
import math
import numpy as np
import matplotlib.pyplot as plt

qcut = 0.2
b0 = (11 * 3 - 2 * 5) / (12 * math.pi)  # CA = 3, nf = 5
Cf = 4/3
const = Cf/(2*math.pi)
qhard = 1000
qmax = 2*qhard

# all the function we need
fun = lambda t,z: (1+pow(z,2))/(t*(1-z))
gun = lambda z: 2/(1-z)
Gunz = lambda z: -2*math.log(1-z)
InvGunz = lambda y: 1-math.exp(-y/2)
maket = lambda q: math.log(pow(q/0.2,2))
findq = lambda t: 0.2*math.sqrt(math.exp(t))
getrb = lambda ma,mb,mc: (ma ** 2 + mc ** 2 - mb ** 2 - math.sqrt((ma ** 2 - mc ** 2 - mb ** 2) ** 2 - 4 * (mc ** 2) * (mb ** 2))) / (2 * ma ** 2)
getrc = lambda ma,mb,mc: (ma ** 2 - mc ** 2 + mb ** 2 - math.sqrt((ma ** 2 - mc ** 2 - mb ** 2) ** 2 - 4 * (mc ** 2) * (mb ** 2))) / (2 * ma ** 2)


def main():
    tlist = []
    zlist = []
    for idx in range(10000):
        q1,z1 = veto(qmax)
        q2,z2 = veto(qmax)
        q1,q2,z1,z2, ene1, ene2 = checking(qmax,q1,q2,z1,z2)
        if z1 == 1 or z1 == 0:
            continue
        else:
            theta1 = (q1/(qhard*math.sqrt(z1*(1-z1))))
        tlist.append(q1)
        zlist.append(theta1)
        
    # create the graph of jet mass v.s. angel
    plt.xlabel("jet mass")
    plt.ylabel("angel")
    plt.title("with correction")
    plt.scatter(tlist,zlist)
    plt.show()
    print(sum(tlist)/(idx+1))
    
    # create the graph of jet mass v.s. number of events
    # step = 5
    # x = np.arange(0, 600, step)
    # plt.xlim(0 - step, 600 + step)
    # plt.hist(tlist, bins=x)
    # plt.xlabel("jet mass")
    # plt.ylabel("number of events")
    # plt.title("with correction")
    # plt.show()
    
    # create the graph of z v.s. number of events
    # step = 0.01
    # x = np.arange(0, 1, step)
    # plt.xlim(0 - step, 1 + step)
    # plt.hist(zlist, bins=x)
    # plt.xlabel("z")
    # plt.ylabel("number of events")
    # plt.title("with correction")
    # plt.show()


def zmake(q,e): # create the maximum and minimum of z
    const = 1-pow(q/e,2)
    if const <= 0:
        return 0.5,0.5
    else:
        zmin = 0.5*(1-math.sqrt(1-pow(q/e,2)))
        zmax = 0.5*(1+math.sqrt(1-pow(q/e,2)))
        return zmin, zmax

def veto(qp,zmin=0,zmax=0): # doing Veto Algorithm
    zmin0, zmax0 = zmake(qcut, qhard)
    alpha_0=alpha(qcut)
    Iz = (Gunz(zmax0) - Gunz(zmin0))*const*alpha_0
    while (True):
        ipt = 1 / Iz
        tempq_sqr = pow(qp,2) * pow(random(), ipt)
        if tempq_sqr > pow(qp,2):
            print("warning")
        if tempq_sqr < pow(qcut,2):
            return qcut,0
        mu = 8/3*alpha_0/(2*math.pi)
        # z = 1-(1-zmin0)*math.exp(-Iz*random()/mu)
        z = InvGunz(random() * (Gunz(zmax0) - Gunz(zmin0)) + Gunz(zmin0))
        tempq = math.sqrt(tempq_sqr)
        # tempt = maket(tempq)
        # q = findq(tempt)
        # ene = (pow(q, 2) + pow(qp, 2)) / (2 * qp)
        zmin,zmax = zmake(tempq,qhard)
        if z>zmax or z<zmin:
            qp = tempq
            continue
        alpha_s = alpha(tempq)
        # alpha_s =0.3
        splitting=(pow(zmin,2)-pow(zmax,2)+2*(zmin-zmax)+4*math.log((1-zmin)/(1-zmax)))*alpha_s*(4/3)*0.5/(2*math.pi)
        accept = splitting/(mu*math.log((1-zmin)/(1-zmax)))
        if random()>accept:
            qp = tempq
            continue
        else:
            return tempq,z

def checking(qp, q1, q2, z1, z2): # checking all variables in constrains
    while (True):
        if q1 + q2 > qp:
            if q1 > q2:
                q1, z1 = veto(q1)
                continue
            else:
                q2, z2 = veto(q2)
                continue
        rb = getrb(qp, q1, q2)
        rc = getrc(qp, q1, q2)
        pmom = math.sqrt((pow(qp,4)+pow(q1,4)+pow(q2,4)-2*(pow(qp,2)*pow(q1,2)+pow(qp,2)*pow(q2,2)+pow(q1,2)*pow(q2,2)))/(4*pow(qp,2)))
        ene1 = math.sqrt(pow(pmom,2)+pow(q1,2))
        ene2 = math.sqrt(pow(pmom,2)+pow(q2,2))
        # ene1 = (pow(qp,2)+pow(q1,2))/(2*qp)
        # ene2 = (pow(qp,2)+pow(q2,2))/(2*qp)
        eneb = ene1 - rb * ene1 + rc * ene2
        enec = ene2 + rb * ene1 - rc * ene2
        zminb, zmaxb = zmake(q1, eneb)
        zminc, zmaxc = zmake(q2, enec)
        if z2 == 0 or z1 == 0:
            break
        if q1 > q2:
            if z1 > zmaxb or z1 < zminb:
                q1, z1 = veto(q1)
                continue
            elif z2 > zmaxc or z2 < zminc:
                q2, z2 = veto(q2)
                continue
            elif q1 + q2 > qp:
                q1, z1 = veto(q1)
                continue
        else:
            if z2 > zmaxc or z2 < zminc:
                q2, z2 = veto(q2)
                continue
            elif z1 > zmaxb or z1 < zminb:
                if z1==0:
                    break
                q1, z1 = veto(q1)
                continue
            elif q1 + q2 > qp:
                q2, z2 = veto(q2)
                continue
        break
    return q1, q2, z1, z2, eneb, enec

def alpha(mu): #create alpha constant
    mz = 91.2
    alpha_z = 0.12
    return alpha_z/(1+2*b0*alpha_z*math.log(mu/mz))

main()
