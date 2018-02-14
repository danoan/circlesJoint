import numpy as np
import matplotlib.pyplot as plt

import utils as UT

D1 = 0
D2 = 4
W = 0.000001
f = 1e-2
k = UT.findAdequateK(D2,W,f)
s = UT.findAdequateS(D2)

wm = UT.wMapping(D2,W)
wrm = UT.wMappingReverse(D2,W)


def um(beta):
    return UT.uniformMapping(beta)

def G(theta,beta):
    pum = um(beta)
    return UT.heaviside(k,wm(pum(theta)))

def S(theta,beta):
    pum = um(beta)
    return UT.mirrorHeaviside(k,wrm(pum(theta)),s)


def ccX(theta,centerXList,centerYList,radiusList,betaList):
    np1 = len(centerXList)
    np2 = len(centerYList)
    np3 = len(radiusList)
    np4 = len(betaList)
    if not(np1==np2==np3==np4):
        print("Missing parameters")
        exit

    lp = list(zip(centerXList,centerYList,radiusList,betaList))
    pp = zip(lp,lp[1:]+[ lp[0] ])
    x = 0
    for t1,t2 in pp:
        (cx,cy,r,beta) = t1
        (ncx,ncy,nr,nbeta) = t2
        
        x += cx+r*np.cos(theta)*S(theta,beta)*G(theta,beta-np.pi/2) + ncx+nr*np.cos(theta)*S(theta,beta+np.pi/2)*G(theta,beta)

    return x

def ccY(theta,centerXList,centerYList,radiusList,betaList):
    np1 = len(centerXList)
    np2 = len(centerYList)
    np3 = len(radiusList)
    np4 = len(betaList)
    if not(np1==np2==np3==np4):
        print("Missing parameters")
        exit

    lp = list( zip(centerXList,centerYList,radiusList,betaList) )
    pp = zip(lp,lp[1:]+[ lp[0] ])
    y = 0
    for t1,t2 in pp:
        (cx,cy,r,beta) = t1
        (ncx,ncy,nr,nbeta) = t2
        
        y += cy+r*np.sin(theta)*S(theta,beta)*G(theta,beta-np.pi/2) + ncy+nr*np.sin(theta)*S(theta,beta+np.pi/2)*G(theta,beta)

    return y


def lengthDiff(circleJoint,l1,l2):
    domain = np.arange(0,2*np.pi,0.01)

    fx = lambda theta :     ccX(theta,
                                circleJoint.centerXList,
                                circleJoint.centerYList,
                                circleJoint.radiusList,
                                circleJoint.betaList)
    fy = lambda theta: ccY(theta,
                           circleJoint.centerXList,
                           circleJoint.centerYList,
                           circleJoint.radiusList,
                           circleJoint.betaList)

    rp = UT.rectParameterization(l1,l2)

    rx = lambda theta: rp(theta)[0]
    ry = lambda theta: rp(theta)[1]


    dist = lambda x1,y1,x2,y2: pow( pow(x2-x1,2)+pow(y2-y1,2),0.5 )

    imageDist = [ dist(fx(theta),
                       fy(theta),
                       rx(theta),
                       ry(theta)) for theta in domain ]

    imageRect = [ ry(theta) for theta in domain ]
    imageCJ = [ fy(theta) for theta in domain ]

    domainRect = [ rx(theta) for theta in domain ]
    domainCJ = [ fx(theta) for theta in domain ]

    fig = plt.figure(1)
    plt.subplot(211)
    plt.plot(domainRect,imageRect,'b-',domainCJ,imageCJ,'r-')

    plt.subplot(212)
    plt.plot(domain,imageDist,'g-')
    plt.show()

class CircleJoint:
    def __init__(self):
        self.centerXList = [-0.5,0.5]
        self.centerYList = [0,0]
        self.radiusList = [3,3]
        self.betaList = [np.pi/4.0,3*np.pi/4.0]


def main():
    circleJoint = CircleJoint()    
    lengthDiff(circleJoint,8,8)


if __name__=='__main__':
    main()
