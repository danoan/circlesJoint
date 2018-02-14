#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import newton

import utils as UT

def rectParameterization(l1,l2):
    '''
    Returns a function which receives one parameter theta
    in the range [0,2pi] and for each theta it returns
    a point on a rect of sides l1 and l2.

    -PI_4 <= t < PI_4 ==> (L1,L2*tan(t))
    PI_4 <= t < 3PI_4 ==> (L1/tan(t),L2)
    3PI_4 <= t < 5PI_4 ==> (-L1,L2*tan(t))
    5PI_4 <= t <7PI_4 ==> (L1/tan(t),-L2)

    I can simplify if I use t+PI_4 as parameter.
    '''
    L1 = l1/2.0
    L2 = l2/2.0
    PI_4 = np.pi/4.0
    PI_2 = np.pi/2.0

    
    l = lambda t: [ [L1,L2*np.tan(t)],
                    [L1/np.tan(t),L2] if np.tan(t) != 0 else [0,0],
                    [-L1,-L2*np.tan(t)],
                    [-L1/np.tan(t),-L2] if np.tan(t) != 0 else [0,0] ]
    
    return lambda t: l(t)[ np.int( (t+PI_4)/PI_2)%4 ]

                

def plotRect():
    '''
    Method to test the rectParameterization. This function
    should use the rectParameterization function to plot the
    corresponding rect.
    '''
    parameterDomain = np.arange(-np.pi/4.0,7*np.pi/4.0,0.01)
    rect = rectParameterization(4,8)
    
    curvePoints = np.array( [ rect(t) for t in parameterDomain ] )
    print(curvePoints)
    fig = plt.figure(1)
    plt.subplot(111)

    #curvePoints=np.array( [ [0,0],[1,1],[2,4],[3,9],[4,16] ] )
    plt.plot(curvePoints[:,0],curvePoints[:,1],"b-")
    plt.show()
    

def findAdequateK(D2,W,f):
    '''
    Consider a continuous version of the Heavside step function.
                           H_k(x) = 1/(1+e^(-kx)),
    for k>0.

    Furthermore, consider an alternative version:
                       H_k(x) = 2*(1/(1+e^(-kx)) - 0.5)

    Such that we have H_k(0) = 0. 

    An adequate k is defined as being a real number such that
                              H_k(W) < 1 - f,
    where 0 < f << 1.

    Notice that we have H_k(x) ~=1 for x>W.
    '''
    k=1
    while(1-UT.heaviside(k,W)>f):
        #print(1-UT.heaviside(k,1))
        k+=2
    return k

def findAdequateS(D2):
    '''
    Consider a Mirrowed version of the Heavside step function on
    the vertical plane x=0.

                       Ĥ_{k,s}(x) = 2*(1/(1+e^(k(x-s))-0.5),
    for k>0,s>0.

    The adequate s is simply the value of D2
    '''
    return D2

def testHeavisideFunctions(D2,W,k,s,tolerance):
    print( "H_k(W) > Ĥ_{k,s}(D2-W)",
           np.abs(UT.heaviside(k,W)-UT.mirrorHeaviside(k,D2-W,s)) < tolerance)

    print("H_k(W+tolerance) > Ĥ_{k,s}(D2-W+tolerance)",
          UT.heaviside(k,W+tolerance)> UT.mirrorHeaviside(k,D2-W+tolerance,s))
    print("H_k(W-tolerance) > Ĥ_{k,s}(D2-W-tolerance)",
          UT.heaviside(k,W-tolerance)< UT.mirrorHeaviside(k,D2-W-tolerance,s))

def plotHeavisideFunctions(D1,D2,k,s):
    domain = np.arange(D1,D2,0.01)
    imageH = [ UT.heaviside(k,x) for x in domain]
    imageMH = [ UT.mirrorHeaviside(k,x,s) for x in domain]
    
    fig = plt.figure(1)
    plt.subplot(211)
    plt.plot(domain,imageH,"b-")

    plt.subplot(212)
    plt.plot(domain,imageMH,"r-")
    plt.show()
 
    
def uniformMapping(beta):
    '''
    It maps the intervals:
              [beta-pi,beta) => [0,0.5)
              [beta,beta+pi] => [0.5,1]
    '''
    return lambda x: x/(2*np.pi) + 1.0/(2*np.pi)*(np.pi-beta)


def wMapping(D2,W):
    '''
    It maps the intervals:
             [0,0.5) => [0,w)
             [0.5,1] => [w,D2]

             A*e^(k*0.5) + C = W
             A*e^(k) + C = D2
    '''
    A = pow(W,2)/(D2-2*W)
    k = np.log(W/A+1)*2
    
    C = -A

    return lambda x: A*pow(np.e,k*x)+C


def wMappingReverse(D2,W):
    '''
    It maps the intervals:
             [0,0.5) => [0,D2-w)
             [0.5,1] => [D2-w,D2]

             A*e^(k*0.5) + = D2-W
             A*e^(k) - 1 = D2
    '''
    A = pow(D2-W,2)/(D2-2*(D2-W))
    k = np.log((D2-W)/A+1)*2
    C = -A

    print(A,k)

    return lambda x: A*pow(np.e,k*x)+C    

def plotUniformMapping(beta):
    '''

    Test uniform mapping
    '''
    f = uniformMapping(beta)

    domain = np.arange(beta-np.pi,beta+np.pi,0.01)
    image = [ f(x) for x in domain ]

    fig = plt.figure(1)
    plt.subplot(111)
    plt.axis('equal')
    plt.plot(domain,image,"b-")
    plt.show()
    

def plotWMapping(D2,W,beta):
    '''
    Tests wMapping
    '''
    um = uniformMapping(beta)
    wm = wMapping(D2,W)

    umDomain = np.arange(beta-np.pi,beta+np.pi,0.01)
    wmDomain = [um(x) for x in umDomain]
    
    image = [ wm(x) for x in wmDomain]

    

    fig = plt.figure(1)
    plt.subplot(111)
    plt.plot(umDomain,image,"b-")
    plt.show()
    

def plotWMappingReverse(D2,W,beta):
    '''
    Tests wMappingReverse
    '''
    um = uniformMapping(beta)
    wm = wMappingReverse(D2,W)

    umDomain = np.arange(beta-np.pi,beta+np.pi,0.01)
    wmDomain = [um(x) for x in umDomain]
    
    image = [ wm(x) for x in wmDomain]

    

    fig = plt.figure(1)
    plt.subplot(111)
    plt.plot(umDomain,image,"b-")
    plt.show()
    

def plotTest(D1,D2,W,f,tolerance,beta):
    '''
    For a given beta, plot the function:

                S(theta,beta) x G(theta,beta-pi/2) 
    on the interval 0 <= theta <= 2pi.
    '''
    k = findAdequateK(D2,W,f)
    s = findAdequateS(D2)

    print(k,f,s)
    
    domain = np.arange(beta-np.pi,beta+np.pi,0.01)

    wm = wMapping(D2,W)
    wrm = wMappingReverse(D2,W)
    um = uniformMapping(beta)

    domainWm = [ um(x) for x in domain ]
    
    G = lambda x: UT.heaviside(k,wm(um(x)))
    S = lambda x: UT.mirrorHeaviside(k,wrm(um(x)),s)

    imageG = [ G(x) for x in domain ]
    imageS = [ S(x) for x in domain ]
    imageGS = [ S(x)*G(x) for x in domain ]
    
    fig = plt.figure(1)
    plt.subplot(311)
    plt.plot(domain,imageG,"r-")

    plt.subplot(312)
    plt.plot(domain,imageS,"b-")

    plt.subplot(313)
    plt.plot(domain,imageGS,"g-")
    plt.show()


def main():
    D1 = 0
    D2 = 4
    W = 0.1
    f=1e-4
    tolerance=1e-8
    
    #plotRect()

    '''
    k=findAdequateK(D2,W,f)
    s=findAdequateS(D2)
    print(k)
    print(s)

    testHeavisideFunctions(D2,W,k,s,tolerance)
    plotHeavisideFunctions(D1,D2,k,s)
    '''
    
    #uniformMapping(2.0)
    #plotUniformMapping(1.0)

    #plotWMapping(D2,W,1.0)
    #plotWMappingReverse(D2,W,1.0)

    plotTest(D1,D2,W,f,tolerance,1.0)
    

if __name__=='__main__':
    main()
