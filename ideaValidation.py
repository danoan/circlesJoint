#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

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
                              H_k(1) < 1 - f,
    where 0 < f << 1.

    Notice that we have H_k(x) ~=1 for x>1.
    '''
    k=1
    while(1-UT.heaviside(k,1)>f):
        #print(1-UT.heaviside(k,1))
        k+=2
    return k

def findAdequateS(D2,k,W,tolerance):
    '''
    Consider a Mirrowed version of the Heavside step function on
    the vertical plane x=0.

                       Ĥ_{k,s}(x) = 1/(1+e^(k(x-s))),
    for k>0,s>0.

    Define
                              ff = 1 - H_k(1).

    An adequate s is defined such that

                              1-Ĥ_{k,s}(1)=ff.

    Notice that the following statements hold
                        1)  H_k(1) == Ĥ_{k,s}(1).
                        2)  H_k(x) > H_{k,s}(x), for x > 1.
                        3)  H_k(x) < H_{k,s}(x), for x < 1.
    '''
    ff = 1-UT.heaviside(k,1)
    s = -1.0/k * np.log(ff/(1-ff)) + (D2-W)

    h = 1-UT.heaviside(k,1)
    mh = 1-UT.mirrorHeaviside(k,1,s)
    
    step = 0.1
    while( np.abs(h-mh) > tolerance ):
        s-=step
        mh = 1-UT.mirrorHeaviside(k,1,s)
    
    return s

def testHeavisideFunctions(k,s,tolerance):
    print(UT.heaviside(k,1))
    print(UT.mirrorHeaviside(k,1,s))
    print( np.abs(UT.heaviside(k,1)-UT.mirrorHeaviside(k,1,s)) < tolerance)

    print(UT.heaviside(k,1+tolerance)> UT.mirrorHeaviside(k,1+tolerance,s))
    print(UT.heaviside(k,1-tolerance)< UT.mirrorHeaviside(k,1-tolerance,s))

def plotHeavisideFunctions(k,s):
    domain = np.arange(-2,2,0.01)
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
    return lambda x: x/(2*np.pi) + beta-np.pi


def wMapping(D2,w,beta):
    '''
    It maps the intervals:
             [0,0.5) => [0,w)
             [0.5,1] => [w,D2]
    '''
    pass


def wMappingReverse(D2,w,beta):
    '''
    It maps the intervals:
             [0,0.5) => [0,D2-w)
             [0.5,1] => [D2-w,D2]
    '''
    pass

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
    

def plotWMapping():
    '''
    Tests wMapping
    '''
    pass

def plotWMappingReverse():
    '''
    Tests wMappingReverse
    '''
    pass

def plotTest(beta):
    '''
    For a given beta, plot the function:

                S(theta,beta) x G(theta,beta-pi/2) 
    on the interval 0 <= theta <= 2pi.
    '''
    pass


def main():
    D1 = 0
    D2 = 4
    W = 1
    f=1e-4
    tolerance=1e-10
    
    #plotRect()

    '''
    k=findAdequateK(D2,W,f)
    s=findAdequateS(D2,k,W,tolerance)
    print(k)
    print(s)

    testHeavisideFunctions(k,s,tolerance)
    plotHeavisideFunctions(k,s)
    '''
    
    uniformMapping(2.0)
    plotUniformMapping(1.0)

if __name__=='__main__':
    main()
