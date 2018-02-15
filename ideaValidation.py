#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import newton

import utils as UT


def testHeavisideFunctions(D2,W,f,tolerance):
    k = UT.findAdequateK(D2,W,f)
    s = UT.findAdequateS(D2)
    
    print( "H_k(W) > Ĥ_{k,s}(D2-W)",
           np.abs(UT.heaviside(k,W)-UT.mirrorHeaviside(k,D2-W,s)) < tolerance)

    print("H_k(W+tolerance) > Ĥ_{k,s}(D2-W+tolerance)",
          UT.heaviside(k,W+tolerance)> UT.mirrorHeaviside(k,D2-W+tolerance,s))
    print("H_k(W-tolerance) > Ĥ_{k,s}(D2-W-tolerance)",
          UT.heaviside(k,W-tolerance)< UT.mirrorHeaviside(k,D2-W-tolerance,s))
    

def plotRect():
    '''
    Method to test the rectParameterization. This function
    should use the rectParameterization function to plot the
    corresponding rect.
    '''
    parameterDomain = np.arange(-np.pi/4.0,7*np.pi/4.0,0.01)
    rect = UT.rectParameterization(4,8)
    
    curvePoints = np.array( [ rect(t) for t in parameterDomain ] )

    fig = plt.figure(1)
    plt.subplot(111)

    plt.plot(curvePoints[:,0],curvePoints[:,1],"b-")
    plt.show()
    

def plotHeavisideFunctions(D1,D2,W,f):
    k = UT.findAdequateK(D2,W,f)
    s = UT.findAdequateS(D2)
    
    domain = np.arange(D1,D2,0.01)
    imageH = [ UT.heaviside(k,x) for x in domain]
    imageMH = [ UT.mirrorHeaviside(k,x,s) for x in domain]
    
    fig = plt.figure(1)
    plt.subplot(211)
    plt.plot(domain,imageH,"b-")

    plt.subplot(212)
    plt.plot(domain,imageMH,"r-")
    plt.show()
 
    

def plotUniformMapping(beta):
    '''

    Test uniform mapping
    '''
    f = UT.uniformMapping(beta)

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
    um = UT.uniformMapping(beta)
    wm = UT.wMapping(D2,W)

    umDomain = np.arange(0,2*np.pi,0.01)
    wmDomain = [um(x) for x in umDomain]
    
    image = [ wm(x) for x in wmDomain]

    
    fig = plt.figure(1)
    plt.subplot(111)
    plt.plot(wmDomain,image,"b-")
    plt.show()
    

def plotWMappingReverse(D2,W,beta):
    '''
    Tests wMappingReverse
    '''
    um = UT.uniformMapping(beta)
    wm = UT.wMappingReverse(D2,W)

    umDomain = np.arange(0,2*np.pi,0.01)
    wmDomain = [um(x) for x in umDomain]
    
    image = [ wm(x) for x in wmDomain]

    

    fig = plt.figure(1)
    plt.subplot(111)
    plt.plot(wmDomain,image,"b-")
    plt.show()
    

def plotTest(D1,D2,W,f,tolerance,beta):
    '''
    For a given beta, plot the function:

                S(theta,beta) x G(theta,beta-pi/2) 
    on the interval 0 <= theta <= 2pi.
    '''
    k = UT.findAdequateK(D2,W,f)
    s = UT.findAdequateS(D2)
    
    domain = np.arange(0,2*np.pi,0.01)

    wm = UT.wMapping(D2,W)
    wrm = UT.wMappingReverse(D2,W)
    um = UT.uniformMapping(beta)

    domainWm = [ um(x) for x in domain ]
    domainHs = [ wm(um(x)) for x in domain ]
    
    
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
    W = 0.001
    f=1e-4
    tolerance=1e-8
    
    #plotRect()

    #testHeavisideFunctions(D2,W,f,tolerance)

    plotHeavisideFunctions(D1,D2,W,f)
    
    #plotUniformMapping(1.0)

    plotWMapping(D2,W,1.0)
    plotWMappingReverse(D2,W,1.0)

    plotTest(D1,D2,W,f,tolerance,1.0)
    

if __name__=='__main__':
    main()
