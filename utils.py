import numpy as np
import matplotlib.pyplot as plt

def heaviside(k,x):
    return (2.0/(1.0+pow(np.e,-k*x)) - 1.0 )

def mirrorHeaviside(k,x,s):
    return (2.0/(1.0+pow(np.e,k*(x-s) ) ) - 1.0)

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
    while(1-heaviside(k,W)>f):
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

def uniformMapping(beta):
    '''
    Input: -i*pi0 <= beta <= i*pi for i > 0

    It maps the intervals:
              [beta-pi,beta) => [0,0.5)
              [beta,beta+pi] => [0.5,1]

    Remark: In order to map the lambda input interval
    [0,2pi] to the mapped interval [beta-pi,beta], we 
    apply the transformation:

                 x* = x + (beta-pi)
    '''
    
    return lambda x: ( (x - beta + np.pi + 8*np.pi)%(2*np.pi) )/(2*np.pi)

def wMapping(D2,W):
    '''
    Input: 0 < W < D2
    
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

    return lambda x: A*pow(np.e,k*x)+C    
