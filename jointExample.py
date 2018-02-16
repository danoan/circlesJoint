#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
'''

The goal is to construct a closed curve from arc segments of 8 circles.
I will further compute the squared curvature functional for the resulting curve
and compare with a closed curve of same length build with 4 arc segments and 4
straight lines.

The ultimate goal is to show that the functional at the first curve has lower
value than the second.

'''

def circlePoint(cx,cy,radius,theta):
    return (cx+radius*np.cos(theta),cy+radius*np.sin(theta))

def arcSegmentPoints(cx,cy,radius,startTheta,endTheta):
    '''
         Counterclockwise evaluation
    '''

    domain = np.arange(startTheta,endTheta,0.01)

    x = [ circlePoint(cx,cy,radius,x)[0] for x in domain ]
    y = [ circlePoint(cx,cy,radius,x)[1] for x in domain ]

    return (x,y)

def uniformJoint(n,a):
    k = 4*np.pi/n

    circles = [ [4,0,k+a,0,0] ]

    for i in range(1,n):
        pc = circles[-1]

        r = 4 if i%2==0 else 32
        rs = r + pc[0]
        if i%2==1: #Concave
            tto = np.pi + pc[2]
            tfrom = tto - a

            cx = pc[3] + rs*np.cos(pc[2])
            cy = pc[4] + rs*np.sin(pc[2])
            
        else: #Convex
            tfrom = np.pi + pc[1]
            tto = tfrom + k + a

            cx = pc[3] + rs*np.cos(pc[1])
            cy = pc[4] + rs*np.sin(pc[1])

        circles.append( [r,tfrom,tto,cx,cy] )

    return circles


def testUniformJoint():
    circles = uniformJoint(8,np.pi/12)

    fig = plt.figure(1)
    plt.subplot(111)
    i=0
    for c in circles:
        x,y = arcSegmentPoints(c[3],c[4],c[0],c[1],c[2])
        print(c)

        if i%2==0:
            plt.plot(x,y,'b-')
        else:
            plt.plot(x,y,'r-')

        i+=1

    plt.show()

def getAngleIntersectionLine(x,y):
    '''
    Je veux trouver la longuer l du segment de la 
    ligne L(theta) qui part de l'origine et va jusqu'a 
    l'intersection de L(theta) avec le point (x,y).

         tan(theta) = y/x
         theta = arctan(y/x)
         
         l = x/cos(theta)     
    '''
    theta =( (np.arctan2(y,x))+2*np.pi )%(2*np.pi)
    l = np.abs( x/np.cos(theta) if np.cos(theta)!=0 else y/np.sin(theta) )
    
    return (l,theta)

def getInterval(x,y,a,initiale):
    '''
    x et y sont les coordenees d'un point terminale d'un quart
    de circle de rayon a. La flag indique se se trait d'une terminale
    initiale ou finale.
    '''
    (l,theta) = getAngleIntersectionLine(x,y)
    angleDiff = 2*np.pi*a/(4*l)
    print(l,theta,angleDiff)
    if initiale:
        alpha = theta
        beta = theta+angleDiff
    else:
        alpha = theta-angleDiff
        beta = theta

    return (alpha,beta,l)


class QuartCircle:
    def __init__(self,cx,cy,rayon,interval,thetaAdjust):
        self.cx = cx
        self.cy = cy
        self.rayon = rayon
        
        self.interval = interval
        self.alpha,self.beta,self.ad = interval

        if self.alpha > 2*np.pi or self.alpha < 0:
            raise(BaseException("Alpha is not in the interval [0,2pi]: %.4f" % (self.alpha,)))
        if self.beta > 2*np.pi or self.beta < 0:
            raise(BaseException("Beta is not in the interval [0,2pi]: %.4f" % (self.beta,)))
            
        
        self.thetaAdjust = thetaAdjust
        
    def evaluate(self,theta):
        if theta < self.alpha or theta > self.beta:
            return (0,0)
            #raise("Attempt of evaluation of a quart of circle out of its interval")

            
        return circlePoint(self.cx,self.cy,self.rayon,self.thetaAdjust(theta))

class LineSegment:
    def __init__(self,interval,x,y):
        self.interval = interval
        self.alpha,self.beta = interval
        self.x = x
        self.y = y

        if self.alpha > 2*np.pi or self.alpha < 0:
            raise(BaseException("Alpha is not in the interval [0,2pi]: %.4f" % (self.alpha,)))
        if self.beta > 2*np.pi or self.beta < 0:
            raise(BaseException("Beta is not in the interval [0,2pi]: %.4f" % (self.beta,)))
        

    def evaluate(self,theta):
        if theta < self.alpha or theta > self.beta:
            return (0,0)
            #raise("Attempt of evaluation of a line segment out of its interval")        

        print(self.alpha,self.beta,theta)
        return ( self.x(theta), self.y(theta) )
    
def evaluateJointRoundedCurve(geometricObjects,theta):
    p = [0,0]
    for g in geometricObjects:
        s = g.evaluate(theta)
        p[0]+=s[0]
        p[1]+=s[1]

    return p

def roundedCorners(m,a):
    '''
    On cosidere le carre de cote m.

    On veut derive une curve C semblant a un carre
    avec de carrefours arrondit. Ça veut dire qu'elle
    est composee de 4 lignes de longueur m* et de 4
    quart de circles q1,q2,q3,q4 de rayon a, un par 
    chaque carrefour (on considere le sense counterclock,
    et que q1 represente la quart de circle pour le 
    carrefour le plus a gauche et le plus haut)
    
    q1--------q4
    |         |
    |         |
    |         |
    q2-------q3

    La carre il y a une longueur de 4m

    La nouvelle curve doit avoir la même longueur, donc
                  2pi*a + 4m* = 4m
    donc,
                  m* = (4m-2pi*a)/4

    Il faut trouver la centre de chaque quart de circle.

    En considerant une evaluation de les quart de circles
    dans le sense counterclockwise, on represente chaque 
    quart de circle qi par une pair de points (ui,vi), ou
                  ui: Point initiale
                  vi: Point finale

    Donc, on impose:
        1) v1[x] = -m   ;   u1[y] = m
        2) u2[x] = -m   ;   v2[y] = -m
        3) v3[x] = m    ;   u3[y] = -m
        4) u4[x] = m    ;   v4[y] = m

    Pour trouver le center (cx,cy) du quart de circle q1

         cx + a*cos(pi) = -m
         cy + a*sin(pi/2) = m

    Pour q2

         cx + a*cos(pi) = -m
         cy + a*sin(3pi/2) = -m

    Pour q3

         cx + a*cos(0) = m
         cy + a*sin(3pi/2) = -m

    Pour q4

         cx + a*cos(0) = m
         cy + a*sin(pi/2) = m
    '''
    cx1 = -m - a*np.cos(np.pi)
    cy1 = m - a*np.sin(np.pi/2)

    cx2 = -m - a*np.cos(np.pi)
    cy2 = -m - a*np.sin(3*np.pi/2)

    cx3 = m - a*np.cos(0)
    cy3 = -m - a*np.sin(3*np.pi/2)

    cx4 = m - a*np.cos(0)
    cy4 = m - a*np.sin(np.pi/2)

    '''
    Ensuite, je veux donner une parameterization de 
    R(theta) pour le carre et une parameterization
    C(theta) pour la nouvelle curve.

               |
       ________|________
       |       |  /     |
      P1_      | /     P4_
       |       |/       |
   ----|-------x--------|----
       |       |        |
      P2_      |       P3_
       |_______|________|
               |
               |
    On considere une ligne l(theta) en partant de l'origine.

    Ces deux parameterizations sont tels que
             R(theta) = intersection l(theta) avec R
             C(theta) = intersection l(theta) avec C


    La parameterization de R est deja code sur la fonction 
    UT.rectParameterization(theta)

    Pour la nouvelle curve C, je dois identifier les theta 
    lequelles j'utilise la parameterization d'un quart de
    circle et les thetas que je doit utiliser la parameterization
    d'une recte.

    Pour faire ça, chaque qi nous donne une intervalle 
    Ii = [alpha, beta] lequel on doit plotter la quart de curve 
    qi si theta belong to Ii.

    '''
    P1 = (-m,cy1)
    P2 = (-m,cy2)
    P3 = (m,cy3)
    P4 = (m,cy4)

    I1 = (alpha1,beta1,len1) = getInterval(P1[0],P1[1],a,False)
    I2 = (alpha2,beta2,len2) = getInterval(P2[0],P2[1],a,True)
    I3 = (alpha3,beta3,len3) = getInterval(P3[0],P3[1],a,False)
    I4 = (alpha4,beta4,len4) = getInterval(P4[0],P4[1],a,True)

    print(I1,I2,I3,I4)

    q1 = QuartCircle(cx1,cy1,a,I1, lambda t: np.pi/2 + (t-alpha1)*len1/a )
    q2 = QuartCircle(cx2,cy2,a,I2, lambda t: np.pi + (t-alpha2)*len2/a)
    q3 = QuartCircle(cx3,cy3,a,I3, lambda t: 3*np.pi/2 + (t-alpha3)*len3/a)
    q4 = QuartCircle(cx4,cy4,a,I4, lambda t: 0 + (t-alpha4)*len4/a)

    iz = zip( [I1,I2,I3,I4],[I2,I3,I4,I1] )
    LI=[ (i1[1],i2[0]) for i1,i2 in iz ]

    l1 = LineSegment(LI[0],lambda t: -m, lambda t: -m*np.tan(t))
    l2 = LineSegment(LI[1],lambda t: -m/np.tan(t), lambda t: -m)
    l3 = LineSegment(LI[2],lambda t: m, lambda t: m*np.tan(t))
    l4 = LineSegment(LI[3],lambda t: m/np.tan(t), lambda t: m)

    
    GO = [q1,q2,q3,q4,l1,l2,l3,l4]

    return lambda t: evaluateJointRoundedCurve(GO,t)
    
    
    
def main():
    m = 8
    a = 4
    frc = roundedCorners(m,a)

    domainP = np.arange(0,2*np.pi,0.01)

    domain = [ frc(t)[0] for t in domainP ]
    image = [ frc(t)[1] for t in domainP ]

    fig = plt.figure()
    plt.subplot(111)
    plt.plot(domain,image,'b-')

    print(image)

    plt.show()


if __name__=='__main__':
    main()
