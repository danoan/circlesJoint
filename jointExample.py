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
        

def main():
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



if __name__=='__main__':
    main()
