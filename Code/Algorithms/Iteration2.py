# #######################################################
#
# CESG 506 - Analysis of Nonlinear Structures
#
# file: Iteration1.py
#
# author: Peter Mackenzie-Helnwein
# created: 03/27/2020
#
# #######################################################

'''
This file is to demonstrate the iterative solution of
a nonlinear equation

Version 2: linear interpolation
'''

# import required libraries *****************************

import sys
import numpy as np

# import everything from my plotter class
from Plotter import *

# defining functions and classes ************************

def F(x):
    y = x*x - 10
    return y


def solver2(interval=None, TOL=1.0e-12):
    '''
    This solver uses linear interpolation to find the root of F(x)
    within  a given interval.  If no interval is given,
    the solver attempts to find a suitable interval.

    The solver uses a default tolerance of 10^(-12).
    You may overwrite the tolerance by specifying, e.g.,
    TOL=1.0e-5 when calling this solver.
    '''

    # finding an interval
    if not interval:
        xleft = 0
        yleft = F(xleft)
        xright = 1
        yright = F(xright)

        for i in range(10):
            if yleft*yright <=0:
                break
            xright *= 2
            yright = F(xright)

        if yleft*yright > 0.0:
            xright = 0
            yright = F(xright)
            xleft = -1
            yleft = F(xleft)

            for i in range(10):
                if yleft*yright <=0:
                    break
                xleft *= 2
                yleft = F(xleft)

        if yleft * yright > 0.0:
            print('cannot find suitable interval')
            return {}
        else:
            interval = [xleft,xright]


    # plotting the nonlinear function over that interval
    x = np.linspace(interval[0],interval[1],1000)
    y = F(x)

    # iterating for the root
    xi = []
    yi = []

    activeInterval = interval
    ys  = 1000000.  # anything larger than TOL will do
    cnt = 0         # set counter for emergency exit in case no solution is found

    while abs(ys) > TOL and cnt < 100:
        dx = activeInterval[1] - activeInterval[0]
        dy = F(activeInterval[1]) - F(activeInterval[0])

        s = activeInterval[0] - F(activeInterval[0]) / (dy/dx)
        ys = F(s)

        if F(activeInterval[0])*ys < 0.0:
            activeInterval[1] = s
        else:
            activeInterval[0] = s

        xi.append(s)
        yi.append(ys)

        cnt += 1

        print('{:3d}: x={:16.12f}  F(x)={:16.12e}'.format(cnt,s,ys))

    # collect iteration data
    answer = {
        'x':x,
        'y':y,
        'xi':xi,
        'yi':yi,
        'label':'Linear Interpolation',
        'save':'Iteration2.png'
    }

    # return the answer
    return answer


def main():
    ans = solver2()
    Plotter(**ans)


# main execution ****************************************

main()
sys.exit(0)
