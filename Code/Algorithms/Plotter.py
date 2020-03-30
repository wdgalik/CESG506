# #######################################################
#
# CESG 506 - Analysis of Nonlinear Structures
#
# file: Plotter.py
#
# author: Peter Mackenzie-Helnwein
# created: 03/27/2020
#
# #######################################################

'''
The Plotter class visualizes the progress of an iterative procedure

Usage:
    Plotter(ans)
    
    ans ... dictionary containing:
       x ....... an array of points on the x-axis
       y ....... an array of function values for each point in x
       xi ...... an array of x-values visited during the iteration
       yi ...... an array of y-values visited during the iteration
       label ... a string used as plot label
       saveTo .. a string defining the filename of the image file
'''

# import required libraries *****************************

import matplotlib.pyplot as plt
import numpy as np

# defining functions and classes ************************

# the Plotter class
class Plotter(object):
    '''
    class: Plotter

    variables:

    methods:

    '''

    def __init__(self, x, y, xi=[], yi=[], label='plotter image', save='', tangent=False, path=True):
        # remember input properties
        self.x = x
        self.y = y
        self.xi = xi
        self.yi = yi
        self.label = label
        self.filename = save
        self.plotPath = path
        self.plotTangent = tangent

        # create a figure

        # plot what we've got
        self.replot()

    def __str__(self):
        s = 'Plotter(\n'
        s += '   x:{},\n'.format(self.x)
        s += '   y:{},\n'.format(self.y)
        s += '   xi:{},\n'.format(self.xi)
        s += '   yi:{},\n'.format(self.yi)
        s += '   label:{},\n'.format(self.label)
        s += '   save:{}\n'.format(self.filename)
        s += '   )'
        return s

    def __repr__(self):
        s = 'Plotter('
        s += 'x:{}, '.format(repr(self.x))
        s += 'y:{}, '.format(repr(self.y))
        s += 'xi:{}, '.format(repr(self.xi))
        s += 'yi:{}, '.format(repr(self.yi))
        s += 'label:{}, '.format(repr(self.label))
        s += 'save:{}'.format(repr(self.filename))
        s += ')'
        return s

    def replot(self):
        plt.clf()
        fig, ax = plt.subplots()
        ax.grid(True)

        # draw zero line in black
        ax.plot(self.x,np.zeros_like(self.x),'k-')

        # draw the nonlinear function
        ax.plot(self.x,self.y,'b-')

        # create and plot vertical lines for iteration plot
        for line in [([s[0],s[0]],[0,s[1]]) for s in zip(self.xi, self.yi)]:
            ax.plot(*line,'r:')

        # create and plot tangents
        if self.plotTangent:
            lastPoint = []
            for point in zip(self.xi, self.yi):
                if lastPoint:
                    x = [lastPoint[0], point[0]]
                    y = [lastPoint[1], 0.0]
                    ax.plot(x,y,'r--')
                lastPoint = point

        # plot points visited during the iteration
        if self.plotPath:
            ax.plot(self.xi,self.yi,'ro:', fillstyle='none')
        else:
            ax.plot(self.xi,self.yi,'ro', fillstyle='none')

        if self.label:
            ax.set_title(self.label)

        if self.filename:
            fig.savefig(self.filename)

        plt.show()


# main function
def main():
    '''
    used as test and demo procedure
    if this file is executed directly:

    python3 Plotter.py
    '''
    x = np.linspace(0,2*np.pi,100)
    xi = np.linspace(0,2*np.pi,8)
    demo = {
        'x':x,
        'y':np.sin(x),
        'xi':xi,
        'yi':np.sin(xi),
        'label':'This is demo data',
        'save':'PlotterDemo.png'
    }
    Plotter(**demo)


# main execution ****************************************

if __name__ == '__main__':
    main()
