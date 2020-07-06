'''A test file of Numerical Analysis Lab 1.
Author: Linting
'''

from math import cos,pi
from matplotlib import pyplot
import interp
from linsolver import SUB, ABS
### Initial Settings
options = 1
runge = lambda x: 1/(1+x**2)
nodes = [x for x in range(-5,6)]
nodes2 = [5*cos((2*i + 1)/42*pi) for i in range(21)]
values = [runge(node) for node in nodes]
values2 = [runge(node) for node in nodes2]
derivates = [-2*x/(1+x**2)**2 for x in nodes]


X = [x/200 for x in range(-1000,1001)] #nodes to draw curves"

sys = interp.Interpolation(nodes, values, derivates)
sys2 = interp.Interpolation(nodes2,values2)
Y1 = sys.newton(X)
Y2 = sys2.lagrange(X)
Y3 = sys.linear(X)
Y4 = sys.hermite(X)
Y5 = sys.spline(X)

# output 
if options == 1:
    fig, ax = pyplot.subplots()
    line1 = ax.plot(X,Y1,X,[runge(x) for x in X])
    nodeplot = ax.plot(nodes, values, '.')
    ax.legend(["Lagrange Interpolation", "Exact Function"])
    pyplot.title("Question 1")

if options == 1:
    fig, ax = pyplot.subplots()
    line1 = ax.plot(X,Y2,X,[runge(x) for x in X])
    nodeplot = ax.plot(nodes2, values2, '.')
    ax.legend(["Newton Interpolation", "Exact Function"])
    pyplot.title("Question 2")

if options == 1:
    fig, ax = pyplot.subplots()
    line1 = ax.plot(X,Y3,X,[runge(x) for x in X])
    nodeplot = ax.plot(nodes, values, '.')
    ax.legend(["Piecewise Linear Interpolation", "Exact Function"])
    pyplot.title("Question 3")
if options == 1:
    fig, ax = pyplot.subplots()
    line1 = ax.plot(X,Y4,X,[runge(x) for x in X])
    nodeplot = ax.plot(nodes, values, '.')
    ax.legend(["Hermite Interpolation", "Exact Function"])
    pyplot.title("Question 4")
if options == 1:
    fig, ax = pyplot.subplots()
    line1 = ax.plot(X,Y5,X,[runge(x) for x in X])
    nodeplot = ax.plot(nodes, values, '.')
    ax.legend(["Cubic Spline Interpolation", "Exact Function"])
    pyplot.title("Question 5")

if options == 2:
    fig, ax = pyplot.subplots()
    error4 = SUB(Y4,[runge(x) for x in X])
    line1 = ax.plot(X,error4)
    pyplot.title("error for Hermite: max = %e, avg = %e" % (max(error4), sum(ABS(error4))/len(error4)))
if options == 2:
    fig, ax = pyplot.subplots()
    error5= SUB(Y5,[runge(x) for x in X])
    line1 = ax.plot(X,error5)
    pyplot.title("error for Spline: max = %e, avg = %e" % (max(error5), sum(ABS(error5))/len(error5)))
pyplot.show()