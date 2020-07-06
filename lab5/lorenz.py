from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
from ode import RK4
def lorenz(x, sigma, rho, beta):
    return [sigma*(x[1]-x[0]), rho*x[0]-x[1]-x[0]*x[2], x[0]*x[1]-beta*x[2]]

def plot(xlist):
    X = [x[0] for x in xlist]
    Y = [x[1] for x in xlist]
    Z = [x[2] for x in xlist]
    fig = pyplot.figure()
    ax1 = Axes3D(fig)
    ax1.plot3D(X,Y,Z,'gray')
x0 = [-10,10,25]
#x0 = [0,0,1]
#x0 = [8.5,8.5,27]
#x0 = [1,0,1]
#x0 = [0.01, 0, 1]
#x0 = [1,20,3]
T = 100
h = 0.01
sigma = 10
rho = 28
beta = 100
xlist,_ = RK4(x0,lambda t,x: lorenz(x, sigma, rho, beta), T, h)
#print(xlist)
plot(xlist)
pyplot.show()