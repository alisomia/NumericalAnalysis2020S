'''
This file contains interpolation class. 
For further information, please use `help(interp.Interpolation)`
Implemented by Linting.
Latest Updated: 20200320
'''
import math
from linsolver import *
from functools import reduce


eps = 1e-6  # machine precision declared 


class Interpolation:
    ''' The Interpolation class implements the classical interpolation methods.
    It requires a nodes set, a value set and a derivatives set (optional)
    Use methods `lagrange`, `newton`, `linear`, `hermite`, `spline` to compute the interpolation value .
    '''
    def __init__(self, nodes, values, derivates = []):
        self.nodes = nodes
        self.values = values
        self.derivatives = derivates # if Hermite is adopted.
        
        if len(nodes)<2:
            raise ValueError("Need more than one point!")
    def lagrange(self, X):
        total = INIT(X)
        for i in range(len(self.nodes)):
            node = self.nodes[i]
            local = reduce(MUL, (SUBf(X,f) for f in self.nodes if abs(f-node)>eps))
            d = reduce(lambda x,y:x*y, (node-f for f in self.nodes if abs(f-node)>eps))
            total = ADD(total,MULf(local, self.values[i]/d))
        return total

    def newton(self, X):
        l = self.values
        diff_quoient = [l[0]]
        for k in range(1,len(self.nodes)):
            l = DIV(SUB(l[1:],l[:-1]), SUB(self.nodes[k:],self.nodes[:-k]))
            diff_quoient.append(l[0])
        return reduce(lambda x,op:ADDf(MUL(x,op[0]), op[1]), zip([SUBf(X,f) for f in self.nodes[::-1]],diff_quoient[::-1]), INIT(X))

    def linear(self, X): # assume that both nodes and X are ordered
        cnt = 0
        total = INIT(X)
        for i in range(len(X)):
            while X[i] > self.nodes[cnt+1]:
                cnt = cnt + 1
            h = self.nodes[cnt+1]-self.nodes[cnt]
            x = (X[i] - self.nodes[cnt])/h
            y = 1-x
            total[i] = self.values[cnt+1] * x +self.values[cnt]*y
        return total

    def hermite(self, X):
        cnt = 0
        total = INIT(X)
        for i in range(len(X)):
            while X[i] > self.nodes[cnt+1] :
                cnt = cnt + 1
            # do somthing here
            h = self.nodes[cnt+1]-self.nodes[cnt]
            x = (X[i] - self.nodes[cnt])/h
            y = 1-x
            total[i] = self.values[cnt] * (1 + 2*x)*y*y +\
                self.values[cnt+1] * (1 + 2*y)*x*x+\
                    self.derivatives[cnt]*x*h*y*y - \
                        self.derivatives[cnt+1]*y*h*x*x
        return total
    def spline(self,X):
        x = self.nodes
        y = self.values
        n = len(x)-1
        dx = SUB(x[1:],x[:-1])
        Lam = [1,*DIV(dx[:-1], ADD(dx[1:], dx[:-1])),0]
        Mu = INIT(x)
        for i in range(1,n):
            Mu[i] = 3 * (1-Lam[i])/dx[i-1]*(y[i]-y[i-1]) + 3*Lam[i]/dx[i]*(y[i+1] - y[i])
        Mu[0] = 3*(y[1]-y[0])/dx[0]
        Mu[n] = 3*(y[n] - y[n-1])/dx[n-1]
        L = ADDf(MULf(Lam[1:],-1),1)
        D = [2 for _ in x]
        U = ADDf(Lam[:-1],0)
        Mat = TrilinearMatrix(D,U,L)
        
        if isinstance(self.derivatives, list) and len(self.derivatives)>0:
            print("The derivates value has been overwrite!")
        self.derivatives = Mat.solve(Mu)
        
        return self.hermite(X)

        
