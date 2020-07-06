#!/usr/bin/env python3
'''This file contains basic operations of linear algebra without importing numpy.
It is first used in Numerical analysis course.
Implemented by Linting.
Latest Updated: 20200320
'''
import math


## Some element-wise operations between two vectors
ADD = lambda  A,B : A if B==0 else [x+y for x,y in zip(A,B)] 
SUB = lambda A,B :[x-y for x,y in zip(A,B)]
MUL = lambda  A,B : INIT(A) if B==0 else [x*y for x,y in zip(A,B)]
DIV = lambda A,B: [x/y for x,y in zip(A,B)]

## Some num-vec operations
ADDf = lambda A,f:[x+f for x in A]
SUBf = lambda A,f: [x-f for x in A]
MULf = lambda A,f: 0 if A==0 else [x*f for x in A]
#DIVf = lambda A,f: [x/f for x in A]

## Zero Initiation with the same size as the given vector
INIT = lambda A: [0 for x in A]
ABS = lambda A: [abs(x) for x in A]

def SUM(l):
    if len(l)==0:
        return 0
    res = l[0]
    for x in l[1:]:
        res = ADD(res,x)
    return res

class TrilinearMatrix:
    '''
    This class implements a trilinear matrix with both its representation and LU factorization.
    Once an instance is constructed, use `solve` directly to obtain the result.
    The program will automatically compute LU fact when first calling `solve`.
    '''
    def __init__(self, _d, _u, _l):
        '''initialization
         _d： diagnoal (size n)
         _l: lower diagnoal (size n-1)
         _u: upper diagnoal (size n-1)
         '''
        self.d, self.u, self.l = _d, _u, _l
        self.L, self.D, self.U = [],[],[]
    def factor(self):
        "compute LU factor"
        L = INIT(self.l)
        U = self.u
        D = INIT(self.d)
        D[0] = self.d[0]
        L[0] = self.l[0]/D[0]
        for i in range(len(L)):
            D[i+1] = self.d[i+1] -L[i]*U[i]
            if i<len(L)-1:
                L[i+1] = self.l[i+1]/D[i+1]
        self.L, self.D, self.U = L, D, U
    def solve(self, y):
        "solve"
        if len(self.D)==0:
            self.factor()
        L,D,U = self.L, self.D, self.U
        # solve Lu = y
        for i in range(len(L)):
            y[i+1] -= y[i]*L[i]
        # solve Ux = u
        for i in range(len(L),0,-1):
            y[i] = y[i]/D[i]
            y[i-1] -= y[i]*U[i-1]
        y[0] = y[0]/D[0]
        return y


if __name__ == '__main__':
    print('This is a module file, demo is not established yet.')
