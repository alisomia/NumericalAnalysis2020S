from linsolver import ADD, MULf, SUM, INIT
from math import log10
def RK(x0, f, T, h, a, b,c):
    ''' x0 the initial value
        F(t,x) the RHS
        T the total time
        h the time step
        a,b,c  the butcher tableau
        '''
    t = 0
    x = x0
    xlist = []
    xlist.append(x0)
    M = len(c)
    K = [INIT(x) for i in range(M)]
    
    while t<T-1e-8:
        if t+h>T+1e-8:
            h = T-t
        for i in range(M):
            K[i] = f(t+a[i]*h, ADD(x, MULf(SUM([ MULf(K[j], b[i][j]) for j in range(i)]), h)))
        x = ADD(x, MULf(SUM([MULf(K[i],c[i]) for i in range(M)]), h))
        xlist.append(x)
        t = t+h
    return xlist,x

def RK4(x0, f, T, h):
    a = [0,1/2,1/2,1]
    b = [[],[1/2],[0,1/2],[0,0,1]]
    c = [1/6,1/3,1/3,1/6]
    return RK(x0, f, T, h, a,b,c)

def Euler(x0, f, T, h):
    a = [1]
    b = [[]]
    c = [1]
    return RK(x0, f, T, h, a,b,c)

def iEuler(x0, f, T, h):
    a = [0,1]
    b = [[],[1/2]]
    c = [1/2,1/2]
    return RK(x0, f, T, h, a,b,c)

def Heun(x0, f, T,h):
    a = [0,2/3]
    b = [[], [2/3]]
    c = [1/4, 3/4]
    return RK(x0, f, T, h, a,b,c)

def Kutta3(x0, f, T, h):
    a = [0,1/2,1]
    b = [[],[1/2],[-1,2]]
    c = [1/6,2/3,1/6]
    return RK(x0, f, T, h, a,b,c)

def PredCorr(x0, f, T, h, c_pred, c_corr):
    tlist = [0,h, 2*h, 3*h]
    xlist, x = RK4(x0, f, 3*h, h)
    N = len(c_pred)
    t = tlist[-1]
    while t<T-1e-8:
        if t+h>T+1e-8:
            h = T-t
        #print(t,[tlist[-(i+1)] for i in range(N)])
        xs = ADD(xlist[-1],  MULf(SUM([MULf(f(tlist[-(i+1)], xlist[-(i+1)]), c_pred[i]) for i in range(N)]), h))
        x = ADD(x , MULf( ADD(MULf(f(t+h, xs),c_corr[0]) , SUM([MULf(f(tlist[-(i+1)], xlist[-(i+1)]), c_corr[i+1])for i in range(N-1)])), h))
        t = t+h
        xlist.append(x)
        tlist.append(t)
    return xlist,x



def Adams4(x0,f,T,h):
    c_pred = [55/24, -59/24, 37/24, -9/24]
    c_corr = [9/24, 19/24, -5/24, 1/24]
    return PredCorr(x0, f, T, h, c_pred, c_corr)



if __name__ == "__main__":
    # do accuracy test
    f = lambda t,x:x
    x0 = [1]
    e = 2.718281828459045235360
    for method in [Euler, iEuler, Heun, Kutta3, RK4, Adams4]:
        err2 = abs(e-method(x0, f, 1, 0.01)[1][0])
        err3 = abs(e-method(x0, f, 1, 0.001)[1][0])
        print("%e %e %.2f" % (err2, err3, log10(err2/err3)))
# _, x = Adams4([1,1], lambda t,x:[x[0], 2*x[1]], 1, 0.0001)
# print("%e" % (2.718281828459045235360 - x[0]) )
# print("%e" % (2.718281828459045235360**2 - x[1]) )
# _, x = RK4(1, lambda t,x:x, 1, 0.1)
# print("%e" % (2.718281828459045235360 - x) )

