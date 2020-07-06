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
    K = [0 for i in range(M)]
    
    while t<T-1e-8:
        if t+h>T+1e-8:
            h = T-t
        for i in range(M):
            K[i] = f(t+a[i]*h, x + h*sum([b[i][j]*K[j] for j in range(i)]))
        x = x + h*sum([c[i]*K[i] for i in range(M)])
        xlist.append(x)
        t = t+h
    return xlist,x

def RK4(x0, f, T, h):
    a = [0,1/2,1/2,1]
    b = [[],[1/2],[0,1/2],[0,0,1]]
    c = [1/6,1/3,1/3,1/6]
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
        xs = xlist[-1] + h * sum([c_pred[i]*f(tlist[-(i+1)], xlist[-(i+1)]) for i in range(N)])
        x = x + h* (c_corr[0]*f(t+h, xs) + sum([c_corr[i+1]*f(tlist[-(i+1)], xlist[-(i+1)]) for i in range(N-1)]))
        t = t+h
        xlist.append(x)
        tlist.append(t)
    return xlist,x

def Adams4(x0,f,T,h):
    c_pred = [55/24, -59/24, 37/24, -9/24]
    c_corr = [9/24, 19/24, -5/24, 1/24]
    return PredCorr(x0, f, T, h, c_pred, c_corr)

# _, x = Adams4(1, lambda t,x:x, 1, 0.1)
# print("%e" % (2.718281828459045235360 - x) )
# _, x = RK4(1, lambda t,x:x, 1, 0.1)
# print("%e" % (2.718281828459045235360 - x) )

from matplotlib import pyplot
