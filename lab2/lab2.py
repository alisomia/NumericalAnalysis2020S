from math import exp,sqrt,log2,pi
#f2 = lambda x:x**3/(exp(x)-1) if x>0 else 0 # the original function
f2 = lambda x:x**3/(exp(x)-1)/exp(x) if x>0 else 0 # the function after processing

def simple_int(l,r,w,c,func):
    ''' the simple integration
    l,r : the interval
    w,c : weight array
    func : function
    '''
    h = r-l
    return h*sum([wi*func(l+ci*h) for (wi,ci) in zip(w,c)])

def _getquadrature(method):
    ''' a function supports the arrays `w` and `c`
    the method should be : 'midpoint','trapezoid','simpson','cotes3','cotes4','gauss2','gauss3'
    '''
    if method == "midpoint":
        return ([1],[.5])
    elif method == "trapezoid":
        return ([.5,.5],[0,1])
    elif method == "simpson":
        return ([1/6,4/6,1/6],[0,1/2,1])
    elif method == "cotes3":
        return ([1/8,3/8,3/8,1/8],[0,1/3,2/3,1])
    elif method == "cotes4":
        return ([7/90,32/90,12/90,32/90,7/90],[0,1/4,2/4,3/4,1])
    elif method == "cotes5":
        return (  [19/288, 25/96 , 25/144,25/144,25/96 ,19/288],[0,1/5,2/5,3/5,4/5,1])
    elif method == "cotes6":
        return ([  41/840	 ,  9/35 	,   9/280	 , 34/105	  , 9/280	,   9/35 	 , 41/840],\
            [0,1/6,2/6,3/6,4/6,5/6,1])
    elif method == "gauss2":
        return ([1/2,1/2],[.5-.5/sqrt(3),.5+.5/sqrt(3)])
    elif method == "gauss3":
        return ([5/18,8/18,5/18],[.5-.5*sqrt(.6),.5,.5+.5*sqrt(.6)])
    elif method == "gauss4":
        return ([0.347854845,	0.652145155,	0.652145155,	0.347854845], \
            list(map(lambda x: (x+1)/2, [-0.861136312,-0.339981044,	0.339981044,	0.861136312])))
    elif method == "gauss5":
        return ([0.236926885,	0.47862867,	0.568888889,	0.47862867,	0.236926885],\
             list(map(lambda x:(x+1)/2,[-0.906179846,	-0.53846931,	0	,0.53846931	,0.906179846])))
    else:
        return method
        
def composite_int(l,r,func, nNode = 1, method = "simpson"):
    '''calculate the composite integration.
    l,r: interval
    func: function
    nNode: number of sub-interval
    method: the method used in each sub-interval
    '''
    (w,c) = _getquadrature(method)
    h = (r-l)/nNode
    return sum([simple_int(l+i*h,l+(i+1)*h,w,c,func) for i in range(nNode)])

def adaptive_int(l,r,func, tol = 1e-8, method = "simpson"):
    '''calculate the adaptive integration.
    l,r: interval
    func: function
    tol: tolerance
    method: the method used in each sub-interval

    return a tuple (result, #node)
    '''
    (w,c) = _getquadrature(method)
    m = (l+r)/2
    if abs(simple_int(l,m,w,c,func) + simple_int(m,r,w,c,func) -simple_int(l,r,w,c,func)) < tol * (r-l):
        return (simple_int(l,r,w,c,func),1)
    elif (r-l) < tol:
        return (simple_int(l,r,w,c,func),1)
    else:
        infol = adaptive_int(l,m,func,tol,method)
        infor = adaptive_int(m,r,func,tol,method)
        return (infol[0] + infor[0], infol[1] + infor[1])  # add the total node

if __name__ == "__main__":
    print("Composite Numerical Integral with #interval = 50, 100, 200")
    exact = pi**4/15 - 6 #adaptive_int(0,30,f2, 1e-16,"cotes4")[0]
    print("                      ", 50, "                      ", 100, "                           ", 200,"                           ", 400)
    for method in ["midpoint","trapezoid","simpson","cotes3","cotes4","gauss2","gauss3"]:
        err50 = abs(exact - composite_int(0,30,f2,50,method))
        err100 = abs(exact - composite_int(0,30,f2,100,method))
        err200 = abs(exact - composite_int(0,30,f2,200,method))
        err400 = abs(exact - composite_int(0,30,f2,400,method))
        print(method, "error=",err50, err100, "(",int(log2(err50/err100)*100)/100,")", err200, "(",int(log2(err100/err200)*100)/100,")"\
            ,err400, "(",int(log2(err200/err400)*100)/100,")")

    print("\nAdaptive Procedure with tolerence 1e-16:")
    for method in ["simpson","cotes3","cotes4","gauss2","gauss3"]:
        i2 = adaptive_int(0,30,f2, 1e-16, method)
        print(method,"nNode", i2[1], "value=",i2[0])
    
    print("Romberg Technique for 4th order scheme")
    print("                      ", 50, "                      ", 100, "                           ", 200)
    for method in ["midpoint","trapezoid","simpson","cotes3","gauss2"]:
        err50 = abs((16*composite_int(0,30,f2,100,method) - composite_int(0,30,f2,50,method))/15 - exact)
        err100 = abs((16*composite_int(0,30,f2,200,method) - composite_int(0,30,f2,100,method))/15 - exact)
        err200 = abs((16*composite_int(0,30,f2,400,method) - composite_int(0,30,f2,200,method))/15 - exact)
        print(method, "error=",err50, err100, "(",int(log2(err50/err100)*100)/100,")", err200, "(",int(log2(err100/err200)*100)/100,")")
    print("Romberg Technique for 6th order scheme")
    print("                      ", 50, "                      ", 100, "                           ", 200)
    for method in ["cotes4","gauss3"]:
        err50 = abs((64*composite_int(0,30,f2,100,method) - composite_int(0,30,f2,50,method))/63 - exact)
        err100 = abs((64*composite_int(0,30,f2,200,method) - composite_int(0,30,f2,100,method))/63 - exact)
        err200 = abs((64*composite_int(0,30,f2,400,method) - composite_int(0,30,f2,200,method))/63 - exact)
        print(method, "error=",err50, err100, "(",int(log2(err50/err100)*100)/100,")", err200, "(",int(log2(err100/err200)*100)/100,")")

