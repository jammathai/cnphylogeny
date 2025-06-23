#usage exmple: python DirectedCopyNumberDistance.py 12 21

def CalcPmin(ui,vi):
    return max(ui-vi,0)

def CalcPmax(ui):
    return max(ui-1,0)

def CalcM(ui,vi,ui1,vi1):
    return (vi-ui)-(vi1-ui1)

def limit(x,down,up):
    if x<down:
        return down
    if x>up:
        return up
    return x

class MyFunc:
    """This is my function"""

    def __init__(self,ui,vi,a,b,base):
        self.ui = ui
        self.vi = vi
        self.pmin = CalcPmin(ui,vi)
        self.pmax = CalcPmax(ui)
        self.a = a
        self.b = b
        self.base = base

    def CalcP(self,p):
        if p<self.pmin or p>self.pmax:
            raise("Error - unvalid p")
        if p<=self.a:
            return self.base
        if p<=self.b:
            return self.base + p - self.a
        if p<=self.pmax:
            return self.base - self.b - self.a + 2*p

##    def GetFirst(self,ui,vi):
##        pmin = CalcPmin(ui,vi)
##        pmax = CalcPmax(ui)
##        return MyFunc(ui,vi,pmin,pmin,vi-ui)

    def CalcNext(self,ui,vi,Qi):
        Mi = CalcM(ui,vi,self.ui,self.vi)
        pmin = CalcPmin(ui,vi)
        pmax = CalcPmax(ui)
        #print "Mi=",Mi
        #print "Qi=",Qi

        if Mi>=0:
            if Qi<=self.a:
                nextBase = self.base
                nextA = self.a-Mi
                nextB = self.b
            elif Qi<=self.b:
                nextBase = self.base + Qi - self.a
                nextA = Qi-Mi
                nextB = self.b
            elif Qi>=self.b:
                nextBase = self.base + Qi - self.a
                nextA = self.b-Mi
                nextB = Qi
        elif Mi<=0:
            if Qi<=self.a:
                nextBase = self.base
                nextA = self.a
                nextB = self.b-Mi
            elif Qi<=self.b:
                nextBase = self.base + Qi - self.a
                nextA = Qi
                nextB = self.b-Mi
            elif Qi>=self.b:
                nextBase = self.base + Qi - self.a
                nextA = min(self.b-Mi,Qi)
                nextB = max(Qi,self.b-Mi)

        #print "temp base=",nextBase
        #print "temp a=",nextA
        #print "temp b=",nextB
         
        if pmin > nextA and pmin<=nextB:
            nextBase = nextBase + pmin - nextA
        if pmin > nextB:
            nextBase = nextBase - nextB - nextA + 2*pmin

        nextA = limit(nextA,pmin,pmax)
        nextB = limit(nextB,nextA,pmax)

        return MyFunc(ui,vi,nextA,nextB,nextBase)
        
def DirectedCopyNumberDistanceLinear(u,v):
    n = len(u)
    N = max(max(u),max(v))

    u.append(N+1)
    v.append(N+1)

    u_m = 0
    prev = -1
    Q = {}
    for i in range(n):
        if v[i]==0:
            u_m = max(u_m,u[i])
        else:
            Q[i] = (u_m,prev)
            prev = i
            u_m = 0

    Q[n] = (u_m,prev)

    prevFunc = MyFunc(N+1,N+1,0,0,0)

    for i in range(n):
        if v[i]>0:
            prevFunc = prevFunc.CalcNext(u[i],v[i],Q[i][0])

    u_m,prev = Q[n]
    p_min = CalcPmin(u[prev],v[prev])
    d = prevFunc.CalcP(p_min)+max(u_m-p_min,0)
    for p in range(p_min+1,prevFunc.pmax+1):
        d = min(d,prevFunc.CalcP(p)+max(u_m-p,0))

    return d

# def stringToCNlist(s):
#     l=[]
#     for c in s:
#         l.append(int(float(c)))
#     return l

# def validateProfiles(u,v):
#     for i in range(len(u)):
#         if u[i]==0:
#             if v[i]==0:
#                 print "In i=",i," u[i]==v[i]==0 -> please delete this index"
#                 raise ValueError("Invalid input")
#             else:
#                 print "In i=",i," u[i]==0 but v[i]>0 -> cannot calculate distance"
#                 raise ValueError("Invalid input")

# import sys

# arguments = sys.argv
# # if len(arguments) != 3:
# #     print "Usage: X.py from to"
# #     raise ValueError("Usage: X.py from to")
# # if len(arguments[1])!=len(arguments[2]):
# #     print "Lengths do no match"
# #     raise ValueError("Lengths do no match")
# u =stringToCNlist(arguments[1])
# v=stringToCNlist(arguments[2])
# # validateProfiles(u,v)
# d = DirectedCopyNumberDistanceLinear(u,v)
# # print "The distance is ",d
