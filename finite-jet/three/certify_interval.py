#!/usr/bin/env python3
"""Rigorous dyadic interval arithmetic: every endpoint operation uses integers.
Intervals are outward rounded to multiples of 2**-BITS. No floating point
enters an enclosure. Dual arithmetic encloses the full Jacobian on a box.
"""
from fractions import Fraction
import json
BITS=384
Q=1<<BITS
def ceildiv(a,b): return -((-a)//b)
class I:
    def __init__(self,x=0):
        if isinstance(x,I): self.a,self.b=x.a,x.b; return
        f=Fraction(str(x)) if not isinstance(x,Fraction) else x
        self.a=f.numerator*Q//f.denominator
        self.b=ceildiv(f.numerator*Q,f.denominator)
    @classmethod
    def raw(cls,a,b):
        r=object.__new__(cls);r.a=a;r.b=b;assert a<=b;return r
    def __add__(self,o):
        o=I(o);return I.raw(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return I.raw(-self.b,-self.a)
    def __sub__(self,o):return self+-I(o)
    def __rsub__(self,o):return I(o)+-self
    def __mul__(self,o):
        o=I(o);v=[self.a*o.a,self.a*o.b,self.b*o.a,self.b*o.b]
        return I.raw(min(v)//Q,ceildiv(max(v),Q))
    __rmul__=__mul__
    def __truediv__(self,o):
        o=I(o);assert not o.a<=0<=o.b, "division interval contains zero"
        fs=[Fraction(a*Q,b) for a in (self.a,self.b) for b in (o.a,o.b)]
        lo,hi=min(fs),max(fs)
        return I.raw(lo.numerator//lo.denominator,ceildiv(hi.numerator,hi.denominator))
    def __rtruediv__(self,o):return I(o)/self
    def __pow__(self,n):
        assert n>=0
        r=I(1)
        for _ in range(n):r=r*self
        return r
    def mid(self):return I.raw((self.a+self.b)//2,(self.a+self.b)//2)
    def absbound(self):return max(abs(self.a),abs(self.b))
    def bounds(self):return [str(Fraction(self.a,Q)),str(Fraction(self.b,Q))]
    def approx(self):return [float(Fraction(self.a,Q)),float(Fraction(self.b,Q))]
class Dual:
    def __init__(self,v,d):self.v=I(v);self.d=list(map(I,d))
    def coerce(self,o):return o if isinstance(o,Dual) else Dual(o,[0]*len(self.d))
    def __add__(self,o):
        o=self.coerce(o);return Dual(self.v+o.v,[a+b for a,b in zip(self.d,o.d)])
    __radd__=__add__
    def __neg__(self):return Dual(-self.v,[-a for a in self.d])
    def __sub__(self,o):return self+-self.coerce(o)
    def __rsub__(self,o):return self.coerce(o)+-self
    def __mul__(self,o):
        o=self.coerce(o);return Dual(self.v*o.v,[a*o.v+self.v*b for a,b in zip(self.d,o.d)])
    __rmul__=__mul__
    def __truediv__(self,o):
        o=self.coerce(o);return Dual(self.v/o.v,[(a*o.v-self.v*b)/(o.v*o.v) for a,b in zip(self.d,o.d)])
    def __rtruediv__(self,o):return self.coerce(o)/self
    def __pow__(self,n):
        r=self.coerce(1)
        for _ in range(n):r=r*self
        return r

def equations(z):
    # Fixed common time coordinate tau=t/T; D_tau=T*D_t, T=1/25.
    i,g,p2,p1,p0,nb,nq1,nq0=z
    T=Fraction(1,25); d=Fraction(-59,200)
    zero=z[0]*0;S=[zero]*12;X=[zero]*12;Y=[zero]*12;V=[zero]*12
    S[0]=zero+1;X[0]=i
    for k in range(12):
        Y[k]=sum(S[j]*X[k-j] for j in range(k+1))
        V[k]=(Y[k]-sum(S[j]*V[k-j] for j in range(1,k+1)))/(1+d)
        if k<11:
            S[k+1]=-T*Y[k]/(k+1)
            X[k+1]=T*(Y[k]-g*X[k])/(k+1)
    from math import factorial
    cols=[(V,2),(V,1),(V,0),(Y,2),(Y,1),(Y,0)]
    coeff=[p2,p1,p0,nb,nq1,nq0]
    return [sum(c*Fraction(factorial(k+j),factorial(k))*seq[k+j]
                for c,(seq,j) in zip(coeff,cols))+
            Fraction(factorial(k+3),factorial(k))*V[k+3] for k in range(9)]

def inv_fraction(a):
    n=len(a);m=[list(row)+[Fraction(int(k==j)) for j in range(n)] for k,row in enumerate(a)]
    for k in range(n):
        pivot=max(range(k,n),key=lambda j:abs(m[j][k]));m[k],m[pivot]=m[pivot],m[k]
        v=m[k][k];assert v;m[k]=[x/v for x in m[k]]
        for j in range(n):
            if j!=k:
                v=m[j][k];m[j]=[x-v*y for x,y in zip(m[j],m[k])]
    return [row[n:] for row in m]

def certify(candidate_path,out_path):
    data=json.load(open(candidate_path))
    vals=data["root_center"]
    c=[I(x).mid() for x in vals];n=len(c);assert n==8
    r=I("1e-35").b
    box=[I.raw(x.a-r,x.b+r) for x in c]
    dual=[Dual(box[j],[int(j==k) for k in range(n)]) for j in range(n)]
    res=equations(dual);J=[f.d for f in res[:8]]
    Jmid=[[Fraction(a.mid().a,Q) for a in row] for row in J]
    inverse=inv_fraction(Jmid)
    Y=[[I(a).mid() for a in row] for row in inverse] # fixed dyadic matrix
    fc=equations(c)
    E=[[I(int(j==k))-sum(Y[j][l]*J[l][k] for l in range(n)) for k in range(n)] for j in range(n)]
    rad=max(sum(v.absbound() for v in row) for row in E)
    correction=[-sum(Y[j][k]*fc[k] for k in range(n)) for j in range(n)]
    dx=[I.raw(-r,r) for _ in c]
    image=[c[j]+correction[j]+sum(E[j][k]*dx[k] for k in range(n)) for j in range(n)]
    inside=all(b.a<x.a and x.b<b.b for b,x in zip(box,image))
    contract=rad<Q
    gap=res[8].v
    assert inside and contract
    assert gap.b<0 or gap.a>0
    result={"method":"exact integer dyadic interval Krawczyk/contraction","bits":BITS,
      "candidate":candidate_path,"radius":"1e-35 (rounded upward to dyadic)",
      "strict_inclusion":inside,"jacobian_contraction":contract,
      "infinity_norm_bound_exact":str(Fraction(rad,Q)),
      "infinity_norm_bound_approx":float(Fraction(rad,Q)),
      "root_box_exact":[b.bounds() for b in box],"image_exact":[x.bounds() for x in image],
      "residual8_exact":gap.bounds(),"residual8_approx":gap.approx(),
      "scope":"Existence and uniqueness of root8 in box, residual8 nonzero. Physical constraints separately certified."}
    open(out_path,"w").write(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k not in ("root_box_exact","image_exact","residual8_exact")},indent=2))
    return box,result

if __name__=="__main__":
    import sys
    certify(sys.argv[1],sys.argv[2])
