#!/usr/bin/env python3
"""Exact free-parameter n=5 checks of the all-n boundary argument.

The general index identity is proved in PROOF.md, not inferred from these
checks. Standard-library Fraction Laurent arithmetic; no Argon imports.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import argparse
import json
import resource
import sys
import time

NAMES=('x','u1','u2','u3','u4','u5','s','t','c2','c3','c4')
ZERO=(0,)*len(NAMES)


class Poly:
    def __init__(self,terms=None):
        self.terms={m:F(c) for m,c in (terms or {}).items() if c}

    @staticmethod
    def cast(v): return v if isinstance(v,Poly) else Poly({ZERO:v})

    def __add__(self,v):
        out=dict(self.terms)
        for m,c in self.cast(v).terms.items(): out[m]=out.get(m,F(0))+c
        return Poly(out)

    __radd__=__add__
    def __neg__(self): return Poly({m:-c for m,c in self.terms.items()})
    def __sub__(self,v): return self+-self.cast(v)
    def __rsub__(self,v): return self.cast(v)+-self

    def __mul__(self,v):
        out={}
        for m,c in self.terms.items():
            for n,d in self.cast(v).terms.items():
                k=tuple(a+b for a,b in zip(m,n))
                out[k]=out.get(k,F(0))+c*d
        return Poly(out)

    __rmul__=__mul__

    def __pow__(self,k):
        assert isinstance(k,int)
        if k<0:
            assert len(self.terms)==1
            m,c=next(iter(self.terms.items()))
            return Poly({tuple(k*i for i in m):c**k})
        out=Poly.cast(1)
        for _ in range(k): out=out*self
        return out

    def partial(self,i):
        out={}
        for m,c in self.terms.items():
            if m[i]:
                key=list(m); key[i]-=1
                out[tuple(key)]=c*m[i]
        return Poly(out)

    def coefficient(self,i,k):
        out={}
        for m,c in self.terms.items():
            if m[i]==k:
                key=list(m); key[i]=0
                out[tuple(key)]=c
        return Poly(out)


def var(i):
    m=list(ZERO); m[i]=1
    return Poly({tuple(m):1})


V={name:var(i) for i,name in enumerate(NAMES)}
x=V['x']; u=[V['u'+str(i)] for i in range(1,6)]
s,t,c2,c3,c4=(V[name] for name in ('s','t','c2','c3','c4'))


def derive(poly,velocities):
    out=Poly()
    for i,name in enumerate(NAMES):
        part=poly.partial(i)
        if part.terms:
            assert name in velocities,('missing derivative',name)
            out=out+part*velocities[name]
    return out


def invariants(ratios):
    tc=ratios[1]
    return [sum(((-tc)**j*ratios[i-j]*F(1,factorial(j)) for j in range(i+1)),Poly())
            for i in range(len(ratios))]


def run():
    checks=[]
    def eq(name,left,right):
        assert not (Poly.cast(left)-right).terms,name
        checks.append(name)

    J=sum(u); velocities={'x':-x*J,'u1':x*J-u[0]}
    velocities.update({'u'+str(i+1):u[i-1]-u[i] for i in range(1,5)})
    D=lambda f:derive(f,velocities)
    ratios=[Poly.cast(1)]+[v*u[0]**-1 for v in u[1:]]
    tc=ratios[1]; cs=invariants(ratios); ps=J*u[0]**-1
    eq('invariant c0=1',cs[0],1); eq('invariant c1=0',cs[1],0)
    forces={}
    for i in range(2,5):
        force=tc*cs[i-1]-cs[i]+(-tc)**i*F(1,factorial(i))
        forces[i]=force
        eq('moving invariant c'+str(i),D(cs[i]),x*ps*force)
        assert D(cs[i]).terms, 'must not freeze an invariant in the full system'
    checks.append('all three invariant lifts actually move under the full chain')

    relation=cs[4]-cs[2]**2
    eq('nonlinear relation chain rule',D(relation),x*ps*(forces[4]-2*cs[2]*forces[2]))
    abstract=[Poly.cast(1),Poly(),c2,c3,c4]
    G={i:t*abstract[i-1]-abstract[i]+(-t)**i*F(1,factorial(i)) for i in range(2,5)}
    idx=NAMES.index('t')
    for i in range(2,5):
        for k in range(2,5):
            eq(f'coefficient separation G{i}, t^{k}',G[i].coefficient(idx,k),F((-1)**i,factorial(i)) if i==k else 0)
    determinant=F(1,2)*F(-1,6)*F(1,24)
    assert determinant==F(-1,288)
    checks.append('partial-derivative coefficient matrix determinant=-1/288')
    eq('nonlinear relation has nonzero t^4 forcing', (G[4]-2*c2*G[2]).coefficient(idx,4),F(1,24))

    # Free invariant parameters describe the full x=0 space, a necessary control.
    vv=[sum((abstract[l]*t**(i-l)*F(1,factorial(i-l)) for l in range(i+1)),Poly()) for i in range(5)]
    states=[s*v for v in vv]
    delta=lambda f:derive(f,{'s':-s,'t':1,'c2':0,'c3':0,'c4':0})
    for i in range(5):
        eq('full Jordan boundary state '+str(i+1),delta(states[i]),(states[i-1] if i else 0)-states[i])
    recovered=invariants(vv)
    for i in range(2,5): eq('inverse invariant transform c'+str(i),recovered[i],abstract[i])
    coordinates=[u[0],tc,cs[2],cs[3],cs[4]]
    product=Poly.cast(1)
    for i,coordinate in enumerate(coordinates):
        for j in range(i+1,5): assert not coordinate.partial(NAMES.index('u'+str(j+1))).terms
        product=product*coordinate.partial(NAMES.index('u'+str(i+1)))
    eq('free-coordinate Jacobian determinant=1/u1^4',product,u[0]**-4)

    # Positions of the first nonzero boundary state, including the terminal axis.
    for k in range(2,6):
        sk=u[k-1]; h=u[0]*sk**-1; q=u[k-2]*sk**-1; pk=J*sk**-1
        eq(f'k={k} exact Ds',D(sk),sk*(q-1))
        eq(f'k={k} exact Dh',D(h),x*pk-h*q)
        if k<5:
            tail=[Poly.cast(1)]+[v*sk**-1 for v in u[k:]]
            tk=tail[1]; ck=invariants(tail)
            eq(f'k={k} exact Dt',D(tk),1-tk*q)
            for i in range(2,len(tail)):
                eq(f'k={k} moving tail invariant c{i}',D(ck[i]),q*(tk*ck[i-1]-ck[i]+(-tk)**i*F(1,factorial(i))))
    return checks,str(determinant)


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path)
    args=parser.parse_args(); resource.setrlimit(resource.RLIMIT_CPU,(60,60))
    start=time.perf_counter(); checks,det=run()
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result={'status':'PASS','check_groups':len(checks),'checks':checks,
            'n5_partial_coefficient_determinant':det,'arithmetic':'exact Fraction Laurent polynomials',
            'elapsed_seconds':time.perf_counter()-start,'peak_rss_bytes':rss if sys.platform=='darwin' else rss*1024,
            'cpu_limit_seconds':60,'rss_limit_enforced':False,'cost_usd':0,
            'limitations':'Free-parameter n=5 identities and controls only. The all-n generating identity, field/valuation argument and EL require their prose proofs; no independent/formal certification.'}
    encoded=json.dumps(result,indent=2)+'\n'
    if args.output: args.output.write_text(encoded)
    print(encoded,end='')
