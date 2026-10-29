#!/usr/bin/env python3
"""Own sparse integer-polynomial implementation; no producer code or CAS.

Full rank mod a prime proves full column rank over Q and C. In the one-column
kernel cases the predicted monomial is checked exactly over Z, so rank gives
the entire kernel. The cofactor parameter range comes from FIRST-STAGE.md,
not from this finite computation. Hard wall-clock cap: 60 seconds.
"""
import itertools
import json
from pathlib import Path
import resource
import time

PRIME = 1000003
START = time.monotonic()


def monomials(nvars, degree):
    return [e for e in itertools.product(range(degree+1), repeat=nvars) if sum(e)<=degree]


def add(poly, exponents, value):
    if value:
        poly[exponents] = poly.get(exponents,0)+value
        if not poly[exponents]: del poly[exponents]


def field(n):
    values = [{} for _ in range(n+1)]
    for j in range(1,n+1):
        e=[0]*(n+1); e[0]=1; e[j]=1
        add(values[0],tuple(e),-1); add(values[1],tuple(e),1)
    for j in range(1,n+1):
        e=[0]*(n+1); e[j]=1; add(values[j],tuple(e),-1)
        if j>1:
            e=[0]*(n+1); e[j-1]=1; add(values[j],tuple(e),1)
    return values


def derive_monomial(e, vectors):
    out={}
    for j,power in enumerate(e):
        if power:
            base=list(e); base[j]-=1
            for v,c in vectors[j].items():
                add(out,tuple(a+b for a,b in zip(base,v)),power*c)
    return out


def column(e,vectors,a,r,k):
    out=derive_monomial(e,vectors)
    add(out,e,r)
    v=list(e); v[0]+=1; add(out,tuple(v),-a)
    for j in range(1,len(e)):
        v=list(e); v[j]+=1; add(out,tuple(v),k)
    return out


def rank_mod(columns):
    basis={}
    for col in columns:
        if time.monotonic()-START>60: raise TimeoutError('60-second limit')
        v={e:c%PRIME for e,c in col.items() if c%PRIME}
        while v:
            pivot=max(v)
            if pivot not in basis:
                inv=pow(v[pivot],-1,PRIME)
                basis[pivot]={e:c*inv%PRIME for e,c in v.items()}
                break
            factor=v[pivot]
            for e,c in basis[pivot].items():
                new=(v.get(e,0)-factor*c)%PRIME
                if new: v[e]=new
                else: v.pop(e,None)
    return len(basis)


def run(n,d):
    mon=monomials(n+1,d); vectors=field(n); records=[]
    for a,r,k in itertools.product(range(d+1),repeat=3):
        cols=[column(e,vectors,a,r,k) for e in mon]
        known=None
        if n==1 and a==r and k+r<=d: known=(k,r)
        if n>=2 and a==r==0: known=(k,)+(0,)*n
        if known is not None:
            assert cols[mon.index(known)]=={}, 'exact known kernel failed'
        rank=rank_mod(cols); expected=len(mon)-(known is not None)
        assert rank==expected,(n,d,a,r,k,rank,expected)
        records.append({'a':a,'r':r,'k':k,'rank_mod_prime':rank,
                        'kernel_monomial_exponents':known})
    return {'n':n,'degree_bound':d,'columns':len(mon),'cofactor':'a*x-r-k*J',
            'cases':len(records),'results':records}


if __name__=='__main__':
    result={'arithmetic':'integer coefficient construction; rank over F_1000003',
            'prime':PRIME,'tests':[run(1,3),run(2,4),run(3,3)]}
    result.update(elapsed_seconds=time.monotonic()-START,
                  maxrss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  api_dollars=0,passed=True)
    target=Path(__file__).with_name('independent-results.json')
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='tests'},indent=2))
    print('cases:',sum(t['cases'] for t in result['tests']))
