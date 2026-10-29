#!/usr/bin/env python3
"""Exact n=4 algebra and Lyapunov controls; not a proof of EL or U4.

No Argon code, numerical integration, floating-point eigenvalues or packages.
"""
from fractions import Fraction as F
from itertools import permutations
from pathlib import Path
import argparse
import json
import resource
import sys
import time

NAMES=('x','u1','u2','u3','u4','a','lam')
ZERO=(0,)*len(NAMES)


class Poly:
    def __init__(self,terms=None):
        self.terms={m:F(c) for m,c in (terms or {}).items() if c}

    @staticmethod
    def cast(v):
        return v if isinstance(v,Poly) else Poly({ZERO:v})

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
        assert isinstance(k,int) and k>=0
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


def variable(i):
    m=list(ZERO); m[i]=1
    return Poly({tuple(m):1})


x,u1,u2,u3,u4,a,lam=(variable(i) for i in range(len(NAMES)))


def determinant(M):
    # Works for Fraction and Poly; a separate elimination solver is below.
    out=0
    n=len(M)
    for perm in permutations(range(n)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n))
        term=sign
        for i in range(n): term=term*M[i][perm[i]]
        out=out+term
    return out


def solve(A,b):
    n=len(b); aug=[list(map(F,row))+[F(v)] for row,v in zip(A,b)]
    for c in range(n):
        pivot=next(i for i in range(c,n) if aug[i][c])
        aug[c],aug[pivot]=aug[pivot],aug[c]
        divisor=aug[c][c]; aug[c]=[v/divisor for v in aug[c]]
        for i in range(n):
            if i!=c and aug[i][c]:
                scale=aug[i][c]
                aug[i]=[v-scale*w for v,w in zip(aug[i],aug[c])]
    return [row[-1] for row in aug]


def matrix_B(alpha):
    B=[[F(0) for _ in range(4)] for _ in range(4)]
    B[0]=[alpha-1,alpha,alpha,alpha]
    for i in range(1,4): B[i][i-1]=F(1); B[i][i]=F(-1)
    return B


def run():
    checks=[]
    def eq(name,left,right):
        assert not (Poly.cast(left)-right).terms,name
        checks.append(name)

    J=u1+u2+u3+u4; K=4*u1+3*u2+2*u3+u4; z=K+4*x
    velocity=[-x*J,x*J-u1,u1-u2,u2-u3,u3-u4,0,0]
    D=lambda f:sum((f.partial(i)*v for i,v in enumerate(velocity)),Poly())
    eq('population without removed class: D(x+J)=-u4',D(x+J),-u4)
    eq('weighted infection: DK=(4x-1)J',D(K),(4*x-1)*J)
    eq('Dz=-J',D(z),-J)
    eq('logarithmic conservation after clearing x: x Dz-Dx=0',x*D(z)-D(x),0)

    symbolic_B=matrix_B(a)
    char_matrix=[[(lam if i==j else 0)-symbolic_B[i][j] for j in range(4)] for i in range(4)]
    w=lam+1
    eq('n4 characteristic polynomial, free a and lambda',determinant(char_matrix),w**4-a*sum(w**k for k in range(4)))

    # Solve the 16-entry Lyapunov equation at a=1/8, strictly below 1/4.
    alpha=F(1,8); B=matrix_B(alpha); rows=[]; rhs=[]
    for i in range(4):
        for j in range(4):
            row=[F(0)]*16
            for k in range(4):
                row[4*k+j]+=B[k][i]
                row[4*i+k]+=B[k][j]
            rows.append(row); rhs.append(F(-int(i==j)))
    solution=solve(rows,rhs); P=[solution[4*i:4*i+4] for i in range(4)]
    assert all(P[i][j]==P[j][i] for i in range(4) for j in range(4))
    checks.append('rational Lyapunov solution is symmetric')
    for i in range(4):
        for j in range(4):
            val=sum(B[k][i]*P[k][j]+P[i][k]*B[k][j] for k in range(4))
            assert val==-int(i==j)
    checks.append('all 16 Lyapunov residual entries vanish')
    minors=[determinant([row[:k] for row in P[:k]]) for k in range(1,5)]
    assert all(v>0 for v in minors)
    checks.append('all four leading principal minors are strictly positive')

    # Equality at the threshold must not be covered by the stable argument.
    threshold=matrix_B(F(1,4))
    assert all(sum(row)==0 for row in threshold)
    checks.append('a=1/4 has zero eigenvector (1,1,1,1)')
    assert F(4)-1/F(1,4)==0
    checks.append('a=1/4 also makes H_x zero')
    assert F(4)-1/alpha!=0
    checks.append('a=1/8 has nonzero H_x')
    return checks,{'a':str(alpha),'P':[[str(v) for v in row] for row in P],
                   'leading_principal_minors':[str(v) for v in minors],
                   'threshold_control_a':'1/4','threshold_control':'zero infection eigenvalue and H_x=0'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(60,60))
    start=time.perf_counter(); checks,certificate=run()
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result={'status':'PASS','check_groups':len(checks),'checks':checks,'lyapunov_control':certificate,
            'arithmetic':'exact rational polynomial identities and linear algebra',
            'elapsed_seconds':time.perf_counter()-start,
            'peak_rss_bytes':rss if sys.platform=='darwin' else rss*1024,
            'cpu_limit_seconds':60,'rss_limit_enforced':False,'cost_usd':0,
            'limitations':'Identity and control checks only; no numerical orbit evidence or formal/independent certification of EL/U4.'}
    out=json.dumps(result,indent=2)+'\n'
    if args.output: args.output.write_text(out)
    print(out,end='')
