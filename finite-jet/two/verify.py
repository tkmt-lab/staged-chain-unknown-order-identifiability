#!/usr/bin/env python3
"""Independent Fraction quotient-ring Taylor jets. No SymPy or generator imports."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json,argparse,time,resource,hashlib
class A:
    modulus=None
    def __init__(self,a=0):
        if isinstance(a,A): self.c=a.c;return
        c=[Q(x) for x in a] if isinstance(a,(list,tuple)) else [Q(a)]
        n=len(self.modulus)-1
        while len(c)>n:
            v=c.pop()/self.modulus[-1];j=len(c)-n
            for i in range(n):c[i+j]-=v*self.modulus[i]
        self.c=tuple(c+[Q(0)]*(n-len(c)))
    def __add__(self,b): b=A(b);return A([x+y for x,y in zip(self.c,b.c)])
    __radd__=__add__
    def __neg__(self):return A([-x for x in self.c])
    def __sub__(self,b):return self+-A(b)
    def __rsub__(self,b):return A(b)+-self
    def __mul__(self,b):
        b=A(b);n=len(self.c);v=[Q(0)]*(2*n-1)
        for i,x in enumerate(self.c):
            for j,y in enumerate(b.c):v[i+j]+=x*y
        return A(v)
    __rmul__=__mul__
    def __truediv__(self,b):
        b=A(b)
        if all(x==0 for x in b.c[1:]):return A([x/b.c[0] for x in self.c])
        # Multiplication by b as a rational matrix; exact Gauss elimination.
        n=len(self.c);columns=[(b*A([0]*j+[1])).c for j in range(n)]
        m=[[columns[j][i] for j in range(n)]+[self.c[i]] for i in range(n)]
        for j in range(n):
            pivot=next(i for i in range(j,n) if m[i][j]);m[j],m[pivot]=m[pivot],m[j]
            v=m[j][j];m[j]=[x/v for x in m[j]]
            for i in range(n):
                if i!=j:
                    v=m[i][j];m[i]=[x-v*y for x,y in zip(m[i],m[j])]
        return A([m[i][-1] for i in range(n)])
    def __eq__(self,b):return self.c==A(b).c
    def dump(self):return [str(x) for x in self.c]
def jet(m,left=False):
    beta=A(m['beta']);C=[A(1)] if left else list(map(A,m['C']))
    K=[[A(m['k'])]] if left else [list(map(A,row)) for row in m['K']]
    x=[A(m['X'])] if left else list(map(A,m['X']));n=len(x)
    states=[[A(m['S'])]]+[[v] for v in x]+[[A(m['R'])]];y=[]
    for j in range(10):
        y.append(beta*sum((states[0][i]*sum((C[a]*states[a+1][j-i] for a in range(n)),A()) for i in range(j+1)),A()))
        if j==9:break
        states[0].append(-y[j]/(j+1))
        for a in range(n):states[a+1].append((sum((K[a][b]*states[b+1][j] for b in range(n)),A())+(y[j] if a==0 else A()))/(j+1))
        states[-1].append(-sum((K[a][b]*states[b+1][j] for a in range(n) for b in range(n)),A())/(j+1))
    assert all(sum((u[j] for u in states),A())==0 for j in range(1,10))
    names=['S']+['X'+str(i+1) for i in range(n)]+['R']
    return {'states':{name:[(x*factorial(j)).dump() for j,x in enumerate(v)] for name,v in zip(names,states)},'output':[(x*factorial(j)).dump() for j,x in enumerate(y)]}
def validate(p):
    m=A.modulus
    assert len(m)==4
    lo,hi=map(Q,p['isolating_interval']);assert lo<hi
    def at(x):
        y=Q(0)
        for c in reversed(m):y=y*x+c
        return y
    assert at(lo)*at(hi)<0
    d,c,b,a=m
    discriminant=b*b*c*c-4*a*c**3-4*b**3*d-27*a*a*d*d+18*a*b*c*d
    assert discriminant<0  # exactly one real root for this real cubic
    prime=None
    for pr in (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97):
        if all(v.denominator%pr for v in m) and m[-1].numerator%pr:
            coeff=[v.numerator*pow(v.denominator,-1,pr)%pr for v in m]
            if all(sum(v*pow(x,j,pr) for j,v in enumerate(coeff))%pr for x in range(pr)):
                prime=pr;break
    assert prime is not None  # irreducible cubic reduction, hence irreducible over Q
    signs={}
    def pos(name,x):
        low=high=Q(0)
        for c0 in reversed(A(x).c):
            values=(low*lo,low*hi,high*lo,high*hi)
            low=min(values)+c0;high=max(values)+c0
        assert low>0,(name,str(low),str(high))
        signs[name]=[str(low),str(high)]
    left,right=p['left'],p['right']
    for name in ('S','X','R','beta'):pos('left '+name,left[name])
    pos('left exit',-A(left['k']))
    for name in ('S','R','beta'):pos('right '+name,right[name])
    for i in range(2):
        pos('right X'+str(i+1),right['X'][i]);pos('C'+str(i+1),right['C'][i])
    K=[[A(v) for v in row] for row in right['K']];C=list(map(A,right['C']))
    assert C[0]==1
    pos('K12',K[0][1]);pos('K21',K[1][0])
    for j in range(2):pos('negative column sum '+str(j+1),-K[0][j]-K[1][j])
    trace=K[0][0]+K[1][1];det=K[0][0]*K[1][1]-K[0][1]*K[1][0]
    pos('negative trace',-trace);pos('determinant',det)
    ck=[sum((C[i]*K[i][j] for i in range(2)),A()) for j in range(2)]
    assert K[1][0]!=0 and C[0]*ck[1]-C[1]*ck[0]!=0
    assert -trace==A(right['d1']) and det==A(right['d0']) and -trace+ck[0]==A(right['c'])
    return {'status':'PASS','cubic_irreducible_mod_prime':prime,'cubic_discriminant':str(discriminant),'exact_sign_intervals':signs}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('pair');ap.add_argument('--out',required=True);ap.add_argument('--reference');args=ap.parse_args()
    start=time.monotonic();raw=Path(args.pair).read_bytes();p=json.loads(raw);A.modulus=list(map(Q,p['minimal_polynomial']))
    validation=validate(p)
    left,right=jet(p['left'],True),jet(p['right']);diff=[A(b)-A(a) for a,b in zip(left['output'],right['output'])]
    assert all(x==0 for x in diff[:8]);assert diff[8]!=0
    for side in [left,right]:assert sum((A(v[0]) for v in side['states'].values()),A())==p['N0']
    if args.reference:
        ref=json.loads(Path(args.reference).read_text())
        assert ref['left']==left and ref['right']==right
    out={'status':'PASS','model_validation':validation,'field':'Q[k]/minimal_polynomial','method':'independent Fraction coefficient arithmetic and Taylor convolution','pair_sha256':hashlib.sha256(raw).hexdigest(),'left':left,'right':right,'difference_right_minus_left':[x.dump() for x in diff],'first_difference':8,'elapsed_seconds':time.monotonic()-start,'peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    Path(args.out).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({key:out[key] for key in ('status','first_difference','elapsed_seconds','peak_rss_KiB')}))
if __name__=='__main__':main()
