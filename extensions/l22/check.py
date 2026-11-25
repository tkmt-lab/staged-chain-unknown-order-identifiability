#!/usr/bin/env python3
"""Nonblind, independently written verifier. No research-side Python imports.
Exact QQ polynomials and outward-rounded dyadic intervals; no floating zero tests.
The input root is algebraic, NOT its rational center. See PROOF.md for jet induction.
"""
import argparse
import copy
from fractions import Fraction as Q
from math import factorial
import hashlib
import json
from pathlib import Path
import resource
import sys
import time
import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.domains import QQ

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'b22data'
BITS = 320
SCALE = 1 << BITS

def need(ok, msg):
    if not ok:
        raise ValueError(msg)

def frac(x):
    return Q(str(x))

class IV:
    """Endpoints are integer multiples of 2^-320. Every operation encloses."""
    def __init__(self, lo=0, hi=None):
        lo, hi = frac(lo), frac(lo if hi is None else hi)
        need(lo <= hi, 'reversed interval')
        self.l = (lo.numerator * SCALE) // lo.denominator
        self.h = -((-hi.numerator * SCALE) // hi.denominator)
    @classmethod
    def raw(cls, l, h):
        obj = object.__new__(cls); obj.l, obj.h = l, h
        return obj
    @staticmethod
    def coerce(v):
        return v if isinstance(v, IV) else IV(v)
    def __add__(self, other):
        other = self.coerce(other)
        return IV.raw(self.l + other.l, self.h + other.h)
    __radd__ = __add__
    def __neg__(self):
        return IV.raw(-self.h, -self.l)
    def __sub__(self, other): return self + (-self.coerce(other))
    def __rsub__(self, other): return self.coerce(other) + (-self)
    def __mul__(self, other):
        other = self.coerce(other)
        products = [a*b for a in (self.l,self.h) for b in (other.l,other.h)]
        return IV.raw(min(products)//SCALE, -((-max(products))//SCALE))
    __rmul__ = __mul__
    def reciprocal(self):
        need(self.l > 0 or self.h < 0, 'division interval contains zero')
        return IV(Q(SCALE,self.h), Q(SCALE,self.l))
    def __truediv__(self, other): return self * self.coerce(other).reciprocal()
    def __rtruediv__(self, other): return self.coerce(other) * self.reciprocal()
    def bounds(self): return [str(Q(self.l,SCALE)),str(Q(self.h,SCALE))]
    def mid(self): return Q(self.l+self.h, 2*SCALE)
    def absmax(self): return Q(max(abs(self.l),abs(self.h)),SCALE)
    def positive(self): return self.l > 0
    def negative(self): return self.h < 0
    def nonzero(self): return self.positive() or self.negative()

def det(A):
    """Bareiss determinant, exact division in QQ[a,b,c,u]."""
    A = [list(row) for row in A]
    n = len(A); old = A[0][0].ring.one; sign = 1
    for k in range(n-1):
        if not A[k][k]:
            j = next(j for j in range(k+1,n) if A[j][k])
            A[k],A[j] = A[j],A[k]; sign = -sign
        pivot = A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                A[i][j] = (pivot*A[i][j]-A[i][k]*A[k][j]).exquo(old)
        for i in range(k+1,n): A[i][k] = pivot.ring.zero
        old = pivot
    return sign*A[-1][-1]

def left_series(a,b,c,u,n=10):
    # Ordinary Taylor coefficients, not ordinary derivative values.
    one = a*0+1
    s = [one]; z = [one/2,u/2]; y = []
    for k in range(n+1):
        y.append(sum(s[i]*z[k-i] for i in range(k+1)))
        s.append(-y[k]/(k+1))
        if k<n:
            ynext = sum(s[i]*z[k+1-i] for i in range(k+2))
            z.append((-a*(k+1)*z[k+1]-b*z[k]+(k+1)*ynext+c*y[k])/((k+1)*(k+2)))
    q = list(s); q[0] += one/2
    r = []
    for k in range(n+1):
        r.append((y[k]-sum(q[i]*r[k-i] for i in range(1,k+1)))/q[0])
    return s,z,y,r

def build_polys(source):
    R,a,b,c,u = ring('a,b,c,u', QQ)
    s,z,y,r = left_series(a,b,c,u)
    need(r[0]==R.one/3 and r[1]==(u-R.one/6)/3, 'right initial jet mismatch')
    Y = [v*factorial(i) for i,v in enumerate(y)]
    U = [v*factorial(i) for i,v in enumerate(r)]
    rows = [[U[j+1],U[j],-Y[j+1],-Y[j],U[j+2]] for j in range(9)]
    saved = json.loads(source.read_text())
    parse = lambda expr: R.from_expr(sp.sympify(expr))
    need(len(saved['rows']) == len(rows), 'row count')
    need(all(len(old)==5 for old in saved['rows']), 'row width')
    need(all(x == parse(t) for row,old in zip(rows,saved['rows']) for x,t in zip(row,old)), 'row polynomial mismatch')
    D = det([row[:4] for row in rows[:4]])
    need(D == parse(saved['determinant']), 'pivot polynomial mismatch')
    P = []
    for j in range(4,9):
        p = det(rows[:4]+[rows[j]])
        need(p == parse(saved['conditions'][j-4]), f'P{j} mismatch')
        P.append(p)
        print(f'polynomial P{j}: exact match, {len(p)} terms', flush=True)
    return R,rows,D,P

def evaluate(p, box):
    powers = []
    for i,x in enumerate(box):
        arr = [IV(1)]
        for k in range(p.degree(i)): arr.append(arr[-1]*x)
        powers.append(arr)
    value = IV(0)
    for exps,coeff in p.items():
        term = IV(frac(coeff))
        for i,e in enumerate(exps): term *= powers[i][e]
        value += term
    return value

def linear_solve(A,b):
    A = [list(row)+[v] for row,v in zip(A,b)]
    n = len(b)
    for k in range(n):
        # No guessed interval pivot: division explicitly rejects zero inclusion.
        pivot = A[k][k]
        A[k] = [v/pivot for v in A[k]]
        for i in range(n):
            if i != k:
                v = A[i][k]
                A[i] = [x-v*y for x,y in zip(A[i],A[k])]
    return [row[-1] for row in A]

def certify_root(R,D,P,center,rho):
    box = [IV(x-rho,x+rho) for x in center]
    midbox = [IV(x) for x in center]
    derivatives = [[p.diff(v) for v in R.gens] for p in P[:4]]
    J0 = [[evaluate(p,midbox).mid() for p in row] for row in derivatives]
    inv = sp.Matrix(J0).inv()
    # Any nonsingular rational preconditioner suffices; use a 200-bit rounding.
    C = [[Q((frac(inv[i,j])*(1<<200)).numerator // (frac(inv[i,j])*(1<<200)).denominator, 1<<200) for j in range(4)] for i in range(4)]
    need(sp.Matrix(C).det()!=0, 'singular preconditioner')
    F0 = [evaluate(p,midbox) for p in P[:4]]
    JX = [[evaluate(p,box) for p in row] for row in derivatives]
    A = [[int(i==j)-sum(C[i][k]*JX[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
    norms = [sum(x.absmax() for x in row) for row in A]
    offsets = [-sum(C[i][k]*F0[k] for k in range(4))+sum(A[i][j]*IV(-rho,rho) for j in range(4)) for i in range(4)]
    need(all(v<1 for v in norms), 'not a contraction')
    need(all(x.absmax()<rho for x in offsets), 'root image not strictly internal')
    d = evaluate(D,box); need(d.nonzero(), 'singular contact matrix')
    gap = -Q(3,2)*evaluate(P[4],box)/d
    need(gap.negative(), 'tenth derivative gap not negative')
    return box, {'contraction_row_bounds':list(map(str,norms)), 'offsets':[x.bounds() for x in offsets], 'preconditioner':[[str(v) for v in row] for row in C], 'pivot':d.bounds(),'gap10':gap.bounds()}, gap

def realize(a,b,c,beta,s0,z0,u,d,w,pop):
    one = a*0+1
    k = [[d-a,w*(d*(a-d)-b)/(c-d)],[(c-d)/w,-d]]
    alpha = c-a; B = (d*(a-c)-b)/(c-d); v = u-beta*s0; J = z0/beta
    x = [J*(B-v)/(B-alpha), J*(v-alpha)/(w*(B-alpha))]
    return {'K':k,'C':[one,w],'beta':beta,'S':s0,'X':x,'R':pop-s0-sum(x)}

def model_checks(model):
    k,c,x = model['K'],model['C'],model['X']
    ck = [sum(c[i]*k[i][j] for i in range(2)) for j in range(2)]
    quantities = {'S':model['S'],'X1':x[0],'X2':x[1],'R':model['R'],'beta':model['beta'],
        'C1':c[0],'C2':c[1], 'K12':k[0][1],'K21':k[1][0],
        'negative_trace':-k[0][0]-k[1][1], 'det_K':k[0][0]*k[1][1]-k[0][1]*k[1][0],
        'negative_colsum1':-k[0][0]-k[1][0], 'negative_colsum2':-k[0][1]-k[1][1],
        'controllability_det':k[1][0], 'observability_det':c[0]*ck[1]-c[1]*ck[0]}
    for name,val in quantities.items():
        need(val.nonzero() if name=='observability_det' else val.positive(), f'model fails {name}')
    need(c[0].l==c[0].h==SCALE, 'Ce1 not exactly one')
    return {name:val.bounds() for name,val in quantities.items()}

def symbolic_model_checks():
    a,b,c,beta,s,z,u,d,w,N = sp.symbols('a b c beta s z u d w N')
    model = realize(a,b,c,beta,s,z,u,d,w,N)
    K,C,X = sp.Matrix(model['K']),sp.Matrix([model['C']]),sp.Matrix(model['X'])
    exprs = [K.trace()+a,K.det()-b,(C*K)[0]+a-c,(C*X)[0]-z/beta,
             (C*K*X)[0]-(u-beta*s)*z/beta,model['R']+s+sum(X)-N,
             C[0]-1]
    need(all(sp.cancel(e)==0 for e in exprs), 'symbolic realization mismatch')
    # Cayley-Hamilton supplies the reduced second-order ODE including beta.
    need(all(sp.cancel(e)==0 for e in K*K+a*K+b*sp.eye(2)), 'Cayley Hamilton mismatch')
    return {'exact_rational_identities':len(exprs)+4,'Ce1':1}

def symbolic_transport_checks():
    # Universal difference recurrence, QQ coefficients: eta=T-q, delta=z*-r.
    # d=y*-y=q*delta+r*eta+eta*delta, eta'=-d,
    # delta''=-A*delta'-B*delta+g*d'+v*d-E.
    q=sp.symbols('q0:11'); r=sp.symbols('r0:11')
    A,B,g,v,e8=sp.symbols('A B g v e8')
    eta=[sp.S.Zero]; delta=[sp.S.Zero,sp.S.Zero]
    def difference(n):
        return sp.expand(sum(q[i]*delta[n-i]+r[i]*eta[n-i]+eta[i]*delta[n-i] for i in range(n+1)))
    out=[]
    for n in range(9):
        dn=difference(n); out.append(dn)
        eta.append(-dn/(n+1))
        dn1=difference(n+1)
        forcing=e8 if n==8 else sp.S.Zero
        delta.append(sp.expand((-A*(n+1)*delta[n+1]-B*delta[n]+g*(n+1)*dn1+v*dn-forcing)/((n+1)*(n+2))))
    out.append(difference(9)); eta.append(-out[9]/10)
    out.append(difference(10))
    need(out[:10]==[0]*10, 'contact transport zero orders')
    need(sp.expand(factorial(10)*out[10]+q[0]*factorial(8)*e8)==0, 'contact transport leading coefficient')
    return {'zero_orders':list(range(10)), 'ordinary_derivative_gap':'-q0 * E^(8)(0)', 'exact_symbolic_recurrence':True}

def arithmetic_checks():
    # Exact rational sample containment, including negative divisors and zeros.
    samples=[(Q(-7,3),Q(-1,5)),(Q(-1,3),Q(2,7)),(Q(0),Q(5,11)),(Q(1,7),Q(9,4))]
    count=0
    for a,b in samples:
        x=IV(a,b)
        for c,d in samples:
            y=IV(c,d)
            for op,exact in [(lambda:x+y,lambda p,q:p+q),(lambda:x-y,lambda p,q:p-q),(lambda:x*y,lambda p,q:p*q)]+([] if c<=0<=d else [(lambda:x/y,lambda p,q:p/q)]):
                enclosure=op()
                for p in [a,(a+b)/2,b]:
                    for q in [c,(c+d)/2,d]:
                        val=exact(p,q)
                        need(Q(enclosure.l,SCALE)<=val<=Q(enclosure.h,SCALE), 'rational interval containment test')
                        count+=1
    try:
        IV(1)/IV(-1,1)
    except ValueError:
        pass
    else:
        raise ValueError('zero denominator not rejected')
    return {'rational_containment_samples':count,'zero_divisor_rejected':True}

def full_ode_jets(model,n=10):
    k,c,x0 = model['K'],model['C'],model['X']
    s=[model['S']]; x=[[v] for v in x0]; r=[model['R']]; y=[]
    for j in range(n+1):
        yy=model['beta']*sum(s[i]*sum(c[h]*x[h][j-i] for h in range(2)) for i in range(j+1))
        y.append(yy)
        s.append(-yy/(j+1))
        kx=[sum(k[h][l]*x[l][j] for l in range(2)) for h in range(2)]
        r.append(-sum(kx)/(j+1))
        for h in range(2): x[h].append((kx[h]+(yy if h==0 else 0))/(j+1))
    return [v*factorial(j) for j,v in enumerate(y)]

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path,default=HERE/'result.json');args=parser.parse_args()
    start=time.monotonic()
    pair_path=SOURCE/'lower_PAIR10.json'; system_path=SOURCE/'lower_system_u.json'
    pair=json.loads(pair_path.read_text())
    need(pair['fixed']=={'S0':'1','beta':'1','z0':'1/2','T0':'3/2'},'fixed slice mismatch')
    center=[frac(pair['center'][i]) for i in pair['active']]; rho=frac(pair['radius'])
    need(pair['active']==[0,1,2,4] and len(center)==4,'variable order')
    need(frac(pair['center'][3])==Q(1,2) and frac(pair['center'][5])==Q(3,2),'fixed center entries')
    R,rows,D,P=build_polys(system_path)
    print('certifying root with polynomial Jacobian',flush=True)
    box,root,gap=certify_root(R,D,P,center,rho)
    need(frac(pair['output_gap_coarse'][0])<Q(gap.l,SCALE) and Q(gap.h,SCALE)<frac(pair['output_gap_coarse'][1]), 'claimed gap enclosure fails')
    vals=[[evaluate(v,box) for v in row] for row in rows[:4]]
    ast,bst,g,v=linear_solve([row[:4] for row in vals],[-row[4] for row in vals])
    need(g.positive(),'right beta')
    a,b,c,u=box
    pop=IV(pair['population'])
    dl,dr=[frac(p['d']) for p in pair['realization']]; wl,wr=[frac(p['weight']) for p in pair['realization']]
    left=realize(a,b,c,IV(1),IV(1),IV(Q(1,2)),u,IV(dl),IV(wl),pop)
    right=realize(ast,bst,v/g,g,IV(Q(3,2)),IV(Q(1,3)),u-Q(1,6),IV(dr),IV(wr),pop)
    identities=symbolic_model_checks()
    transport=symbolic_transport_checks()
    arithmetic=arithmetic_checks()
    checks=[model_checks(m) for m in (left,right)]
    print('positive minimal models verified',flush=True)
    # Direct ODE jets: separation can reject; zero containment NEVER accepts equality.
    yl,yr=full_ode_jets(left),full_ode_jets(right)
    direct=[(b-a).bounds() for a,b in zip(yl,yr)]
    need((yr[10]-yl[10]).negative(),'direct ODE tenth derivative separation fails')
    controls=[]
    for kind in ['coefficient','initial_state']:
        bad=copy.deepcopy(right)
        if kind=='coefficient': bad['K'][0][0]-=Q(1,1000); order=1
        else:
            bad['X'][0]+=Q(1,1000); bad['R']-=Q(1,1000); order=0
        badchecks=model_checks(bad)
        badjets=full_ode_jets(bad)
        delta=badjets[order]-yl[order]
        need(delta.nonzero(),f'{kind} mutation not rejected')
        controls.append({'kind':kind,'mutation':'K11 -= 1/1000' if kind=='coefficient' else 'X1 += 1/1000; R -= 1/1000 (population balance)',
                         'status':'REJECT','first_separated_order':order,'difference':delta.bounds(),'still_positive_minimal':True,'model_checks':badchecks})
    result={'status':'PASS','claim_L':'成立','independence':'single-agent nonblind; own checker; no research Python executed/imported',
        'root':root,'rho':str(rho),'root_box':[x.bounds() for x in box],
        'polynomial_comparison':{'rows':len(rows),'row_entries':sum(map(len,rows)), 'determinants':1+len(P),'all_exact':True},
        'symbolic_model_checks':identities,'symbolic_transport':transport,'arithmetic_checks':arithmetic,'model_checks':checks,
        'exact_contact':{'zero_orders':list(range(10)),'basis':'certified P4=P5=P6=P7=0, nonsingular M, Schur complement and Taylor induction in PROOF.md; not interval containment'},
        'direct_ode_difference_enclosures_diagnostic_only':direct,'negative_controls':controls,
        'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (pair_path,system_path)},
        'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'python':sys.version,'sympy':sp.__version__,'dyadic_bits':BITS,
        'elapsed_seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print('PASS: (L), two positive minimal models, 0..9 exact by root/induction, order 10 separated; two controls REJECT',flush=True)

if __name__=='__main__':
    main()
