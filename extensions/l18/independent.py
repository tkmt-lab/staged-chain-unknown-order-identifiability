#!/usr/bin/env python3
"""Independent stage-1 verifier. Only research inputs: INPUT.json, CANDIDATE.json.
Uses mpmath's outward-rounded interval context; no research-side code imported.
Run from any directory: python3 -B independent.py > stage1-results.json
"""
import hashlib, json, math, platform, resource, time
from pathlib import Path
import mpmath
from mpmath import iv, mp
resource.setrlimit(resource.RLIMIT_AS, (6*1024**3, 6*1024**3))
resource.setrlimit(resource.RLIMIT_CPU, (7000, 7000))
START=time.monotonic()
iv.dps=240
mp.dps=260
BASE=Path(__file__).resolve().parent
EXPECTED={'INPUT.json':'3e8c9e4d5247650fbece3681834925069a745653699b919d01ddb8fa8b7216c6',
'CANDIDATE.json':'9e68a52db275b1aa4f0b42ee1ebcc6a523b3fa80930508d98ce3e7af39ec5c10'}
inputs={}
for name,digest in EXPECTED.items():
    raw=(BASE/name).read_bytes()
    assert hashlib.sha256(raw).hexdigest()==digest
    inputs[name]=json.loads(raw)
cand=inputs['CANDIDATE.json']
I=iv.mpf
ZERO=I(0)
ONE=I(1)
def excludes(x): return bool(x>0) or bool(x<0)
def upper_abs(x): return abs(x).b
def midmp(x): return (mp.mpf(x._mpi_[0])+mp.mpf(x._mpi_[1]))/2
def brief(x): return [mp.nstr(mp.mpf(t),16) for t in x._mpi_]
def encode(x):
    # Endpoint tuples give exact binary rational endpoints (sign,mantissa,exponent,bits).
    return {'approx':brief(x),'binary_endpoints':x._mpi_}
checks={}
def note(name,x): checks[name]=encode(x)
class Jet:
    """First-order forward differentiation in (h,k,r), over intervals."""
    def __init__(self,v,g=None):
        self.v=v if hasattr(v,'_mpi_') else I(v)
        self.g=g if g is not None else [ZERO]*3
    def __add__(self,o):
        o=asjet(o); return Jet(self.v+o.v,[a+b for a,b in zip(self.g,o.g)])
    __radd__=__add__
    def __neg__(self): return Jet(-self.v,[-a for a in self.g])
    def __sub__(self,o): return self+-asjet(o)
    def __rsub__(self,o): return asjet(o)+-self
    def __mul__(self,o):
        o=asjet(o); return Jet(self.v*o.v,[a*o.v+self.v*b for a,b in zip(self.g,o.g)])
    __rmul__=__mul__
    def __truediv__(self,o):
        o=asjet(o); assert excludes(o.v)
        return Jet(self.v/o.v,[(a-self.v/o.v*b)/o.v for a,b in zip(self.g,o.g)])
    def __rtruediv__(self,o): return asjet(o)/self

def asjet(x): return x if isinstance(x,Jet) else Jet(x)
def value(x): return x.v if isinstance(x,Jet) else x

def solve(A,b):
    """Interval Gaussian elimination: all pivots exclude zero, proving invertibility.
    For Jet entries the same arithmetic differentiates the unique solution.
    """
    a=[list(row)+[rhs] for row,rhs in zip(A,b)]
    n=len(a)
    for k in range(n):
        p=max(range(k,n),key=lambda j: abs(midmp(value(a[j][k]))))
        a[k],a[p]=a[p],a[k]
        assert excludes(value(a[k][k])), ('pivot',k)
        for i in range(k+1,n):
            t=a[i][k]/a[k][k]
            for j in range(k+1,n+1): a[i][j]=a[i][j]-t*a[k][j]
            a[i][k]=0
    x=[0]*n
    for i in range(n-1,-1,-1):
        x[i]=(a[i][n]-sum(a[i][j]*x[j] for j in range(i+1,n)))/a[i][i]
    return x

def determinant(A):
    a=[list(row) for row in A]; det=ONE; n=len(a)
    for k in range(n):
        p=max(range(k,n),key=lambda j: abs(midmp(a[j][k])))
        if p!=k: a[k],a[p]=a[p],a[k]; det=-det
        assert excludes(a[k][k])
        det=det*a[k][k]
        for i in range(k+1,n):
            t=a[i][k]/a[k][k]
            for j in range(k+1,n): a[i][j]=a[i][j]-t*a[k][j]
    return det

def left(h,k,r):
    S=[1]; X=[r]; Y=[]
    for j in range(19):
        y=sum(math.comb(j,l)*S[l]*X[j-l] for l in range(j+1))
        Y.append(y); S.append(-y); X.append(k*X[j]+y)
    q=[S[0]+h]+S[1:]
    z=[]
    for j in range(19):
        z.append((Y[j]-sum(math.comb(j,l)*q[l]*z[j-l] for l in range(1,j+1)))/q[0])
    return S,X,Y,z

def contact(coords,ad=False):
    if ad: coords=[Jet(x,[ONE if i==j else ZERO for j in range(3)]) for i,x in enumerate(coords)]
    S,X,Y,z=left(*coords)
    A=[[z[j+i] for i in range(5)]+[-Y[j+i] for i in range(5)] for j in range(10)]
    dv=solve(A,[-z[j+5] for j in range(10)])
    d,v=dv[:5],dv[5:]
    R=[z[j+5]+sum(d[i]*z[j+i]-v[i]*Y[j+i] for i in range(5)) for j in range(14)]
    return R,d,v,(S,X,Y,z)

center=[I(s) for s in cand['center']]
radius=I(cand['radius'])
box=[x+I([-1,1])*radius for x in center]
R0,_,_,_=contact(center,True)
RB,dj,vj,_=contact(box,True)
J=[[RB[i+10].g[j] for j in range(3)] for i in range(3)]
Jmid=mp.matrix([[midmp(x) for x in row] for row in J])
pre=Jmid**-1
# Finite decimals are exact rational preconditioner definitions (enclosed by iv).
C=[[I(mp.nstr(pre[i,j],230)) for j in range(3)] for i in range(3)]
E=[[(ONE if i==j else ZERO)-sum(C[i][k]*J[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
q=max([sum(upper_abs(x) for x in row) for row in E],key=midmp)
beta=max([upper_abs(sum(C[i][j]*R0[j+10].v for j in range(3))) for i in range(3)],key=midmp)
assert q<1
assert beta+q*radius<radius
assert excludes(determinant(C))
note('contraction_q',q); note('center_correction_beta',beta)
note('self_map_radius_ratio',(beta+q*radius)/radius)
# Thus x -> x - C F(x) contracts the entire supplied rational box into itself.
# Every fixed point lies in the certified contraction image about the original center.
root_radius=beta+q*radius
note('contracted_root_radius',root_radius)
box=[x+I([-1,1])*root_radius.b for x in center]
R,d,v,(S,X,Y,z)=contact(box)
note('R13',R[13]); assert excludes(R[13])
for i in range(5): note('d'+str(i),d[i]); note('v'+str(i),v[i])

# Uniform existence and uniqueness of three real poles by sign change + derivative.
def poly(x):
    t=ONE
    for a in reversed(d): t=t*x+a
    return t
def dp(x): return 5*x**4+sum(i*d[i]*x**(i-1) for i in range(1,5))
poles=[]
for idx,s in enumerate(cand['real_pole_centers']):
    m=I(s); eps=I(cand['pole_radius']); interval=m+I([-1,1])*eps
    lo,hi=poly(m-eps),poly(m+eps)
    assert (lo>0 and hi<0) or (lo<0 and hi>0)
    derivative=dp(interval); assert excludes(derivative)
    poles.append(interval)
    note('pole'+str(idx+1),interval)
    note('pole'+str(idx+1)+'_left_value',lo); note('pole'+str(idx+1)+'_right_value',hi)
    note('pole'+str(idx+1)+'_derivative',derivative)
l1,l2,l3=poles
assert l1.b<l2.a and l2.b<l3.a
s1=l1+l2+l3; s2=l1*l2+l1*l3+l2*l3
q1=d[4]+s1; q0=d[3]+s1*q1-s2
alpha=-q1/2; eta=q0-alpha**2
assert eta>0 and alpha<0
note('alpha',alpha); note('eta',eta)
delta=iv.sqrt(eta/3)
a=(2*alpha+l3)/3; b=(l3-alpha)/3+delta; c=(l3-alpha)/3-delta
K=[[a,b,c,ZERO,ZERO],[c,a,b,ZERO,ZERO],[b,c,a,ZERO,ZERO],
   [I(3)/8,I(17)/40,I(1)/5,l1,ZERO],[ZERO,ZERO,ZERO,ONE,l2]]
for i in range(5):
    for j in range(5):
        if i!=j: assert K[i][j]>=0
cols=[sum(K[i][j] for i in range(5)) for j in range(5)]
assert all(x<0 for x in cols)
assert all(x<0 for x in poles)
for i,x in enumerate(cols): note('column_sum'+str(i),x)
for name,x in [('cycle_a',a),('cycle_b',b),('cycle_c',c)]: note(name,x)

def mv(A,x): return [sum(a*b for a,b in zip(row,x)) for row in A]
def transpose(A): return [list(row) for row in zip(*A)]
def dot(x,y): return sum(a*b for a,b in zip(x,y))
g=v[4]; assert g>0; note('g',g)
p=[x/g for x in v]
H=[]
for j in range(5): H.append(p[4-j]-sum(d[5-t]*H[j-t] for t in range(1,j+1)))
powers=[[ONE,ZERO,ZERO,ZERO,ZERO]]
for j in range(4): powers.append(mv(K,powers[-1]))
Qc=transpose(powers)
Cv=solve(powers,H)
O=[Cv]
for j in range(4): O.append(mv(transpose(K),O[-1]))
Z=solve(O,[z[j]/g-sum(H[j-1-l]*Y[l] for l in range(j)) for j in range(5)])
for name,x in [('det_Qc',determinant(Qc)),('det_O',determinant(O))]:
    assert excludes(x); note(name,x)
assert all(x>0 for x in Cv+Z)
for j,x in enumerate(Cv): note('C'+str(j+1),x)
for j,x in enumerate(Z): note('Z'+str(j+1),x)
# C e1 = H0 = v4/g = 1 is an exact identity, not an interval-zero test.
assert Cv[0].a<=1 and Cv[0].b>=1
h,k,r=box
assert 1+h>0 and r>0 and k<0
Rleft=100-1-r; Rright=100-(1+h)-sum(Z)
assert Rleft>0 and Rright>0
note('left_R0',Rleft); note('right_R0',Rright)
# Independent direct nonlinear ODE derivative recurrence (diagnostic corroboration).
def rhs_jets(K,Cv,beta,S0,Z0,n=18):
    ss=[S0]; xx=[Z0]; yy=[]
    for j in range(n+1):
        y=beta*sum(math.comb(j,l)*ss[l]*dot(Cv,xx[j-l]) for l in range(j+1))
        yy.append(y); ss.append(-y)
        dz=mv(K,xx[j]); dz[0]+=y; xx.append(dz)
    return yy
yright=rhs_jets(K,Cv,g,1+h,Z)
diffs=[Y[j]-yright[j] for j in range(19)]
for j,x in enumerate(diffs): note('direct_jet_difference_'+str(j),x)
assert all(x.a<=0 and x.b>=0 for x in diffs[:18]) # corroboration ONLY
assert excludes(diffs[18])
expected=(1+h)*R[13]; assert excludes(expected)
note('exact_formula_jet18_difference',expected)
# Reject altered transfer coefficient and altered initial state, with separated residuals.
epsilon=I('1e-6')
vbad=v[:]; vbad[0]+=epsilon
bad_R0=z[5]+sum(d[i]*z[i]-vbad[i]*Y[i] for i in range(5))
assert excludes(bad_R0); note('negative_control_v0_plus_1e-6_residual',bad_R0)
Zbad=Z[:]; Zbad[0]+=epsilon
bad_y0=Y[0]-g*(1+h)*dot(Cv,Zbad)
assert excludes(bad_y0); note('negative_control_Z1_plus_1e-6_y0_difference',bad_y0)
badcenter=center[:]; badcenter[0]+=I('1e-8')
badR,_,_,_=contact(badcenter)
assert excludes(badR[10]); note('negative_control_h_plus_1e-8_R10',badR[10])
result={'status':'PASS','input_sha256':EXPECTED,'precision_decimal':iv.dps,
'python':platform.python_version(),'mpmath':mpmath.__version__,
'checks':checks,'elapsed_seconds':time.monotonic()-START,
'peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
'proof_dependency':'STAGE1.md; low-order zero containment is not used to prove equality'}
print(json.dumps(result,indent=2))
