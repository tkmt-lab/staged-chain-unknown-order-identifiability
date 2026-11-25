#!/usr/bin/env python3
"""Exact resultant certificates checked by integer determinants; custom
finite-field Euclid and rational Sturm arithmetic, independent of producer code.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import json,resource,time,hashlib
from fractions import Fraction as Q
from pathlib import Path
resource.setrlimit(resource.RLIMIT_AS,(7*1024**3,7*1024**3))
resource.setrlimit(resource.RLIMIT_CPU,(7000,7000))
import sympy as s
HERE=Path(__file__).resolve().parent;SOURCE=HERE.parent/'b22data'
h,u,lam=s.symbols('h u lam')
def prim(z,*vv):
    p=s.Poly(z,*vv).clear_denoms()[1].primitive()[1]
    return p if p.LC()>0 else -p

def determinant(matrix):
    """Integer fraction-free Gaussian elimination, with row pivoting."""
    m=[list(row) for row in matrix];n=len(m);prev=1;sign=1
    for k in range(n-1):
        pivot=next((i for i in range(k,n) if m[i][k]),None)
        if pivot is None:return 0
        if pivot!=k:m[k],m[pivot]=m[pivot],m[k];sign=-sign
        pk=m[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                val=m[i][j]*pk-m[i][k]*m[k][j]
                quo,rem=divmod(val,prev);assert rem==0
                m[i][j]=quo
            m[i][k]=0
        prev=pk
    return sign*m[-1][-1]
def sylvester(f,g):
    n,m=len(f)-1,len(g)-1
    return [[0]*i+f+[0]*(m-i-1) for i in range(m)]+[[0]*i+g+[0]*(n-i-1) for i in range(n)]
def evaluate(coefficients,value):
    out=0
    for z in coefficients:out=out*value+z
    return out

def trim(p):
    while p and not p[-1]:p.pop()
    return p
def remainder(p,q,prime=None):
    p=list(p);q=list(q)
    while len(p)>=len(q) and p:
        shift=len(p)-len(q)
        v=p[-1]*pow(q[-1],-1,prime)%prime if prime else p[-1]/q[-1]
        for j,z in enumerate(q):
            p[shift+j]-=v*z
            if prime:p[shift+j]%=prime
        trim(p)
    return p
def modgcd(p,q,prime):
    p=trim([int(z)%prime for z in p]);q=trim([int(z)%prime for z in q])
    while q:p,q=q,remainder(p,q,prime)
    return p
def variations(signs):return sum(a!=b for a,b in zip(signs,signs[1:]))

START=time.monotonic()
assert determinant([[0,2],[3,4]])==-6
assert determinant([[1,2],[2,4]])==0
assert determinant([[0,0],[0,1]])==0
assert remainder([Q(-1),Q(0),Q(0),Q(1)],[Q(-1),Q(1)])==[]
ours=json.loads((HERE/'constraints.json').read_text())
pol={k:s.Poly(s.sympify(z),u) for k,z in ours['reduced'].items()}
producer=json.loads((SOURCE/'upper_degree8_resultant.json').read_text())
resultants={};checks={}
for name in ('Q','R'):
    f,g=pol['E'],pol[name];n,m=f.degree(),g.degree()
    hf=max(s.degree(z,h) for z in f.all_coeffs());hg=max(s.degree(z,h) for z in g.all_coeffs())
    bound=int(m*hf+n*hg)
    fc=[[int(z) for z in s.Poly(co,h).all_coeffs()] for co in f.all_coeffs()]
    gc=[[int(z) for z in s.Poly(co,h).all_coeffs()] for co in g.all_coeffs()]
    candidate=prim(s.sympify(producer[name]['primitive']),h)
    candidate_coeff=[int(z) for z in candidate.all_coeffs()]
    assert candidate.degree()<=bound
    scale=None;points=[];changed_rejected=False
    for value in range(-(bound//2),bound-bound//2+1):
        matrix=sylvester([evaluate(co,value) for co in fc],[evaluate(co,value) for co in gc])
        det=determinant(matrix);expect=evaluate(candidate_coeff,value)
        if scale is None and expect:scale=Q(det,expect);assert scale
        if scale is not None:
            assert det==scale*expect
            if det!=scale*(expect+1):changed_rejected=True
        else:assert det==0
        points.append(value)
    assert len(points)==bound+1 and scale and changed_rejected
    resultants[name]=candidate
    checks[name]=dict(matrix_size=n+m,degree_bound=bound,points=len(points),first=points[0],last=points[-1],
        scale=str(scale),actual_degree=candidate.degree(),negative_control='primitive resultant +1 rejected')
    print('Resultant',name,'verified at',len(points),'exact points',flush=True)
# Reconstruct the exceptional octic from F(p7) and the lambda pivot polynomial.
f7=s.sympify(ours['raw']['F']).subs(lam,0)
pivot=760293059706180000*u*u-343854447300*u+167
p=s.Poly(f7,u);q=s.Poly(pivot,u)
mat=s.Matrix(sylvester(p.all_coeffs(),q.all_coeffs()))
G8=prim(mat.det(method='domain-ge'),h)
G4=s.Poly(2443344151248000*h**4-15437822400*h**3+3665200*h*h-4851*h+36,h)
known=prim(G4.as_expr()*G8.as_expr()**2,h)
cofactors={k:resultants[k].exquo(known) for k in resultants}
prime=1009
for cf in cofactors.values():assert int(cf.LC())%prime!=0
mg=modgcd(list(reversed(cofactors['Q'].all_coeffs())),list(reversed(cofactors['R'].all_coeffs())),prime)
assert len(mg)==1
# Positive quartic decomposition covers every real h, including h=0.
q1=2443344151248000*h*h-15437822400*h+1832600
q2=1832600*h*h-4851*h+36
assert s.expand(h*h*q1+q2-G4.as_expr())==0
quad_discriminants=[int(s.discriminant(z,h)) for z in (q1,q2)]
assert all(d<0 for d in quad_discriminants) and all(s.Poly(z,h).LC()>0 for z in (q1,q2))
# Sturm chain in ascending coefficient order, entirely using fractions.Fraction.
a0=[Q(int(z)) for z in reversed(G8.all_coeffs())]
a1=[i*a0[i] for i in range(1,len(a0))]
chain=[a0,a1]
while len(chain[-1])>1:
    r=remainder(chain[-2],chain[-1]);assert r
    chain.append([-z for z in r])
plus=[1 if z[-1]>0 else -1 for z in chain]
minus=[sg*((-1)**(len(z)-1)) for sg,z in zip(plus,chain)]
assert variations(minus)==variations(plus)==4
assert [len(z)-1 for z in chain]==list(range(8,-1,-1))
# Check recurrence once more, including no premature zero remainder.
for i in range(2,len(chain)):assert remainder(chain[i-2],chain[i-1])==[-z for z in chain[i]]
assert not remainder(chain[-2],chain[-1])
assert prim(s.sympify(json.loads((SOURCE/'upper_degree8_verify.json').read_text())['G8']),h)==G8
output=dict(status='PASS',resultant_checks=checks,G4=str(G4.as_expr()),G8=str(G8.as_expr()),
    gcd_certification=dict(common_factor=str(known.as_expr()),cofactor_degrees={k:z.degree() for k,z in cofactors.items()},
        modulus=prime,cofactor_degrees_preserved=True,modular_cofactor_gcd_degree=len(mg)-1),
    G4_discriminants=quad_discriminants,sturm_degrees=[len(z)-1 for z in chain],
    sturm_minus=minus,sturm_plus=plus,sturm_variations=[variations(minus),variations(plus)],
    sturm_chain=[[str(co) for co in z] for z in chain],
    seconds=time.monotonic()-START,rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(HERE/'elimination.json').write_text(json.dumps(output,indent=2)+'\n')
print('Elimination PASS',output['seconds'],'seconds',flush=True)
