#!/usr/bin/env python3
"""Nonblind B verification; no producer code imported or executed.
All mathematical decisions use exact polynomial/rational arithmetic.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse, json, math, resource, time, hashlib
from pathlib import Path
resource.setrlimit(resource.RLIMIT_AS,(7*1024**3,7*1024**3))
resource.setrlimit(resource.RLIMIT_CPU,(7000,7000))
import sympy as sp
HERE=Path(__file__).resolve().parent
S,J,F,L,A,a,b,c,d,e,g,v,h=sp.symbols('S J F L A a b c d e g v h')
PARAM=(S,F,a,b,c,d,e,g,v,h)
def equal(x,y=0):
    assert sp.cancel(x-y)==0, str(x-y)[:300]
def save(name,data):
    data.update(status='PASS',seconds=time.monotonic()-START,rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (HERE/(name+'.json')).write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k:v for k,v in data.items() if k not in ('P','H','edges','raw','reduced','resultants')},ensure_ascii=False),flush=True)
def dt(expr):
    return -S*J*sp.diff(expr,S)+((S-a)*J+F)*sp.diff(expr,J)+(c*S-b)*J*sp.diff(expr,F)
def residual():
    z=S*J/(S+h)
    R=dt(dt(z))+d*dt(z)+e*z-g*dt(S*J)-v*S*J
    return sp.Poly(sp.cancel((S+h)**3*R/S),J).as_expr()
def bridge():
    P=residual(); V=-dt(P); q=S+h
    M=3*h*J+q*((g-1)*S+a-d+g*h)
    equal(sp.diff(P,F),-q*M)
    equal(sp.Poly(P,J).LC(),h*(h-S))
    equal(sp.Poly(V,J).LC(),-h*S)
    equal(P.subs(S,-h),2*h*h*J**3)
    f=A+b*L-c*S
    p=d-a+(1-g)*S
    k=sp.symbols('k')
    U=(b-e)*L+(v-c)*S+k
    H=p*J+U
    partial=lambda z:S*sp.diff(z,S)+sp.diff(z,L)
    flow=lambda z:partial(z)+(-(S-a)-f/J)*sp.diff(z,J)
    Q=sp.expand(U*(p*partial(U)-partial(p)*U)-(S-a)*p*p*U+f*p**3)
    identity=-p*(H-U)*flow(H)+partial(p)*(H**2-2*H*U)+p*partial(U)*H-(S-a)*p*p*H
    equal(Q,identity)
    equal(-J*flow(H),P.subs({h:0,F:f})/S**2)
    poly=sp.Poly(Q,S,L)
    equal(poly.coeff_monomial(S*L**2),-(b-e)**2*(1-g))
    equal(poly.coeff_monomial(S**3*L),(1-g)**2*(-(b-e)+(1-g)*b))
    support=[list(z) for z in poly.monoms()]
    assert set(poly.monoms()) <= {(0,0),(0,1),(1,0),(1,1),(1,2),(2,0),(2,1),(3,0),(3,1),(4,0)}
    assert sp.expand(Q+1-identity)!=0
    # CH and first integral are valid without controllability/observability.
    k11,k12,k21,k22,c2,x1,x2,yy=sp.symbols('k11 k12 k21 k22 c2 x1 x2 yy')
    K=sp.Matrix([[k11,k12],[k21,k22]]); C=sp.Matrix([[1,c2]]); X=sp.Matrix([x1,x2])
    aa=-sp.trace(K);bb=K.det();cc=aa+(C*K)[0]
    assert (K*K+aa*K+bb*sp.eye(2)).applyfunc(sp.expand)==sp.zeros(2)
    remainder=(C*(K+aa*sp.eye(2))*X)[0]
    equal(remainder,-k22*x1+k12*x2+c2*k21*x1-c2*k11*x2)
    equal((C*(K+aa*sp.eye(2))*(K*X+sp.Matrix([yy,0])))[0],-bb*(C*X)[0]+cc*yy)
    save('bridge',dict(P=str(P),h0_support=support,h0_identity=True,CH_without_minimality=True,
        negative_control='Q0+1 rejected',residual_J_degrees=[sp.degree(P,J).__int__(),sp.degree(V,J).__int__()]))
def universal():
    P=residual();V=sp.expand(-dt(P))
    from sympy.polys.rings import ring
    from functools import lru_cache
    ring0,*_=ring(PARAM,sp.ZZ)
    pc=[ring0.from_expr(z) for z in sp.Poly(P,J).all_coeffs()]
    vc=[ring0.from_expr(z) for z in sp.Poly(V,J).all_coeffs()]
    matrix=[[ring0.zero]*k+pc+[ring0.zero]*(3-k) for k in range(4)]
    matrix += [[ring0.zero]*k+vc+[ring0.zero]*(2-k) for k in range(3)]
    print('Computing universal Sylvester determinant by cached Laplace minors',flush=True)
    @lru_cache(None)
    def minor(columns):
        row=7-len(columns)
        if not columns:return ring0.one
        value=ring0.zero
        for pos,col in enumerate(columns):
            if matrix[row][col]:
                term=matrix[row][col]*minor(columns[:pos]+columns[pos+1:])
                value=value+term if pos%2==0 else value-term
        return value
    det_ring=minor(tuple(range(7)))
    print('Universal determinant computed:',len(det_ring),'monomials',flush=True)
    divisor=ring0.from_expr((S+h)**3*F)
    quotient,rem=det_ring.div(divisor);assert not rem
    assert quotient*divisor==det_ring
    print('Exact quotient and reconstruction computed',flush=True)
    terms=list(quotient.items())
    maxima={}
    for mon,co in terms: maxima[mon[1]]=max(maxima.get(mon[1],-1),mon[0])
    # Substitution F=A+bL-cS: bound all contributions, not just existing terms.
    support_bound={ell:max(i+k-ell for k,i in maxima.items() if k>=ell) for ell in range(max(maxima)+1)}
    assert support_bound=={0:19,1:17,2:16,3:12,4:8,5:4}
    edges={}
    for ell in (2,3,4,5):
        indexes=range(5) if ell==5 else range(9) if ell==4 else (10,11,12) if ell==3 else (14,15,16)
        for i in indexes:
            value=0
            for mon,co in terms:
                extraS=i-mon[0];extraA=mon[1]-ell-extraS
                if min(extraS,extraA)>=0:
                    value+=co*sp.prod(z**n for z,n in zip(PARAM[2:],mon[2:]))*sp.binomial(mon[1],ell)*sp.binomial(mon[1]-ell,extraS)*b**ell*(-c)**extraS*A**extraA
            edges[f'{i},{ell}']=str(sp.factor(value))
    top=sum(co*sp.prod(z**n for z,n in zip(PARAM,mon))/F**5 for mon,co in terms if mon[1]==5)
    equal(top,108*h**7*(S*S+h*S+h*h)**2)
    for i,ell,value in [(8,4,864*h**5*(g-1)**2*b**4),(12,3,652*h**3*(g-1)**4*b**3),(16,2,-8*h*(g-1)**6*b*b)]:
        equal(sp.sympify(edges[f'{i},{ell}']),value)
    assert sp.expand(top-108*h**7*(S*S+h*S+h*h)**2+1)!=0
    save('universal',dict(support_bound=support_bound,quotient_terms=len(terms),edges=edges,
        exact_division=True,nonzero_top_for_h_b_nonzero=True,negative_control='top coefficient +1 rejected'))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['bridge','universal']);args=parser.parse_args()
    START=time.monotonic();globals()[args.stage]()
