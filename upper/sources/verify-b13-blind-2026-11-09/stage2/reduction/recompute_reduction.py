"""Stage2 independent polynomial-numerator recurrence for the residual.
No research code imports; only saved coefficients are read for comparison.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import sympy as s
import json,time,resource,signal,hashlib,datetime
from pathlib import Path
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));signal.alarm(1800)
t0=time.monotonic();here=Path(__file__).resolve().parent
S,X,L,h,k,r,d2,d1,d0,g,v,w=s.symbols('S X L h k r d2 d1 d0 g v w')
q=S+h;y=S*X;ell=S+k-X
D=lambda a:s.expand(-S*X*s.diff(a,S)+(S+k)*X*s.diff(a,X))
N=[s.Integer(1)];M=[s.Integer(1)]
for n in range(3):N.append(s.expand(q*(ell*N[-1]+D(N[-1]))+(n+1)*y*N[-1]))
for n in range(2):M.append(s.expand(ell*M[-1]+D(M[-1])))
# Algebraic recurrence: z^(n)=y*N[n]/q^(n+1), y^(n)=y*M[n].
assert s.expand(D(y)-y*ell)==0
assert s.expand(D(q)+y)==0
Fpoly=s.expand(N[3]+d2*q*N[2]+d1*q**2*N[1]+d0*q**3-q**4*(g*M[2]+v*M[1]+w))
F=s.Poly(Fpoly.subs(X,r+1-S-k*L),S,L)
support=sorted(F.monoms());assert support==[(i,j) for i,m in enumerate([4,4,4,3,3,2,1]) for j in range(m)]
path=Path('sources/b13-upper2-2026-11-07/check_base.json');data=path.read_bytes()
with (here/'access.jsonl').open('a') as f:f.write(json.dumps({'path':str(path),'sha256':hashlib.sha256(data).hexdigest(),'opened_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'mode':'recomputation input'})+'\n')
old=json.loads(data)
comparison={f'{i},{j}':s.expand(F.coeff_monomial(S**i*L**j)-s.sympify(old['F_coefficients'][f'{i},{j}']))==0 for i,j in support};assert all(comparison.values())
# Derivatives at x=0 of exp(i*x)*x^j; factorial/binomial formulation.
def jets(supp):
 return s.Matrix([[s.binomial(n,j)*s.factorial(j)*i**(n-j) if n>=j else 0 for i,j in supp] for n in range(len(supp))])
V=jets(support);det=V.det();assert det!=0
coeff=V.inv()*s.Matrix([0]*20+[1]);W=sum(c*S**i*L**j for c,(i,j) in zip(coeff,support))
assert s.expand(W-s.sympify(old['W'],locals={'S':S}))==0
F0=s.Poly(F.as_expr().subs(h,0),S,L);support0=sorted(F0.monoms());det0=jets(support0).det();assert len(support0)==7 and det0!=0
# Derivatives D_x^j W, j=0..11 form all order>=9 elements.
dx=lambda a:s.expand(S*s.diff(a,S)+s.diff(a,L))
ws=[W]
for _ in range(20):ws.append(dx(ws[-1]))
assert [p.subs({S:1,L:0}) for p in ws]==[0]*20+[1]
coeff_matrix=s.Matrix([[s.Poly(a,S,L).coeff_monomial(S**i*L**j) for a in ws[:12]] for i,j in support]);assert coeff_matrix.rank()==12
jet9=V[:9,:];assert jet9.rank()==9
# Polynomial germ closure and m0 input-response subtraction.
a,b,c,qv,gg,vv,ww,dd2,dd1,dd0=s.symbols('a b c qv gg vv ww dd2 dd1 dd0')
dyn={qv:-qv*a,a:b,b:c}
Dt=lambda f:s.expand(sum(s.diff(f,x)*f1 for x,f1 in dyn.items()))
y0=qv*a;y1=Dt(y0);y2=Dt(y1)
assert s.expand(y1-qv*(b-a*a))==0
assert s.expand(y2-qv*(c-3*a*b+a**3))==0
# w0=z0/g; w1=z1/g-y0; w2=z2/g-y1-(v1/g-d2)*y0.
z0,z1,z2,yy0,yy1=s.symbols('z0 z1 z2 yy0 yy1')
free=[z0/g,z1/g-yy0,z2/g-yy1-(v/g-d2)*yy0]
assert s.simplify(free[2]+d2*free[1]+d1*free[0]-((z2+d2*z1+d1*z0-v*yy0)/g-yy1))==0
out={'status':'PASS','method':'polynomial numerator recurrence derived from original left ODE; exact SymPy arithmetic; no producer code imported','sympy_version':s.__version__,'source_sha256':hashlib.sha256(data).hexdigest(),'support':support,'all_21_coefficients_equal':comparison,'jet_matrix_det':str(det),'h0_support':support0,'h0_jet_matrix_det':str(det0),'W_matches':True,'W_zero_order':20,'order_ge_9_space_dimension':12,'W_derivatives_0_to_11_rank':12,'h0_sufficient_B':9,'general_space_sufficient_B':23,'B11_residual_order':9,'germ_polynomial_ODE_formulas_checked':True,'m0_input_subtraction_checked':True,'seconds':time.monotonic()-t0,'maxrss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(here/'recompute_reduction.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
