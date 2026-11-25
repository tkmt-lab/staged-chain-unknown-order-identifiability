#!/usr/bin/env python3
"""Non-blind separate implementation; reads no researcher code or data.
Exact rational differentiation, explicit Leibniz determinant enumeration,
and sparse polynomial arithmetic. Run from anywhere; output is JSON stdout.
"""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
import sys, time, json, resource, itertools, hashlib
from pathlib import Path
resource.setrlimit(resource.RLIMIT_AS, (6*1024**3, 6*1024**3))
resource.setrlimit(resource.RLIMIT_CPU, (7000, 7000))
import sympy as s
from sympy.polys.rings import ring
START = time.monotonic()
results = {}
def record(name, value=True):
    results[name] = value
    print(name, 'OK', round(time.monotonic()-START, 3), file=sys.stderr, flush=True)
def eq(name, got, want):
    if s.cancel(got-want) != 0:
        raise AssertionError(name)
    record(name)
def reject(name, got, wrong):
    if s.cancel(got-wrong) == 0:
        raise AssertionError('negative control accepted: '+name)
    record(name)
S,F,J,a,b,c,h,d2,d1,d0,g,v,w = s.symbols('S F J a b c h d2 d1 d0 g v w')
q=S+h
flow=(-S*J,(c*S-b)*J,(S-a)*J+F)
def D(p):
    return s.cancel(sum(s.diff(p,x)*dx for x,dx in zip((S,F,J),flow)))
y=S*J
z=y/q
y1=D(y); y2=D(y1)
z1=D(z); z2=D(z1); z3=D(z2)
R=z3+d2*z2+d1*z1+d0*z-g*y2-v*y1-w*y
P=s.cancel(q**4*R/S).expand()
assert s.denom(P)==1
Q=-D(P).expand()
assert s.denom(Q)==1
record('residual_degrees', [s.degree(P,J),s.degree(Q,J),s.degree(P,F)])
eq('F_quadratic',s.Poly(P,F).coeff_monomial(F**2),-3*h*q**2)
eq('q_primitive',P.subs(S,-h),-6*h**3*J**4)
alpha=s.Poly(P,F).coeff_monomial(F**2)
beta=s.Poly(P,F).coeff_monomial(F)
gamma=P.subs(F,0)
disc=s.Poly((beta**2-4*alpha*gamma).expand(),J)
assert disc.degree()==4
eq('discriminant_top',disc.nth(4),24*h**2*q**2*(S*S-S*h+h*h))
eq('T_discriminant',s.discriminant(S*S-S*h+h*h,S),-3*h*h)
eq('T_at_minus_h',(S*S-S*h+h*h).subs(S,-h),3*h*h)
eq('J_top_P',s.Poly(P,J).nth(4),-h*(S*S-4*S*h+h*h))
eq('J_top_Q',s.Poly(Q,J).nth(5),2*S*h*(-S+2*h))
eq('nondivisibility_remainder',s.rem(2*S*(S-2*h),S*S-4*S*h+h*h,S),4*S*h-2*h*h)
reject('wrong_discriminant',disc.nth(4),24*h**2*q**2*(S*S+S*h+h*h))
reject('wrong_J_top',s.Poly(Q,J).nth(5),2*S*h*(-S+h))
# q-valuations are computed by shifting S=q-h and inspecting every coefficient.
qq=s.symbols('qq')
vals=[]
for poly in (P,Q):
    row=[]
    for coeff in reversed(s.Poly(poly,J).all_coeffs()):
        shifted=s.Poly(coeff.subs(S,qq-h).expand(),qq)
        row.append(min(k[0] for k,val in shifted.terms() if val))
    vals.append(row)
assert vals==[[2,2,1,1,0],[2,1,1,0,0,0]]
record('q_valuations',vals)
assert min(2*k+i for i,k in enumerate(vals[0]))>=4
assert min(2*k+i for i,k in enumerate(vals[1]))>=3
eq('F_common_root_P',P.subs({F:0,J:0}),0)
eq('F_common_root_Q',Q.subs({F:0,J:0}),0)
record('q6F_factor_exponents', {'x_rescaling':4*5+3*4-4*5,'q':6,'F':1})
# Explicit Sylvester matrix; descending coefficients, all permutations via DFS.
def determinant_leibniz(mat, one, multiply, add, zero):
    n=len(mat); total=zero; leaves=0
    def visit(row,cols,product,sign):
        nonlocal total,leaves
        if row==n:
            total=add(total,product if sign==1 else -product); leaves+=1; return
        for col in range(n):
            if col in cols or not mat[row][col]: continue
            p=multiply(product,mat[row][col])
            if p:
                visit(row+1,cols+(col,),p,sign*(-1 if sum(k>col for k in cols)%2 else 1))
    visit(0,(),one,1)
    return total,leaves
pc=s.Poly(P,J).all_coeffs(); qc=s.Poly(Q,J).all_coeffs()
# Tropical maximum convolution, explicitly enumerate every nonzero permutation.
def bound_coeff(expr):
    out={}
    for (si,fi),val in s.Poly(expr,S,F).terms():
        if val: out[fi]=max(si,out.get(fi,-1))
    return out
def matrix_of(pc,qc,zero):
    m,n=len(pc)-1,len(qc)-1
    return [[zero]*i+pc+[zero]*(n-1-i) for i in range(n)]+[[zero]*i+qc+[zero]*(m-1-i) for i in range(m)]
bm=matrix_of(list(map(bound_coeff,pc)),list(map(bound_coeff,qc)),{})
def maxprod(p,r):
    out={}
    for j,x in p.items():
        for k,y0 in r.items(): out[j+k]=max(out.get(j+k,-1),x+y0)
    return out
bounds={}; support_leaves=0
def support_visit(row,cols,p):
    global support_leaves
    if row==len(bm):
        support_leaves+=1
        for k,val in p.items(): bounds[k]=max(bounds.get(k,-1),val)
        return
    for col in range(len(bm)):
        if col not in cols and bm[row][col]:
            support_visit(row+1,cols+(col,),maxprod(p,bm[row][col]))
support_visit(0,(),{0:0})
assert [bounds[k] for k in range(1,11)]==[44,42,40,38,36,34,30,26,22,18]
assert max(bounds)==10
record('raw_support',{'bounds':bounds,'permutations':support_leaves})
# Highest weight layers: (eps,F)-truncated ring with S and all parameters.
# No parameter specialization or divisions in determinant enumeration.
ep=s.symbols('ep')
params=(S,a,b,c,h,d2,d1,d0,g,v,w)
pr,*gens=ring((ep,F)+params,s.QQ)
def topcoeff(expr,degree,jpower):
    out=pr.zero
    for (si,fi),val in s.Poly(expr,S,F).terms():
        if not val: continue
        deficit=degree-jpower-si-2*fi
        assert deficit>=0
        if deficit<3 and fi<4:
            out+=pr.from_expr(val*S**si)*gens[0]**deficit*gens[1]**fi
    return out
pm=[topcoeff(x,7,4-i) for i,x in enumerate(pc)]
qm=[topcoeff(x,8,5-i) for i,x in enumerate(qc)]
def truncmul(p,r):
    out={}
    for u,cu in p.items():
        for vv,cv in r.items():
            if u[0]+vv[0]>=3 or u[1]+vv[1]>=4: continue
            mon=tuple(x+y for x,y in zip(u,vv))
            out[mon]=out.get(mon,0)+cu*cv
    return pr.from_dict({k:v0 for k,v0 in out.items() if v0})
# Small exact determinant with a zero diagonal pivot; compare with SymPy's
# generic determinant, then project to the same quotient ring.
control=s.Matrix([[0,1+ep,F],[ep*F,S,1],[1,F+ep,a]])
cm=[[pr.from_expr(control[i,j]) for j in range(3)] for i in range(3)]
cd,_=determinant_leibniz(cm,pr.one,truncmul,lambda p,r:p+r,pr.zero)
full=pr.from_expr(control.det())
projected=pr.from_dict({mon:cc for mon,cc in full.items() if mon[0]<3 and mon[1]<4})
assert cd==projected and cd
assert cd!=projected+pr.one
record('zero_pivot_determinant_control')
record('top_matrix_ready')
top,leaves=determinant_leibniz(matrix_of(pm,qm,pr.zero),pr.one,truncmul,lambda p,r:p+r,pr.zero)
expected=-288*ep**2*F**3*S**39*h*(g-1)**7*(a*g-7*c*g+6*c-d2+v)
assert top==pr.from_expr(expected)
assert top!=pr.from_expr(-expected)
assert top!=pr.from_expr(expected+ep**2*F**2*S**41)
record('top_negative_controls')
record('universal_top',{'formula':str(expected),'nonzero_permutations':leaves,'monomials':len(top)})
rows=[bounds[k]-6 for k in range(1,11)]
for k in (1,2,3):
    # For each F degree, remove zero monomials within computed top layers.
    while (e:=47-rows[k-1]-6-2*k) in range(3):
        if any(m[0]==e and m[1]==k for m in top): break
        rows[k-1]-=1
assert rows==[36,34,33,32,30,28,24,20,16,12]
substituted=[max(rows[m]+m-ell for m in range(ell,len(rows))) for ell in range(len(rows))]
assert substituted==rows and sum(x+1 for x in rows)==275
record('final_support',{'rows':rows,'dimension':sum(x+1 for x in rows),'B':sum(x+1 for x in rows)+3})
# h=0 is reconstructed with an independent integrated residual in (S,L,J).
L,A,H,t,u=s.symbols('L A H t u')
f=A+b*L-c*S
def Dt(p):
    return s.expand(-S*J*s.diff(p,S)-J*s.diff(p,L)+((S-a)*J+f)*s.diff(p,J))
y=S*J
W=Dt(Dt(J))+(a+u)*Dt(J)+d1*J-d0*L-(1-t)*Dt(y)-v*y+w*S-H
W=s.expand(W)
W_claim=-t*S*J**2+(t*S**2+(-a*t+u+c-v)*S+d1-b-a*u)*J+(t*S+u)*f-d0*L+w*S-H
eq('integrated_W',W,W_claim)
R0=Dt(Dt(Dt(J)))+(a+u)*Dt(Dt(J))+d1*Dt(J)+d0*J-(1-t)*Dt(Dt(y))-v*Dt(y)-w*y
eq('W_derivative',Dt(W),R0)
reject('wrong_W_d0_sign',Dt(W+2*d0*L),R0)
eq('W_L_coefficient',s.Poly(W,L).nth(1),b*t*S+b*u-d0)
eq('W_J2',s.Poly(W,J).nth(2),-t*S)
VW=-Dt(W)
eq('VW_J3',s.Poly(VW,J).nth(3),-t*S)
# Full h=0 determinant, 5 by 5, directly by permutations (no 13-term formula).
hr,*hg=ring((S,L,a,b,c,A,H,t,u,d1,d0,v,w),s.QQ)
hpc=[hr.from_expr(x) for x in s.Poly(W,J).all_coeffs()]
hqc=[hr.from_expr(x) for x in s.Poly(VW,J).all_coeffs()]
E0,h_leaves=determinant_leibniz(matrix_of(hpc,hqc,hr.zero),hr.one,lambda p,r:p*r,lambda p,r:p+r,hr.zero)
assert E0
assert all(mon[0]>=1 and mon[7]>=1 for mon in E0)
quot=E0.exquo(hr.from_expr(t*S))
hrows={}
for mon,val in quot.items(): hrows[mon[1]]=max(hrows.get(mon[1],-1),mon[0])
assert hrows=={0:8,1:7,2:5,3:4}
record('h0_resultant',{'quotient_rows':hrows,'dimension':sum(v0+1 for v0 in hrows.values()),'permutations':h_leaves,'E_monomials':len(E0)})
# g=1, including vanishing linear pivot: division-free identity.
x,j0,k,lam,nu=s.symbols('x j0 k lam nu')
p=x*S+j0; U=k*L+lam*S+nu
part=lambda expr:s.expand(S*s.diff(expr,S)+s.diff(expr,L))
Q0=s.expand(U*(p*part(U)-part(p)*U)-(S-a)*p**2*U+f*p**3)
Wlin=p*J+U
partialW=part(Wlin)+s.diff(Wlin,J)*(-(S-a)-f/J)
identity=-p*(Wlin-U)*partialW+part(p)*(Wlin**2-2*Wlin*U)+p*part(U)*Wlin-(S-a)*p**2*Wlin
eq('g1_division_free_identity',identity,Q0)
eq('g1_matching_W',W.subs(t,0),Wlin.subs({x:u+c-v,j0:d1-b-a*u,k:b*u-d0,lam:w-u*c,nu:u*A-H}))
support0=set(s.Poly(Q0,S,L).monoms())
allowed={(0,0),(0,1),(1,0),(1,1),(1,2),(2,0),(2,1),(3,0),(3,1),(4,0)}
assert support0<=allowed
eq('g1_SL2',s.Poly(Q0,S,L).coeff_monomial(S*L**2),-k*k*x)
eq('g1_S3L',s.Poly(Q0,S,L).coeff_monomial(S**3*L),x*x*(-k+b*x))
reject('g1_wrong_SL2',s.Poly(Q0,S,L).coeff_monomial(S*L**2),k*k*x)
record('g1_space',{'support':sorted(support0),'dimension_bound':len(allowed)})
# Independent CH matrix checks use symbolic 2x2 and 3x3 matrices.
for n in (2,3):
    K=s.Matrix(n,n,s.symbols('k0:'+str(n*n)))
    C=s.Matrix([[1]+list(s.symbols('c1:'+str(n)))])
    e=s.eye(n)[:,0]; zz=s.symbols('zz')
    coeff=K.charpoly(zz).all_coeffs()
    power=s.zeros(n)
    for aa in coeff: power=power*K+aa*s.eye(n)
    assert power.applyfunc(s.expand)==s.zeros(n)
    numerator=s.expand((C*(zz*s.eye(n)-K).adjugate()*e)[0])
    assert s.Poly(numerator,zz).LC()==1
    # Derivative recurrence for arbitrary initial state X and arbitrary input jets.
    X=s.Matrix(s.symbols('X0:'+str(n)))
    yy=s.symbols('input0:'+str(n))
    Xjets=[X]
    for i in range(n): Xjets.append(K*Xjets[-1]+e*yy[i])
    lhs=sum(coeff[i]*(C*Xjets[n-i])[0] for i in range(n+1))
    rhs=sum(s.Poly(numerator,zz).nth(i)*yy[i] for i in range(n))
    eq('CH_arbitrary_initial_'+str(n),s.expand(lhs),s.expand(rhs))
    if n==2:
        aa=coeff[1]
        eq('positive_F_matrix',s.expand((C*(K+aa*s.eye(n))*X)[0]),
           (C*s.Matrix([[-K[1,1],K[0,1]],[K[1,0],-K[0,0]]])*X)[0])
    record('CH_'+str(n), {'degree':n,'numerator_degree':s.degree(numerator,zz)})
# The closed RHS is checked directly by differentiating y=q*z.
qq,zz,zz1,zz2,zz3=s.symbols('qq zz zz1 zz2 zz3')
Dc=lambda p:s.expand(-qq*zz*s.diff(p,qq)+zz1*s.diff(p,zz)+zz2*s.diff(p,zz1)+zz3*s.diff(p,zz2))
yy=qq*zz
eq('closed_y1',Dc(yy),qq*(zz1-zz*zz))
eq('closed_y2',Dc(Dc(yy)),qq*(zz**3-3*zz*zz1+zz2))
record('jet_orders',{'P_min_at_B278':278-2,'P_max_h_nonzero':275,'W_min_at_B30':30-1,'W_max_E_nonzero':28,'W_min_at_B12':12-1})
results.update(status='PASS',python=sys.version,sympy=s.__version__,seconds=time.monotonic()-START,peak_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
print(json.dumps(results,ensure_ascii=False,indent=2,default=str))
