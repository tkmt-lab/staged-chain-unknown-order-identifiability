"""Exact checker for the degree-eleven certificate, no generator imports.
Model columns use 1D logarithmic differentiation; W uses local series.
Resultants use fixed-size integer Sylvester matrices, with degree bounds.
The latter routine is adapted from the previous m5-upper/check.py.
"""
import argparse,json,math,itertools,time,resource,signal,hashlib,datetime
from pathlib import Path
import sympy as s
ap=argparse.ArgumentParser(); ap.add_argument('certificate',type=Path); ap.add_argument('--out',type=Path,required=True); args=ap.parse_args()
signal.alarm(1800); resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3)); start=time.monotonic()
C=json.loads(args.certificate.read_text()); H,u,A,c,e,x=s.symbols('H u A c e x'); S,L,h,k,r,d,d3,v3=s.symbols('S L h k r d4 d3 v3')
def require(test,msg):
 if not test: raise ValueError('REJECT: '+msg)
def dec(rec,vs):
 require(len({tuple(e) for e,v in rec})==len(rec),'duplicate monomial')
 require(all(len(e)==len(vs) and all(isinstance(a,int) and a>=0 for a in e) for e,v in rec),'exponents')
 return sum(s.Rational(v)*math.prod(t**i for t,i in zip(vs,e)) for e,v in rec)
def same(a,b,msg): require(s.cancel(a-b)==0,msg)
def prim(z,*vs): return s.Poly(z,*vs).clear_denoms()[1].primitive()[1].as_expr()
def associated(a,b,vs,msg):
 pa,pb=s.Poly(a,*vs),s.Poly(b,*vs); require(not pa.is_zero and not pb.is_zero,msg); same(a,pb.as_expr()*pa.LC()/pb.LC(),msg)
def nonzero_monomial(z):
 num,den=s.cancel(z).as_numer_denom()
 return all(len(s.Poly(f,H,u).terms())==1 and f!=0 and all(co.is_Rational for co in s.Poly(f,H,u).coeffs()) for f in [num,den])
def strip_hu(z):
 p=s.Poly(z,H,u,A); powers=[min(e[i] for e in p.monoms()) for i in range(2)]
 return prim(s.cancel(z/H**powers[0]/u**powers[1]),H,u,A)
# Generate all eleven columns from the first integral, not a recurrence in S,X.
y=S*(r+1-S-k*L); q=S+h
partial=lambda f:s.diff(f,S)+s.diff(f,L)/S
bb=[s.Integer(1)]; aa=[s.Integer(1)]
for j in range(5):
 bb.append(s.expand(-q*(partial(y)*bb[-1]+y*partial(bb[-1]))+(j+1)*y*bb[-1]))
 aa.append(s.expand(-partial(y)*aa[-1]-y*partial(aa[-1])))
cols=[bb[5]]+[s.expand(q**(5-j)*bb[j]) for j in range(5)]+[s.expand(-q**6*aa[j]) for j in range(5)]
polys=[s.Poly(z,S,L) for z in cols]
allowed={(i,j) for i in range(11) for j in range(6) if i+j<=10 and (i,j)!=(5,5)}
small={(5,0)}|{(i,j) for i in range(6,11) for j in range(11-i)}
require(set().union(*(set(p.monoms()) for p in polys))==allowed and len(allowed)==50,'support')
require(set().union(*(set(s.Poly(z.subs(h,0),S,L).monoms()) for z in cols))==small and len(small)==16,'h0 support')
mu=[max(j for ii,j in allowed if ii==i)+1 for i in range(11)]
w={}
for i,m in enumerate(mu):
 series=[s.Integer(1)]+[s.Integer(0)]*(m-1)
 for j,n in enumerate(mu):
  if i==j: continue
  fac=[(-1)**ell*s.binomial(n+ell-1,ell)*s.Rational(i-j)**(-n-ell) for ell in range(m)]
  series=[sum(series[a]*fac[ell-a] for a in range(ell+1)) for ell in range(m)]
 for j in range(m): w[i,j]=series[m-1-j]/s.factorial(j)
for n in range(50): same(sum(co*s.factorial(n)/s.factorial(n-j)*s.Integer(i)**(n-j) for (i,j),co in w.items() if j<=n),int(n==49),'W jets')
p=dec(C['layers']['p'],[x,H,u,A,c,e]); require(s.degree(p,x)<=11,'degree eleven')
for i in list(range(5))+list(range(6,11)):
 j=5 if i<5 else 10-i
 actual=polys[0].coeff_monomial(S**i*L**j)+(1-A/H**5)*polys[10].coeff_monomial(S**i*L**j)
 for f in polys[1:10]: same(f.coeff_monomial(S**i*L**j),0,'omitted top columns')
 same(w[i,j]*p.subs(x,i),actual.subs({h:1/H,k:1/u})*w[0,5]*H**5*u**5,'interpolation values')
same(p.subs(x,0),1,'normalization')
same(s.diff(p,c),s.prod(x-i for i in list(range(5))+list(range(6,11))),'free coefficient')
same(s.diff(p,e),x*s.diff(p,c),'second free coefficient')
# Apply p(D) to W; all derivatives computed from the saved polynomial.
def target(i,j): return s.expand(sum(s.diff(p,x,ell).subs(x,i)*s.binomial(j+ell,ell)*w[i,j+ell] for ell in range(mu[i]-j)))
rows=[]
for i in range(6):
 for idx in list(range(1,5))+list(range(6,10)): same(polys[idx].coeff_monomial(S**i*L**4),0,'omitted L4 columns')
 actual=(polys[0]+d*polys[5]+(1-A/H**5)*polys[10]).coeff_monomial(S**i*L**4)
 diff=s.cancel(actual.subs({h:1/H,k:1/u})*w[0,5]*H**5*u**5-target(i,4)); num,den=diff.as_numer_denom(); require(nonzero_monomial(den),'L4 denominator')
 row=[dec(z,[H,u,A]) for z in C['layers']['rows'][i]]; saved=row[0]*r+row[1]*d+row[2]*c+row[3]*e+row[4]
 ratio=s.cancel(saved.coeff(c)/s.diff(num,c)); require(nonzero_monomial(ratio),'L4 row scaling'); same(saved,num*ratio,'L4 row')
 rows.append(row)
M=s.Matrix(rows); pivot=M[[0,1,5],:3].det(); require(nonzero_monomial(pivot),'L4 pivot')
base_sub=[s.cancel(z) for z in M[[0,1,5],:3].inv()*(-M[[0,1,5],4]-M[[0,1,5],3]*e)]
el=C['elim']; erows=[]
for i in range(2,5):
 diff=s.cancel((M[i,:3]*s.Matrix(base_sub))[0]+M[i,3]*e+M[i,4]); num,den=diff.as_numer_denom()
 require(nonzero_monomial(den),'e row denominator')
 saved=dec(el['e_rows'][i-2],[H,u,A,e]); ratio=s.cancel(s.diff(saved,e)/s.diff(num,e))
 require(nonzero_monomial(ratio),'e row unit scaling'); same(saved,num*ratio,'e row')
 erows.append(saved)
locals0=dict(H=H,u=u,A=A)
bez=[s.sympify(v,locals=locals0) for v in el['bezout']]
gcd=sum(z*s.diff(row,e) for z,row in zip(bez,erows)); require(nonzero_monomial(gcd),'e coefficient Bezout unit')
es=s.cancel(-sum(z*row.subs(e,0) for z,row in zip(bez,erows))/gcd)
sub=[s.cancel(z.subs(e,es)) for z in base_sub]+[es]
for z,saved in zip(sub,el['sub']):
 same(z,s.sympify(saved,locals=locals0),'saved substitution'); require(nonzero_monomial(z.as_numer_denom()[1]),'substitution denominator')
ts=[]
for row,saved in zip(erows,el['T']):
 diff=s.cancel(row.subs(e,es)); num,den=diff.as_numer_denom(); require(nonzero_monomial(den),'T denominator')
 expected=strip_hu(num); same(expected,dec(saved,[H,u,A]),'T');ts.append(expected)
cr=C['cross']; t4=dec(cr['T'],[H,u,A]); require(s.degree(t4,A)==1,'core linear in A')
bez2=[s.sympify(v,locals=locals0) for v in cr['cofactor_bezout']]
for v in bez2: require(s.denom(v).is_Rational,'core Bezout polynomial')
same(bez2[0]*ts[0]+bez2[1]*ts[1],t4,'core Bezout')
pproj=dec(cr['P'],[H,u]); same(strip_hu(t4.coeff(A)*ts[2].subs(A,0)-ts[2].coeff(A)*t4.subs(A,0)),pproj,'P cross')
print('PASS model, W, top layers, L4 and global pivot coverage',flush=True)
rows3=[]
for i in range(7):
 for idx in [1,2,3,6,7,8]: same(polys[idx].coeff_monomial(S**i*L**3),0,'omitted L3 columns')
 actual=(polys[0]+d*polys[5]+d3*polys[4]+(1-A/H**5)*polys[10]+v3*polys[9]).coeff_monomial(S**i*L**3)
 diff=s.cancel(actual.subs({h:1/H,k:1/u,r:sub[0],d:sub[1]})*w[0,5]*H**5*u**5-target(i,3).subs({c:sub[2],e:sub[3]})); num,den=diff.as_numer_denom(); require(nonzero_monomial(den),'L3 denominator')
 row=[dec(z,[H,u,A]) for z in C['layer3']['rows'][i]]; saved=row[0]*d3+row[1]*v3+row[2]
 associated(num,saved,[d3,v3,H,u,A],'L3 row'); rows3.append(row)
M3=s.Matrix(rows3); pivot3=M3[:2,:2].det(); require(nonzero_monomial(pivot3),'L3 pivot')
q3=[]
for i in range(2,7):
 t=strip_hu(M3[[0,1,i],:].det(method='domain-ge')); same(t,dec(C['cross']['minors'][i-2],[H,u,A]),'L3 minor'); require(s.degree(t,A)<=2,'L3 A degree')
 a=t4.coeff(A); b=t4.subs(A,0)
 qq=strip_hu(t.coeff(A,2)*b**2-t.coeff(A)*a*b+t.subs(A,0)*a**2)
 expected=dec(C['cross']['Z'][i-2],[H,u]); same(qq,expected,'L3 cross'); q3.append(expected)
print('PASS L3 rows, minors and division-free A elimination',flush=True)
def bareiss(mat):
 mat=[row[:] for row in mat]; sign=1; prev=1; n=len(mat)
 for k0 in range(n-1):
  if mat[k0][k0]==0:
   row=next((i for i in range(k0+1,n) if mat[i][k0]),None)
   if row is None: return 0
   mat[k0],mat[row]=mat[row],mat[k0]; sign=-sign
  pivot=mat[k0][k0]
  for i in range(k0+1,n):
   for j in range(k0+1,n):
    val=mat[i][j]*pivot-mat[i][k0]*mat[k0][j]; require(val%prev==0,'Bareiss exact division'); mat[i][j]=val//prev
   mat[i][k0]=0
  prev=pivot
 return sign*mat[-1][-1]
def resultant_identity(f,g,expected,label):
 pf=s.Poly(f,H,u); pg=s.Poly(g,H,u); n=pf.degree(H); m=pg.degree(H)
 degree=m*pf.degree(u)+n*pg.degree(u); require(s.degree(expected,u)<=degree,label+' degree')
 def coefficients(p,N,at):
  arr=[0]*(N+1)
  for (eh,eu),co in p.terms(): require(co.is_Integer,'integer resultant data'); arr[N-eh]+=int(co)*at**eu
  return arr
 ep=s.Poly(expected,u); require(all(v.is_Integer for v in ep.all_coeffs()),'integer resultant')
 for ii in range(degree+1):
  at=0 if ii==0 else (ii+1)//2*(1 if ii%2 else -1)
  fc=coefficients(pf,n,at); gc=coefficients(pg,m,at)
  mat=[[0]*j+fc+[0]*(m-1-j) for j in range(m)]+[[0]*j+gc+[0]*(n-1-j) for j in range(n)]
  require(bareiss(mat)==int(ep.eval(at)),label+' determinant at '+str(at))
 print('PASS',label,'identity from',degree+1,'exact values',flush=True)
 return degree+1
require(42-4==38 and 49-38==11 and 38>=16 and 42>=4,'contact indices')
result=dict(status='PASS',L4_pivot=str(s.factor(pivot)),L3_pivot=str(s.factor(pivot3)),e_Bezout_gcd=str(s.factor(gcd)),scope='model, W, degree11 interpolation, all L4/L3 conditions and pivot coverage; projection separately checked',seconds=time.monotonic()-start,peak_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,certificate_sha256=hashlib.sha256(args.certificate.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
args.out.write_text(json.dumps(result,indent=2)+'\n');print(result,flush=True)
