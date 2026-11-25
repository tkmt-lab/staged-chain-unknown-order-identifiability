"""Verify integer resultant identities by fixed-size Sylvester determinants.
Every equality is certified by a degree bound and D+1 distinct exact values.
The final modular computation proves coprimality of already verified ZZ polynomials.
"""
import sys,os,json,time,resource,signal,argparse,hashlib
from pathlib import Path
O=Path(__file__).resolve().parent;sys.path.insert(0,str(O/'.deps'))
import sympy as s
from flint import fmpz_mat,fmpz_poly
ap=argparse.ArgumentParser();ap.add_argument('--certificate',type=Path,default=O/'certificate42.json');ap.add_argument('--projection',type=Path,default=O/'projection42.json');ap.add_argument('--alternate',type=Path,default=O/'alternate42.json');ap.add_argument('--out',type=Path,default=O/'PROJECTION_CHECK42.json');args=ap.parse_args()
signal.alarm(3600);resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));start=time.monotonic()
H,u,A=s.symbols('H u A')
def require(v,msg):
 if not v:raise ValueError('REJECT: '+msg)
def dec(z,vs):
 require(len({tuple(e) for e,v in z})==len(z),'duplicate terms')
 require(all(len(e)==len(vs) and all(isinstance(i,int) and i>=0 for i in e) for e,v in z),'exponents')
 return sum(s.Rational(v)*s.prod(t**i for t,i in zip(vs,e)) for e,v in z)
def strip(z):
 p=s.Poly(z,H,u);powers=[min(e[i] for e in p.monoms()) for i in range(2)]
 return s.Poly(s.cancel(z/H**powers[0]/u**powers[1]),H,u).clear_denoms()[1].primitive()[1].as_expr()
C=json.loads(args.certificate.read_text());cr=C['cross'];proj=json.loads(args.projection.read_text());alt=json.loads(args.alternate.read_text());P=dec(cr['P'],[H,u]);checks=[];Rs=[]
# Three resultants suffice if their gcd is one; no trust in the saved gcd field.
selected=[('base',row) for row in proj['results'][:2]]+[('alternate',row) for row in alt['results']]
for kind,row in selected:
 i=row['i'];label=kind+str(i)
 if kind=='base':Z=dec(cr['Z'][i],[H,u])
 else:
  t=dec(C['elim']['T'][2],[H,u,A]);a=t.coeff(A);b=t.subs(A,0);z=dec(cr['minors'][i],[H,u,A])
  Z=dec(row['Z'],[H,u]);require(s.expand(Z-strip(z.coeff(A,2)*b**2-z.coeff(A)*a*b+z.subs(A,0)*a**2))==0,'alternate cross')
 pp=s.Poly(P,H,u);zz=s.Poly(Z,H,u);n=pp.degree(H);m=zz.degree(H);D=m*pp.degree(u)+n*zz.degree(u)
 r=fmpz_poly([int(v) for v in row['resultant']]);content=int(row['content']);require(content!=0,'content nonzero');require(0<=r.degree()<=D,'resultant degree bound')
 def coeff(poly,N,at):
  ans=[0]*(N+1)
  for (eh,eu),co in poly.terms():require(co.is_Integer,'integer coefficient');ans[N-eh]+=int(co)*at**eu
  return ans
 for j in range(D+1):
  at=0 if j==0 else (j+1)//2*(1 if j%2 else -1)
  a0=coeff(pp,n,at);b0=coeff(zz,m,at)
  M=[[0]*i+a0+[0]*(m-1-i) for i in range(m)]+[[0]*i+b0+[0]*(n-1-i) for i in range(n)]
  require(fmpz_mat(M).det()==content*r(at),label+' exact determinant at '+str(at))
 checks.append(dict(label=label,degree_bound=D,points=D+1,matrix_size=n+m,resultant_degree=r.degree(),seconds=time.monotonic()-start))
 Rs.append([int(c) for c in r.coeffs()]);print('PASS',checks[-1],flush=True)
p=1000003
# Independent ordinary modular polynomial Euclid, ascending coefficient order.
def trim(a):
 while a and a[-1]==0:a.pop()
 return a
def rem(a,b):
 a=a[:];inv=pow(b[-1],-1,p)
 while len(a)>=len(b):
  v=a[-1]*inv%p;shift=len(a)-len(b)
  for j,x in enumerate(b):a[shift+j]=(a[shift+j]-v*x)%p
  trim(a)
 return a
g=[]
for r in Rs:
 rr=trim([c%p for c in r]);require(len(rr)==len(r),'modular degree preserved')
 if not g:g=rr
 else:
  while rr:g,rr=rr,rem(g,rr)
require(len(g)==1,'QQ coprimality')
result=dict(status='PASS',scope='global finite-H exclusion by QQ resultant identities and coprimality',checks=checks,prime=p,modular_degrees_preserved=True,gcd_degree=0,seconds=time.monotonic()-start,rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,inputs={str(f.name):hashlib.sha256(f.read_bytes()).hexdigest() for f in (args.certificate,args.projection,args.alternate)},checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
args.out.write_text(json.dumps(result,indent=2)+'\n');print('PASS global exclusion',flush=True)
