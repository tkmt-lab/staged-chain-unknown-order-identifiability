"""Independent exact sparse-rational verification of a polynomial identity and jets."""
import json,math,random,hashlib,time
from fractions import Fraction as F
from pathlib import Path
P=Path(__file__).resolve().parent
N=8;Z=(0,)*N

def clean(a):return {k:v for k,v in a.items() if v}
def add(a,b):
 c=a.copy()
 for k,v in b.items():c[k]=c.get(k,F(0))+v
 return clean(c)
def mul(a,b):
 c={}
 for x,u in a.items():
  for y,v in b.items():
   z=tuple(i+j for i,j in zip(x,y));c[z]=c.get(z,F(0))+u*v
 return clean(c)
def neg(a):return {k:-v for k,v in a.items()}
def var(i):
 e=[0]*N;e[i]=1;return {tuple(e):F(1)}
def diff(a,i):
 c={}
 for e,v in a.items():
  if e[i]:
   k=list(e);k[i]-=1;c[tuple(k)]=v*e[i]
 return c
def decode(data):return {tuple(m):F(n,d) for m,n,d in data}
b,k,s,x,B,K,S,X=map(var,range(N))
zero={};left=mul(mul(b,s),x);right=mul(mul(B,S),X)
field=[zero,zero,neg(left),add(mul(k,x),left),zero,zero,neg(right),add(mul(K,X),right)]
def delta(a):
 z={}
 for i,v in enumerate(field):z=add(z,mul(diff(a,i),v))
 return z
p=[add(left,neg(right))]
for _ in range(5):p.append(delta(p[-1]))
c=json.loads((P/'CERTIFICATE.json').read_text());coeff=list(map(decode,c['certificate']))
total={}
for a,g in zip(coeff,p[:5]):total=add(total,mul(a,g))
assert total==p[5], 'coefficient identity failed'
def jet(beta,k,s0,x0,n=5):
 S=[F(s0)];X=[F(x0)];Y=[]
 for i in range(n+1):
  yi=F(beta)*sum(F(math.comb(i,j))*S[j]*X[i-j] for j in range(i+1));Y.append(yi)
  S.append(-yi);X.append(F(k)*X[i]+yi)
 return Y
def evalpoly(a,vals):return sum(v*math.prod(z**e for z,e in zip(vals,key)) for key,v in a.items())
L=(F(1),F(-1),F(1),F(2));R=(F(3,2),F(-4),F(2),F(2,3));yl,yr=jet(*L),jet(*R)
assert yl[:4]==yr[:4] and yl[4]!=yr[4]
assert sum([L[2],L[3],F(1)])==sum([R[2],R[3],F(4,3)])==4
assert all(t>0 for t in [L[0],L[2],L[3],F(1),R[0],R[2],R[3],F(4,3)])
random.seed(4405)
for _ in range(40):
 a=[F(random.randint(1,5),random.randint(1,4)),F(-random.randint(1,5),random.randint(1,4)),F(random.randint(1,5),random.randint(1,4)),F(random.randint(1,5),random.randint(1,4))]
 b2=[F(random.randint(1,5),random.randint(1,4)),F(-random.randint(1,5),random.randint(1,4)),F(random.randint(1,5),random.randint(1,4)),F(random.randint(1,5),random.randint(1,4))]
 j1,j2=jet(*a),jet(*b2)
 for i in range(6):assert evalpoly(p[i],a+b2)==j1[i]-j2[i]
result={'status':'PASS','polynomial_identity_terms':len(total),'certificate_coefficient_terms':[len(v) for v in coeff],'certificate_sha256':hashlib.sha256((P/'CERTIFICATE.json').read_bytes()).hexdigest(),'left_jet':[str(v) for v in yl],'right_jet':[str(v) for v in yr],'positive_initial_components':True,'common_population':4,'rational_pair_checks':40,'checks_per_pair':6,'sympy_used_in_verifier':False}
(P/'VERIFY.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
