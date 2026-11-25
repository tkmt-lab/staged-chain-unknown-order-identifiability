from independent import *
from fractions import Fraction
import sympy as sp

def qp(a):return fmpq_poly([fmpq(x) for x in a])
def strip(a):
 while a and not a[-1]:a.pop()
 return a
def gcd_ext(a,b,q):
 while b:
  inv=b[-1].xgcd(q)[1]%q
  a=list(a)
  while len(a)>=len(b):
   c=a[-1]*inv%q;j=len(a)-len(b)
   for i in range(len(b)):a[j+i]=(a[j+i]-c*b[i])%q
   strip(a)
  a,b=b,a
 inv=a[-1].xgcd(q)[1]%q
 return [(c*inv)%q for c in a]
def eval_hk(rows,q,K):
 # Evaluate using nested Horner in k, with polynomial coefficients in h.
 deg=max(e[-1] for e,c in rows) if len(rows[0][0])==3 else 0
 out=[]
 for z in range(deg+1):
  grouped={}
  for e,c in rows:
   if len(e)==3 and e[2]!=z:continue
   a=grouped.setdefault(e[1],{})
   a[e[0]]=int(c)
  val=qp([])
  for j in range(max(grouped,default=-1),-1,-1):
   d=grouped.get(j,{});coeff=qp([d.get(a,0) for a in range(max(d,default=-1)+1)])%q
   val=(val*K+coeff)%q
  out.append(val)
 return strip(out)
def reduce_h(rows,q,constructor):
 # polynomial in k over finite extension
 out=[]
 for j in range(max(e[1] for e,c in rows)+1):
  d={e[0]:int(c) for e,c in rows if e[1]==j}
  out.append(constructor([d.get(i,0) for i in range(max(d,default=-1)+1)])%q)
 return strip(out)
H=data('projection/H50.json');factors=[qp(x['coefficients']) for x in H['factors']]
Frows=[data('exact/'+s+'.json')['coefficients'] for s in ['F4','W1','W2','W3']]
Rrows=[data(f'audit/CRT_R0{j}.json.gz')['coefficients'] for j in [1,2,3]]
Grows=data('audit/G82.json.gz')['coefficients']
records=[]
x=sp.Symbol('h'); rootdata=data('projection/REAL_ROOTS.json')['roots'];allboxes=[]
p=104729;constructor=lambda a:nmod_poly([int(c)%p for c in a],p)
for index,q in enumerate(factors):
 st=time.monotonic(); fi=data(f'projection/QQ_FIBER_{index}.json');K=qp(fi['k_value_coefficients'])
 spq=sp.Poly.from_list(list(reversed([int(c) for c in q.coeffs()])),x);assert spq.is_irreducible
 fs=[eval_hk(a,q,K) for a in Frows];assert all(not a[0] for a in fs)
 g=fs[0]
 for f in fs[1:]:g=gcd_ext(g,f,q)
 assert len(g)==2 and not g[0] and g[1]==1
 Q=constructor(q.coeffs());assert Q.gcd(Q.derivative()).degree()==0
 bridge=[]
 for v,m in Q.factor()[1]:
  assert m==1
  rs=[reduce_h(a,v,constructor) for a in Rrows];assert len(rs[0])==182
  gg=rs[0]
  for a in rs[1:]:gg=gcd_ext(gg,a,v)
  km=constructor([int(c.numerator)*pow(int(c.denominator),-1,p)%p for c in K.coeffs()])%v
  assert len(gg)==2 and gg[0]==-km and gg[1]==1
  bridge.append(v.degree())
 qa=data(f'projection/QADIC_R02_H80_{index}.json');qm=q**qa['multiplicity'];section=qp(qa['section_coefficients'])
 assert not eval_hk(Grows,qm,section) and not eval_hk(Rrows[1],qm,section)
 boxes=[b for b in rootdata if b['factor_index']==index]
 assert spq.count_roots(-sp.oo,sp.oo)==len(boxes)
 for b in boxes:
  lo,hi=map(sp.Rational,b['h_interval']);assert lo<=hi
  assert (spq.eval(lo)==0 if lo==hi else spq.count_roots(lo,hi)==1)
  assert bool(lo>-1)==b['h_gt_minus1'];assert (lo>-1 or hi<-1)
  # interval Horner, exact rationals
  al,ah=Fraction(str(lo)),Fraction(str(hi));l=u=Fraction(0)
  for c in reversed(K.coeffs()):
   prod=[l*al,l*ah,u*al,u*ah];l=min(prod)+Fraction(str(c));u=max(prod)+Fraction(str(c))
  assert u<0
  allboxes.append(dict(factor=index,h=[str(lo),str(hi)],k_interval=[str(l),str(u)],h_gt_minus1=bool(lo>-1)))
 records.append(dict(factor=index,degree=q.degree(),t_gcd='t',modular_factors=bridge,qadic_power=qa['multiplicity'],roots=len(boxes)))
 print('fiber',index,records[-1],'time',time.monotonic()-st,flush=True)
 save('FIBERS.json',dict(status='RUNNING',factors=records,boxes=allboxes))
assert len(allboxes)==16 and sum(b['h_gt_minus1'] for b in allboxes)==10
save('FIBERS.json',dict(status='PASS',factors=records,boxes=allboxes,real_roots=16,admissible_h=10))
