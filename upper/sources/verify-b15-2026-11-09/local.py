from independent import *
def qp(a):return fmpq_poly([fmpq(x) for x in a])
def trunc(a,N):return a.truncate(N)
def mul(a,b,N):return (a*b).truncate(N)
def inv(a,N):
 assert a[0]
 # Newton inverse in Q[q]/q^N
 b=qp([1/a[0]]);n=1
 while n<N:
  n=min(2*n,N);b=mul(b,2-mul(a,b,n),n)
 assert mul(a,b,N)==1
 return b
def read_coeff(p,N,shift=False):
 rows=data(p)['coefficients'];d=max(e[1] for e,c in rows);out=[]
 for j in range(d+1):
  dic={e[0]:fmpq(c) for e,c in rows if e[1]==j}
  a=qp([dic.get(i,0) for i in range(max(dic,default=-1)+1)])
  if shift:
   a=qp([sum(a[j]*math.comb(j,i)*(-1)**(j-i) for j in range(i,a.degree()+1)) for i in range(N)])
  out.append(trunc(a,N))
 return out
def mobius(a,c,N):
 # w^degree a(c+1/w); an integral GL2 transformation
 m=len(a)-1;out=[qp([]) for _ in a]
 for i,v in enumerate(a):
  for j in range(i+1):out[m-j]+=v*(math.comb(i,j)*c**(i-j))
 return [trunc(v,N) for v in out]
def rem(a,b,N):
 a=list(a);ib=inv(b[-1],N)
 for j in range(len(a)-len(b),-1,-1):
  c=mul(a[j+len(b)-1],ib,N)
  for i in range(len(b)):a[j+i]-=mul(c,b[i],N)
 return [trunc(v,N) for v in a[:len(b)-1]]
def valuation(a,N):return min((i for v in a for i,c in enumerate(v.coeffs()) if c),default=N)
def bound(a,b,N,target):
 total=0;trace=[]
 while total<target:
  m,n=len(a)-1,len(b)-1
  va,vb=valuation(a,N),valuation(b,N)
  total+=n*va+m*vb;trace.append([m,n,N,va,vb,total])
  if total>=target:return total,trace
  assert va<N and vb<N
  N-=max(va,vb)
  a=[trunc(v.right_shift(va),N) for v in a];b=[trunc(v.right_shift(vb),N) for v in b]
  if n==0:return total,trace
  if m<n:a,b=b,a;m,n=n,m
  if not a[-1][0] or not b[-1][0]:
   c=0
   while not sum(v[0]*c**i for i,v in enumerate(a)) or not sum(v[0]*c**i for i,v in enumerate(b)):c+=1
   a,b=mobius(a,c,N),mobius(b,c,N)
  a,b=b,rem(a,b,N)
 return total,trace

def main(branch):
 setup=data(f'projection/CLUSTER_G{branch}.json'); cert=data(f'projection/CLUSTER_CERT_G{branch}.json')
 results=[]
 for item in cert['cluster_values']:
  i=item['index'];N=item['precision'];target=item['lower_bound']
  F=read_coeff(f'projection/cluster_G{branch}_{i}_F.json',N)
  savedG=read_coeff(f'projection/cluster_G{branch}_{i}_G.json',N)
  # Saved cluster may have rational denominator cleared. Normalize by its unit LC.
  assert F[-1][0];unit=inv(F[-1],N);F=[mul(a,unit,N) for a in F]
  q0=qp(setup['clusters'][i]['factor_coefficients'])**setup['clusters'][i]['multiplicity']
  assert qp([a[0] for a in F])==q0
  orig=mobius(read_coeff(f'audit/G{branch}.json.gz',N,True),setup['mobius_c'],N)
  assert orig[-1][0]
  assert not any(rem(orig,F,N)),('factor does not divide',branch,i)
  # Constant complementary factor coprime guarantees unique Hensel branch.
  orig0=qp([a[0] for a in orig]);comp,r=divmod(orig0,q0);assert not r and comp.gcd(q0).degree()==0
  R=mobius(read_coeff('audit/CRT_R02.json.gz',N,True),setup['mobius_c'],N)
  G=rem(R,F,N)
  v,trace=bound(F,G,N,target);assert v>=target
  # Also independently audit the stored local pair; no trust in status fields.
  vs,ts=bound(F,savedG,N,target);assert vs>=target
  results.append(dict(index=i,precision=N,lower_bound=v,trace=trace,factor_divides=True,complement_coprime=True))
  print('cluster',branch,i,v,'elapsed',time.monotonic()-START,flush=True)
  save(f'LOCAL_{branch}.json',dict(status='RUNNING',clusters=results))
 assert sum(r['lower_bound'] for r in results)>=cert['required']
 save(f'LOCAL_{branch}.json',dict(status='PASS',clusters=results,lower_bound=sum(r['lower_bound'] for r in results)))
if __name__=='__main__':main(int(sys.argv[1]))
