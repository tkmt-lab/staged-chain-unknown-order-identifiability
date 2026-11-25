from independent import *
import numpy as np
from scipy.optimize import linear_sum_assignment
from sympy import primitive_root
p=104729;N=p-1;g=int(primitive_root(p));gi=pow(g,-1,p)
# Geometric multipoint evaluation via chirp convolution, not interpolation.
b=[1];ib=[1];ratio=iratio=1
for i in range(2*N+1):
 b.append(b[-1]*ratio%p);ib.append(ib[-1]*iratio%p);ratio=ratio*g%p;iratio=iratio*gi%p
points=np.array([pow(g,i,p) for i in range(N)],dtype=np.int64)
assert len(set(map(int,points)))==N

def evaluate(f):
 d=f.degree()
 if d<0:return np.zeros(N,dtype=np.int64)
 aa=nmod_poly([int(f[j])*ib[j]%p for j in range(d,-1,-1)],p)
 conv=aa*nmod_poly(b[:N+d],p)
 vals=np.fromiter((int(conv[d+i])*ib[i]%p for i in range(N)),dtype=np.int64,count=N)
 for i in [0,1,17,N-1]:assert vals[i]==int(f(int(points[i])))
 return vals

def table(rows):
 deg=max(e[1] for e,c in rows);A=np.empty((deg+1,N),dtype=np.int64)
 for j in range(deg+1):
  dic={e[0]:int(c)%p for e,c in rows if e[1]==j}
  f=nmod_poly([dic.get(i,0) for i in range(max(dic,default=-1)+1)],p)
  A[j]=evaluate(f)
 return A.T.copy()
def support_bound(fr,gr,which):
 m=max(e[1] for e,c in fr);n=max(e[1] for e,c in gr)
 M=np.full((m+n,m+n),10**8,dtype=np.int64)
 for start,count,rows,degree in [(0,n,fr,m),(n,m,gr,n)]:
  weights=[]
  for j in range(degree+1):
   hs=[e[0] for e,c in rows if e[1]==j]
   weights.append(min(hs) if which=='valuation' else -max(hs))
  for i in range(count):M[start+i,i:i+degree+1]=weights
 ri,ci=linear_sum_assignment(M);answer=int(M[ri,ci].sum())
 return answer if which=='valuation' else -answer

def main(branch):
 grows=data(f'audit/G{branch}.json.gz')['coefficients'];G=table(grows);print('G evaluated',branch,time.monotonic()-START,flush=True)
 results=[];records=[];hh=nmod_poly([0,1],p);w={82:479,99:500}[branch]
 for jj in [2,3]:
  rows=data(f'audit/CRT_R0{jj}.json.gz')['coefficients'];bound=support_bound(grows,rows,'degree');v=support_bound(grows,rows,'valuation')
  path=f'sources/b14-pilot-2026-11-09/upper/104729/interp_branch{branch}_H01_0{jj}_quotient.json.gz'
  obj=json.loads(gzip.decompress(read(path)));dic={e[0]:int(c) for e,c in obj['coefficients']};candidate=nmod_poly([dic.get(i,0) for i in range(max(dic)+1)],p)*hh**v*(hh+1)**w
  assert candidate.degree()<=bound
  expected=evaluate(candidate);R=table(rows);print('R evaluated',branch,jj,time.monotonic()-START,flush=True)
  count=0;scalar=None
  for i in range(N):
   if not G[i,-1] or not R[i,-1]:continue
   actual=int(nmod_poly(G[i].tolist(),p).resultant(nmod_poly(R[i].tolist(),p)))
   y=int(expected[i])
   if scalar is None and y:scalar=actual*pow(y,-1,p)%p;assert scalar
   if scalar is not None:assert actual==scalar*y%p,(branch,jj,i)
   else:assert actual==0
   count+=1
  assert count>bound and scalar
  results.append(candidate);records.append(dict(pair=f'0{jj}',degree_bound=bound,candidate_degree=candidate.degree(),matching_points=count,scalar=scalar,h_valuation_bound=v))
  print('all evaluations matched',branch,jj,count,time.monotonic()-START,flush=True)
  del R
 U,V=results;assert U.degree()==records[0]['degree_bound']
 D=hh**records[0]['h_valuation_bound']*(hh+1)**w
 if branch==82:
  for f in data('projection/H50.json')['factors']:D*=nmod_poly(list(map(int,f['coefficients'])),p)**(f['multiplicity']+1)
 Q,rem=divmod(U,D);assert not rem and Q.gcd(V).degree()==0 and Q(0) and Q(p-1)
 save(f'PROJECTION_{branch}.json',dict(status='PASS',method='deterministic full finite-field evaluation at more than degree-bound points using geometric chirp convolution',prime=p,records=records,degree_upper=records[0]['degree_bound'],degree_actual=U.degree(),h_valuation_bound=records[0]['h_valuation_bound'],hplus1_divisor=w,U_quotient_degree=Q.degree(),gcd_degree=0))
if __name__=='__main__':main(int(sys.argv[1]))
