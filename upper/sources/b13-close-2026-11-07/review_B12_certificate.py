"""Short audit of the degree-lift prerequisites; no large resultant computation."""
from pathlib import Path
import json,time
from math import isqrt
from flint import fmpz_mpoly_ctx,fmpz_mpoly,nmod_mpoly_ctx,nmod_mpoly,fmpq
P=Path(__file__).parent;start=time.monotonic();prime=32003
assert all(prime%d for d in range(2,isqrt(prime)+1))
Z=fmpz_mpoly_ctx.get(['h','k']);M=nmod_mpoly_ctx.get(['h','k'],prime)
raw=json.loads((P/'lower_projection.json').read_text());new=json.loads((P/'upper_contact10.json').read_text());cert=json.loads((P/'lower_degree_certificate.json').read_text());out={'status':'PASS','prime':prime,'prime_trial_divisor_limit':isqrt(prime),'k_degree_preservation':{},'assignment_certificates':{}}
def parse(x):return fmpz_mpoly(x.replace('**','^'),Z)
def weights(F):
 d=F.to_dict();m=int(F.degrees()[1]);co=[{} for _ in range(m+1)]
 for (i,j),v in d.items():co[m-int(j)][int(i)]=v
 return [None if not row else {'degree':max(row),'h':min(row)} for row in co]
for tag in ['27','22']:
 G=parse(next(f['polynomial'] for f in raw['factors']['01']['factors'] if f['degrees'][1]==int(tag)))
 for key in ['02','03']:
  H=parse(next(f['polynomial'] for f in raw['factors']['02']['factors'] if f['degrees'][1]>1)) if key=='02' else next(parse(f) for f,m in new['factors']['03'] if parse(f).degrees()[1]>1)
  degs=[int(F.degrees()[1]) for F in [G,H]];moddegs=[int(nmod_mpoly(str(F),M).degrees()[1]) for F in [G,H]];assert degs==moddegs
  label=tag+':'+key;out['k_degree_preservation'][label]=degs
  a,b=weights(G),weights(H);m,n=degs;mat=[[None]*i+a+[None]*(n-1-i) for i in range(n)]+[[None]*i+b+[None]*(m-1-i) for i in range(m)]
  checked={}
  for name,sg in [('degree',-1),('h',1)]:
   C=cert['branches'][tag][key][name];ass=C['assignment'];u=C['row_dual'];v=C['column_dual'];N=m+n
   assert sorted(ass)==list(range(N)) and len(u)==len(v)==N
   assert all(mat[i][ass[i]] is not None for i in range(N))
   # Only finite entries are relevant to determinant monomials.
   assert all(u[i]+v[j]<=sg*mat[i][j][name] for i in range(N) for j in range(N) if mat[i][j] is not None)
   value=sum(sg*mat[i][ass[i]][name] for i in range(N));assert value==sum(u)+sum(v)==C['value'];assert sg*value==C['bound']
   checked[name]=C['bound']
  assert checked['degree']==cert['branches'][tag][key]['mod_degree'];out['assignment_certificates'][label]=checked
# Verify the algebraic condition for second-order vanishing from stored exact derivatives.
T=json.loads((P/'contact10_tangent.json').read_text());Gh=fmpq(T['G_h']);Gk=fmpq(T['G_k']);assert Gk
for row in T['pairs'].values():
 Rh=fmpq(row['H_h']);Rk=fmpq(row['H_k']);assert Rk and Gh*Rk-Gk*Rh==0
out['tangent_derivative_identity']=True;out['seconds']=time.monotonic()-start
(P/'review_B12_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS short B12 lift-prerequisite review')
