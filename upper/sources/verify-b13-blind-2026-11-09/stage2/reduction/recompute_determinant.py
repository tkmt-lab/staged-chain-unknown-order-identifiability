"""Regenerate all 70 contact entries and determinants A,P6 from model polynomials."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import sympy as s,flint
from flint import fmpz_mpoly_ctx,fmpz_mpoly
from pathlib import Path
import json,time,resource,signal,hashlib,datetime
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));signal.alarm(1800)
t0=time.monotonic();here=Path(__file__).resolve().parent
S,X,h,k,r=s.symbols('S X h k r');q=S+h;y=S*X;ell=S+k-X
D=lambda a:s.expand(-S*X*s.diff(a,S)+(S+k)*X*s.diff(a,X))
Dx=lambda a:s.expand(S*s.diff(a,S)-(S+k)*s.diff(a,X))
N=[s.Integer(1)];M=[s.Integer(1)]
for n in range(3):N.append(s.expand(q*(ell*N[-1]+D(N[-1]))+(n+1)*y*N[-1]))
for n in range(2):M.append(s.expand(ell*M[-1]+D(M[-1])))
B=[q*N[2],q**2*N[1],q**3,-q**4*M[2],-q**4*M[1],-q**4,N[3]]
C=fmpz_mpoly_ctx.get(['h','k','r']);zero=C.constant(0)
def convert(p):return C.from_dict({mon:int(c) for mon,c in s.Poly(p,h,k,r).terms()})
rows=[]
for j in range(10):
 rows.append([convert(p.subs({S:1,X:r}).expand()) for p in B]);B=[Dx(p) for p in B]
print('generated rows',time.monotonic()-t0,flush=True)
p=Path('sources/b13-lower3-2026-11-07/elim_general.json');raw=p.read_bytes()
with (here/'access.jsonl').open('a') as f:f.write(json.dumps({'path':str(p),'sha256':hashlib.sha256(raw).hexdigest(),'opened_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'mode':'full 70 entries and A/P6 computation comparison'})+'\n')
old=json.loads(raw)
for i,row in enumerate(rows):
 for j,v in enumerate(row):assert v==fmpz_mpoly(old['rows'][i][j].replace('**','^'),C),(i,j)
print('all 70 rows equal',time.monotonic()-t0,flush=True)
def bareiss(mat):
 a=[row[:] for row in mat];last=C.constant(1);sgn=1;n=len(a)
 for j in range(n-1):
  if not a[j][j]:
   p=next(i for i in range(j+1,n) if a[i][j]);a[p],a[j]=a[j],a[p];sgn=-sgn
  pivot=a[j][j]
  for i in range(j+1,n):
   for l in range(j+1,n):
    numer=pivot*a[i][l]-a[i][j]*a[j][l]
    a[i][l]=numer//last
    assert a[i][l]*last==numer
   a[i][j]=zero
  last=pivot
  print('Bareiss',n,'step',j,'sec',time.monotonic()-t0,flush=True)
 return sgn*a[-1][-1]
checks={}
for name,n in [('A',6),('P6',7)]:
 actual=bareiss([row[:n] for row in rows[:n]])
 expected=fmpz_mpoly(old[name].replace('**','^'),C)
 assert actual==expected,name
 checks[name]={'exact_equal':True,'terms':len(actual.to_dict()),'degrees':list(map(int,actual.degrees())),'normalized_sha256':hashlib.sha256(str(actual).encode()).hexdigest()}
 print(name,'PASS',time.monotonic()-t0,flush=True)
# Formal exterior-power recurrence for derivatives of det(R0,...,R6).
from collections import defaultdict
terms={tuple(range(7)):1};expansions=[]
for order in range(4):
 expansions.append({'order':order,'terms':[[list(mon),coeff] for mon,coeff in sorted(terms.items())]})
 nxt=defaultdict(int)
 for mon,coeff in terms.items():
  for j in range(7):
   new=list(mon);new[j]+=1
   if len(set(new))<7:continue
   inversions=sum(new[i]>new[l] for i in range(7) for l in range(i+1,7));nxt[tuple(sorted(new))]+=coeff*(-1)**inversions
 terms={mon:coeff for mon,coeff in nxt.items() if coeff}
assert expansions[3]['terms']==[[[0,1,2,3,4,5,9],1],[[0,1,2,3,4,6,8],2],[[0,1,2,3,5,6,7],1]]
out={'status':'PASS','method':'original ODE numerator recurrence; exact polynomial Bareiss with every division checked by multiplication','flint':flint.__version__,'source_sha256':hashlib.sha256(raw).hexdigest(),'all_70_contact_entries_equal':True,'determinants':checks,'formal_determinant_derivatives':expansions,'seconds':time.monotonic()-t0,'maxrss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(here/'recompute_determinant.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
