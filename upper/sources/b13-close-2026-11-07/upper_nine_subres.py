"""Native multivariate coefficients, Brown's exact subresultant PRS in t."""
import resource,sys
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));sys.set_int_max_str_digits(0)
import flint,json,time
from pathlib import Path
p=Path(__file__).parent;start=time.monotonic();C=flint.fmpz_mpoly_ctx.get(('h','k','t'));h,k,t=C.gens();J=json.loads((p/'lower_projection.json').read_text());E=[flint.fmpz_mpoly(s,C) for s in J['polynomials']]
def degree(F):return int(F.degrees()[2])
def lc(F):
 d=degree(F);return C.from_dict({(hh,kk,0):v for (hh,kk,tt),v in F.to_dict().items() if tt==d})
def prem(A,B):
 n=degree(A);m=degree(B);left=n-m+1;L=lc(B)
 while A and degree(A)>=m:
  d=degree(A)-m;A=L*A-lc(A)*t**d*B;left-=1
 assert left>=0
 return L**left*A

def exact(A,B):
 q,r=divmod(A,B);assert r==0;return q
f,g=E[:2];n,m=degree(f),degree(g);R=[f,g];d=n-m;b=C.constant((-1)**(d+1));z=prem(f,g)*b;L=lc(g);cc=L**d;scalar=[C.constant(1),cc];cc=-cc
while z:
 dd=degree(z);R.append(z);print('SR',dd,z.degrees(),len(z.to_dict()),time.monotonic()-start,flush=True)
 out={'variables':['h','k','t'],'degrees':[degree(x) for x in R],'subresultants':[str(x) for x in R],'seconds':time.monotonic()-start};(p/'upper_nine_subres.json').write_text(json.dumps(out,indent=2)+'\n')
 f,g,m,d=g,z,dd,m-dd;b=-L*cc**d;z=exact(prem(f,g),b);L=lc(g)
 cc=exact((-L)**d,cc**(d-1)) if d>1 else -L;scalar.append(-cc)
expected=flint.fmpz_mpoly(J['resultants']['01'],C);assert R[-1]==expected or R[-1]==-expected;print('constant resultant matches native PASS',time.monotonic()-start,flush=True)
