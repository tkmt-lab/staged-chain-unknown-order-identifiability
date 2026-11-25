"""G22 generic common root has A=0, checked by homogeneous substitution."""
import resource,signal,time,json,hashlib
from pathlib import Path
from flint import fmpz_mpoly_ctx,fmpz_mpoly
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));signal.alarm(1750)
P=Path(__file__).parent;start=time.monotonic();C=fmpz_mpoly_ctx.get(('h','k'));C3=fmpz_mpoly_ctx.get(('h','k','t'));conv=lambda f:fmpz_mpoly(f.replace('**','^'),C);raw=json.loads((P/'lower_projection.json').read_text());lin=json.loads((P/'upper_nine_linear.json').read_text());G=conv(next(f['polynomial'] for f in raw['factors']['01']['factors'] if f['degrees'][1]==22));a,b=map(conv,[lin['a'],lin['b']]);A=fmpz_mpoly(json.loads((P/'structure_nine_pivots.json').read_text())['transformed']['A'],C3)
gcd=a.gcd(G);assert gcd==1;print('gcd(a,G22)=1',time.monotonic()-start,flush=True)
d=A.degrees()[2];cs=[C.from_dict({(i,j):co for (i,j,n),co in A.to_dict().items() if n==v}) for v in range(d+1)]
val=C.constant(0)
for j,c in enumerate(cs):val+=c*(-b)**j*a**(d-j);print('added',j,'terms',len(val),'sec',time.monotonic()-start,flush=True)
quot,rem=divmod(val,G);assert not rem
print('PASS A(-b/a)*a^4 divisible by G22',quot.degrees(),len(quot),time.monotonic()-start,flush=True)
out={'status':'PASS','A_t_degree':int(d),'gcd_a_G22':'1','homogeneous_A_divisible_by_G22':True,'homogeneous_A_degree':list(map(int,val.degrees())),'quotient_degree':list(map(int,quot.degrees())),'quotient_terms':len(quot),'quotient_sha256':hashlib.sha256(str(quot).encode()).hexdigest(),'G22_multiplicity_in_F_WF_resultant':next(f['power'] for f in raw['factors']['01']['factors'] if f['degrees'][1]==22),'seconds':time.monotonic()-start}
(P/'structure_G22_A.json').write_text(json.dumps(out,indent=2)+'\n')
