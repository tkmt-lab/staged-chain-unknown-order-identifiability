"""Exact normalization, differential recurrence, and factorization checks."""
from pathlib import Path
import json,hashlib,sys,resource,signal
from flint import fmpz_mpoly_ctx,fmpz_mpoly
sys.set_int_max_str_digits(0);resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,)*2);signal.alarm(1750)
O=Path(__file__).resolve().parent;P=O.parent/'b13-lower3-2026-11-07';raw=json.loads((P/'elim_general.json').read_text());data=json.loads((O/'lower_projection.json').read_text());C=fmpz_mpoly_ctx.get(['h','k','t']);h,k,t=C.gens();q=1+h
p6=fmpz_mpoly(raw['P6'].replace('**','^').replace('r','t'),C);f=p6//(-32*h**3);assert p6==(-32*h**3)*f
ps=[fmpz_mpoly(v,C) for v in data['polynomials']]
# q*f(h,k,r) = c*sum(a_abc h^a k^b r^c q^(7-c)).
mapped=sum((int(v)*h**a*k**b*t**c*q**(7-c) for (a,b,c),v in ps[0].to_dict().items()),C.constant(0))
def proportional(a,b):
 da=a.to_dict();db=b.to_dict();mon=next(iter(db));assert da[mon]%db[mon]==0;c=da[mon]//db[mon];assert a==int(c)*b;return str(c)
constants={'normalization':proportional(q*f,mapped)}
W=lambda p:-h*q*p.derivative('h')-k*q*p.derivative('k')-(t+1+k)*p.derivative('t')
for j in range(2):constants['W'+str(j)+'_to_'+str(j+1)]=proportional(W(ps[j]),ps[j+1])
factor_checks={}
for key,pair in [('01',(0,1)),('02',(0,2)),('12',(1,2))]:
 expected=fmpz_mpoly(data['resultants'][key],C);actual=ps[pair[0]].resultant(ps[pair[1]],'t');assert actual==expected
 record=data['factors'][key];product=C.constant(int(record['constant']))
 for item in record['factors']:product*=fmpz_mpoly(item['polynomial'],C)**item['power']
 assert product==expected;factor_checks[key]={'resultant_recomputed_equal':True,'factor_product_equal':True,'degree':list(map(int,actual.degrees()))}
out={'status':'PASS','normalization_constants':constants,'checks':factor_checks,'inputs':{str(p.relative_to(O.parent.parent.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [P/'elim_general.json',O/'lower_projection.json']}}
(O/'lower_PROJECTION_CHECK.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS',constants)
