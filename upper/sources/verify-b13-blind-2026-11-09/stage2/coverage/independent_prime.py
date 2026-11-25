"""A fresh prime and implementation of all projected coprimality conclusions.
Exact divisibility in characteristic zero is supplied by replayed certificates;
this file verifies degree preservation and quotient coprimality independently.
"""
from pathlib import Path
from fractions import Fraction
from math import isqrt
import flint,json,resource,signal,time,sys
sys.set_int_max_str_digits(0)
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,)*2);signal.alarm(1750)
O=Path(__file__).resolve().parent;P=O.parents[2]/'b13-close-2026-11-07';start=time.monotonic();p=65537
assert all(p%d for d in range(2,isqrt(p)+1))
read=lambda n:json.loads((P/n).read_text())
raw=read('lower_projection.json');new=read('upper_contact10.json');lin=read('upper_nine_linear.json');rur=read('lower_RUR_CANDIDATE.json');bound=read('lower_degree_certificate.json')
Z=flint.fmpz_mpoly_ctx.get(('h','k'));C=flint.nmod_mpoly_ctx.get(('h','k'),p);h,k=C.gens()
def parse(expr):
 z=flint.fmpz_mpoly(expr.replace('**','^'),Z);v=C.from_dict({m:int(c)%p for m,c in z.to_dict().items()});assert z.degrees()[1]==v.degrees()[1];return v
Gs={str(n):parse(next(f['polynomial'] for f in raw['factors']['01']['factors'] if f['degrees'][1]==n)) for n in [22,27]}
H02=parse(next(f['polynomial'] for f in raw['factors']['02']['factors'] if f['degrees'][1]>1));H03=next(parse(expr) for expr,e in new['factors']['03'] if flint.fmpz_mpoly(expr,Z).degrees()[1]>1)
R={};report={'prime':p,'resultant_degrees':{},'B12':{}}
for n,hh,qq in [('22',62,92),('27',90,98)]:
 D=h**hh*(h+1)**qq
 if n=='22':D*=(3*h+5)**2*(3*h*h+8*h-4)*(14*h*h+40*h+29)
 qs=[]
 for key,F in [('02',H02),('03',H03)]:
  a=Gs[n].resultant(F,'k');R[(n,key)]=a;assert a.degrees()[0]==bound['branches'][n][key]['degree']['bound'];q,r=divmod(a,D);assert not r;qs.append(q);report['resultant_degrees'][n+':'+key]=int(a.degrees()[0])
 assert qs[0].gcd(qs[1]).is_constant();report['B12'][n]='quotient gcd 1'
HC=list(map(Fraction,rur['H_coefficients_low_first']));assert HC[-1]==1;assert all(c.denominator%p for c in HC)
H=sum((c.numerator%p)*pow(c.denominator%p,-1,p)*h**i for i,c in enumerate(HC));assert H.degrees()[0]==289 and H.gcd(H.derivative('h')).is_constant()
G=Gs['27'];Euler=h*G.derivative('h')+k*G.derivative('k');S=G.resultant(Euler,'k');D=h**90*(h+1)**98*H;q,r=divmod(R[('27','02')],D);assert not r and q.gcd(S).is_constant();RH,rem=divmod(R[('27','02')],H);assert not rem and H.gcd(RH).is_constant();LC=C.from_dict({(i,0):c for (i,j),c in G.to_dict().items() if j==27});assert H.gcd(LC).is_constant()
A=parse(lin['a']);Ra=G.resultant(A,'k');U=h**90*(h+1)**98;Q1,e1=divmod(Ra,U);Q2,e2=divmod(R[('27','02')],U);assert not e1 and not e2 and Q1.gcd(Q2).is_constant()
report.update(status='PASS',RUR_coverage='quotient gcd 1',H_occurs_once=True,G27_degree_preserved_over_H=True,linear_a_nonzero_on_G27=True,seconds=time.monotonic()-start)
(O/'independent_prime.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
