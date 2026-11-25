"""Exact independent substitution of the H2 finite projection branch."""
import resource,signal,time,json
from pathlib import Path
import sympy as s
from flint import fmpz_mpoly_ctx,fmpz_mpoly,fmpq_poly
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));signal.alarm(1750)
P=Path(__file__).parent;start=time.monotonic();ctx=fmpz_mpoly_ctx.get(('h','k','t'));h,k,t=s.symbols('h k t');HP=s.Poly(3*h*h+8*h-4,h);HM=fmpq_poly([-4,8,3]);hh=fmpq_poly([0,1]);kk=fmpq_poly([-14,-3])/2;raw=json.loads((P/'lower_projection.json').read_text());new=json.loads((P/'upper_contact10.json').read_text())
def parse(f):return fmpz_mpoly(f.replace('**','^'),ctx)
def subst(F):
 d=F.to_dict();ph=[fmpq_poly([1])];pk=[fmpq_poly([1])]
 for i in range(F.degrees()[0]):ph.append(ph[-1]*hh%HM)
 for i in range(F.degrees()[1]):pk.append(pk[-1]*kk%HM)
 out={}
 for (i,j,n),co in d.items():out[n]=(out.get(n,fmpq_poly([]))+int(co)*ph[i]*pk[j])%HM
 return out
Gs={'G22':parse(next(f['polynomial'] for f in raw['factors']['01']['factors'] if f['degrees'][1]==22)), 'H02':parse(next(f['polynomial'] for f in raw['factors']['02']['factors'] if f['degrees'][1]>1)), 'H03':next(parse(f) for f,m in new['factors']['03'] if parse(f).degrees()[1]>1)}
for name,F in Gs.items():assert all(not cc for cc in subst(F).values()),name;print('projection substitution zero',name,flush=True)
K=s.QQ.algebraic_field(s.CRootOf(HP,0))
def toK(p):return K([s.QQ.convert(s.Rational(str(p[j]))) for j in reversed(range(p.degree()+1))]) if p else K.zero
def fromK(a):
 co=a.to_list();return sum(s.Rational(str(v))*h**(len(co)-1-j) for j,v in enumerate(co))
polys=[parse(f) for f in raw['polynomials']]+[parse(new['W3F'])];special=[]
for j,F in enumerate(polys):
 vals=subst(F);pol=s.Poly.from_dict({(i,):toK(v) for i,v in vals.items()},t,domain=K);special.append(pol);print('original degree',j,pol.degree(),flush=True)
G=special[0];degrees=[]
for pp in special[1:]:G=G.gcd(pp);degrees.append(G.degree())
enc=lambda pp:sum(fromK(v)*t**mon[0] for mon,v in pp.rep.to_dict().items())
print('t gcd',enc(G),'degrees',degrees,flush=True)
A=special[0].exquo(s.Poly(t,t,domain=K));B=special[1].exquo(s.Poly(t,t,domain=K));SU,SV,SD=A.gcdex(B);assert SD==s.Poly(1,t,domain=K) and SU*A+SV*B==SD

out={'status':'PASS','H2':'3*h**2+8*h-4','k':'-3*h/2-7','three_projection_substitutions_zero':True,'original_t_degrees':[pp.degree() for pp in special],'successive_t_gcd_degrees':degrees,'original_t_gcd':str(enc(G)),'bezout_after_dividing_t':[str(enc(SU)),str(enc(SV))],'bezout_identity_verified':True,'original_specialized_t_polynomials':[str(enc(pp)) for pp in special],'seconds':time.monotonic()-start}
(P/'structure_contact10_H2.json').write_text(json.dumps(out,indent=2)+'\n')
