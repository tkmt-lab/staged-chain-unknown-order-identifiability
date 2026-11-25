"""Verify F=WF=W2F exactly using k=beta_k/H', t=beta_t/H'."""
import resource,signal,sys
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));signal.alarm(1750);sys.set_int_max_str_digits(0)
import flint,json,time,hashlib,math
from fractions import Fraction
from pathlib import Path
p=Path(__file__).parent;start=time.monotonic();source=p/'lower_RUR_CANDIDATE.json';J=json.loads(source.read_text());old=json.loads((p/'lower_SHAPE_BETA_CANDIDATE.json').read_text());HC=J['H_coefficients_low_first'];Hq=[Fraction(x) for x in HC];bk=[Fraction(x) for x in J['beta_k_coefficients_low_first']];bt=[Fraction(x) for x in J['beta_t_coefficients_low_first']];assert bk==list(map(Fraction,old['coefficients_low_first']));assert Hq==list(map(Fraction,old['H_coefficients_low_first']));assert Hq[-1]==1;aq=[(i+1)*Hq[i+1] for i in range(len(Hq)-1)]
def ip(vals,d):
 assert all((v*d).denominator==1 for v in vals)
 return flint.fmpz_poly([int(v*d) for v in vals])
HD=math.lcm(*(x.denominator for x in Hq));H=ip(Hq,HD);H=H//math.gcd(*(int(c) for c in H.coeffs()));D=math.lcm(*(x.denominator for vals in [aq,bk,bt] for x in vals));A=ip(aq,D);B=ip(bk,D);T=ip(bt,D);cc=math.gcd(*(int(x) for P in [A,B,T] for x in P.coeffs()));A=A//cc;B=B//cc;T=T//cc;assert H.gcd(A).degree()==0
C=flint.fmpz_mpoly_ctx.get(('h','k','t'));raw=json.loads((p/'lower_projection.json').read_text());out={'status':'RUNNING','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'source':source.name,'H_squarefree_and_common_denominator_coprime':True,'checks':{}}
def powers(P,n):
 out=[flint.fmpz_poly([1])]
 for _ in range(n):out.append(out[-1]*P)
 return out
# All three original contact polynomials have k-degree9 and t-degree7.
AP=powers(A,9);TP=powers(T,7);TM=[TP[j]*AP[7-j] for j in range(8)]
for tag,expr in zip(['F','WF','W2F'],raw['polynomials']):
 F=flint.fmpz_mpoly(expr,C);dk=int(F.degrees()[1]);dt=int(F.degrees()[2]);assert (dk,dt)==(9,7);coeff={}
 for (hh,kk,tt),v in F.to_dict().items():coeff.setdefault((int(kk),int(tt)),{})[int(hh)]=int(v)
 CT=[]
 for j in range(dk+1):
  P=flint.fmpz_poly()
  for ell in range(dt+1):
   co=coeff.get((j,ell),{});ch=flint.fmpz_poly([co.get(i,0) for i in range(max(co)+1)]) if co else flint.fmpz_poly();P+=ch*TM[ell]
  CT.append(P)
 N=CT[-1]
 for j in range(dk-1,-1,-1):N=N*B+CT[j]*AP[dk-j]
 print('built',tag,'degree',N.degree(),'maxbits',max(abs(int(c)).bit_length() for c in N.coeffs()),time.monotonic()-start,flush=True)
 Q,R=divmod(N,H);assert R==0;assert Q*H==N
 out['checks'][tag]={'status':'PASS','numerator_degree':N.degree(),'quotient_degree':Q.degree(),'homogenizing_exponent_k_plus_t':dk+dt,'seconds':time.monotonic()-start};(p/'upper_rur_contact_verified.json').write_text(json.dumps(out,indent=2)+'\n');print('exact divides',tag,'PASS',time.monotonic()-start,flush=True)
out['status']='PASS';out['seconds']=time.monotonic()-start;(p/'upper_rur_contact_verified.json').write_text(json.dumps(out,indent=2)+'\n');print('ALL CONTACT RUR CHECKS PASS',flush=True)
