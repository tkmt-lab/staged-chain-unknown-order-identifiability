"""Exact RUR check by integer homogeneous numerators, with no quotient-field inverses."""
import resource,signal,sys
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));signal.alarm(1750);sys.set_int_max_str_digits(0)
import flint,json,time,hashlib,math
from fractions import Fraction
from pathlib import Path
p=Path(__file__).parent;start=time.monotonic();candidate=p/(sys.argv[1] if len(sys.argv)>1 else 'lower_SHAPE_BETA_CANDIDATE.json');J=json.loads(candidate.read_text());HC=J.get('H_coefficients_low_first') or json.loads((p/'structure_H289_CANDIDATE.json').read_text())['coefficients_low_first'];Hq=[Fraction(x) for x in HC];bk=[Fraction(x) for x in (J.get('beta_k_coefficients_low_first') or J['coefficients_low_first'])];bt=[Fraction(x) for x in J['beta_t_coefficients_low_first']] if 'beta_t_coefficients_low_first' in J else None
assert Hq[-1]==1;aq=[(i+1)*Hq[i+1] for i in range(len(Hq)-1)]
def integer_poly(vals,den):return flint.fmpz_poly([int(x*den) for x in vals])
HD=math.lcm(*(x.denominator for x in Hq));H=integer_poly(Hq,HD);hc=math.gcd(*(int(x) for x in H.coeffs()));H=H//hc
D=math.lcm(*(x.denominator for vals in [aq,bk]+([bt] if bt else []) for x in vals));A=integer_poly(aq,D);B=integer_poly(bk,D);T=integer_poly(bt,D) if bt else None
common=math.gcd(*(int(x) for P in [A,B]+([T] if T is not None else []) for x in P.coeffs()));A=A//common;B=B//common;T=T//common if T is not None else None
assert H.gcd(A).degree()==0
out={'status':'RUNNING','H_degree':H.degree(),'H_input_sha256':hashlib.sha256(json.dumps(HC).encode()).hexdigest(),'candidate_path':candidate.name,'candidate_sha256':hashlib.sha256(candidate.read_bytes()).hexdigest(),'A_degree':A.degree(),'B_degree':B.degree(),'coordinate_convention':'k=B/A, optionally t=T/A; A is a nonzero rational multiple of derivative of monic H','checks':{}}

def digest(P):
 hh=hashlib.sha256()
 for cc in P.coeffs():
  v=int(cc);bb=abs(v).to_bytes(max(1,(abs(v).bit_length()+7)//8),'big');hh.update(bytes([v<0]));hh.update(len(bb).to_bytes(8,'big'));hh.update(bb)
 return hh.hexdigest()
def certify(tag,N):
 print('built',tag,'degree',N.degree(),'maxbits',max(abs(int(c)).bit_length() for c in N.coeffs()),time.monotonic()-start,flush=True)
 Q,R=divmod(N,H);assert R==0,tag+' remainder not zero';assert Q*H==N
 out['checks'][tag]={'status':'PASS','numerator_degree':N.degree(),'quotient_degree':Q.degree(),'numerator_hash':digest(N),'quotient_hash':digest(Q),'hash_encoding':'coefficients low first: sign byte, 8byte big-endian absolute-byte-length, absolute big-endian bytes','seconds':time.monotonic()-start};(p/'upper_rur_verified.json').write_text(json.dumps(out,indent=2)+'\n');print('exact divides',tag,'PASS',time.monotonic()-start,flush=True)

def bivariate(tag,expr):
 C=flint.fmpz_mpoly_ctx.get(('h','k'));F=flint.fmpz_mpoly(expr.replace('**','^'),C);d=int(F.degrees()[1]);coeff={}
 for (hh,kk),v in F.to_dict().items():coeff.setdefault(int(kk),{})[int(hh)]=int(v)
 cs={kk:flint.fmpz_poly([cc.get(i,0) for i in range(max(cc)+1)]) for kk,cc in coeff.items()};N=cs.get(d,flint.fmpz_poly());apow=flint.fmpz_poly([1])
 for i in range(d-1,-1,-1):
  apow*=A;N=N*B+cs.get(i,flint.fmpz_poly())*apow
  if i%5==0:print('Horner',tag,i,N.degree(),time.monotonic()-start,flush=True)
 certify(tag,N)

raw=json.loads((p/'lower_projection.json').read_text());G=next(f['polynomial'] for f in raw['factors']['01']['factors'] if f['degrees'][1]==27);HH=next(f['polynomial'] for f in raw['factors']['02']['factors'] if f['degrees'][1]>1)
bivariate('G27',G);bivariate('H02',HH)
# The Euler condition is tested separately and is not assumed merely from projection.
bivariate('E_G27',json.loads((p/'upper_nine_radial.json').read_text())['E_G27'])
out['status']='PASS';out['seconds']=time.monotonic()-start;(p/'upper_rur_verified.json').write_text(json.dumps(out,indent=2)+'\n');print('ALL INTEGER RUR CHECKS PASS',flush=True)
