"""Independent exact Descartes interval cover; no logarithms or floating decisions."""
import resource,signal,json,time,hashlib
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));signal.alarm(1750)
from pathlib import Path
from flint import fmpz_poly,fmpq_poly,fmpq,fmpz
p=Path(__file__).parent;start=time.monotonic()
C=json.loads((p/'structure_H289_CANDIDATE.json').read_text());Hq=fmpq_poly([fmpq(x) for x in C['coefficients_low_first']]);H=Hq.numer();H=H//abs(H.content());n=H.degree();assert H[n]>0
cs=[int(H[i]) for i in range(n+1)];e=max([0]+[(abs(cs[i]).bit_length()-abs(cs[n]).bit_length()+1+n-i-1)//(n-i) for i in range(n) if cs[i]])
M=1<<e;assert all(abs(cs[i])<=abs(cs[n])*M**(n-i) for i in range(n));bound=fmpq(2*M)
def sign_at(x):
 N=fmpz(x.numerator);D=fmpz(x.denominator);v=fmpz(cs[n]);power=fmpz(1)
 for c in reversed(cs[:-1]):power*=D;v=v*N+c*power
 assert v!=0
 return 1 if v>0 else -1

def variations(a,b):
 den=int(a.denominator)*int(b.denominator);A=int(a.numerator)*(den//int(a.denominator));B=int(b.numerator)*(den//int(b.denominator))
 N=fmpz_poly([A,B]);D=fmpz_poly([den,den]);R=fmpz_poly([cs[n]]);power=fmpz_poly([1])
 for c in reversed(cs[:-1]):power*=D;R=R*N+c*power
 signs=[1 if R[i]>0 else -1 for i in range(R.degree()+1) if R[i]]
 return sum(signs[i]!=signs[i-1] for i in range(1,len(signs)))
trace=[];ones=[]
def cover(a,b,kind):
 if a==b:return 0
 assert a<b;stack=[(a,b)];count=0
 while stack:
  lo,hi=stack.pop();v=variations(lo,hi);entry={'lo':str(lo),'hi':str(hi),'variations':v,'kind':kind};trace.append(entry)
  if v<=1:
   if v:ones.append([str(lo),str(hi)]);count+=v
  else:
   mid=(lo+hi)/2;entry['mid_sign']=sign_at(mid);stack.extend([(mid,hi),(lo,mid)])
  if len(trace)%50==0:print('nodes',len(trace),'roots',len(ones),'seconds',time.monotonic()-start,flush=True)
 return count
J=json.loads((p/'structure_H289_ISOLATION.json').read_text());intervals=J['factors'][0]['intervals'];previous=fmpq(-1);assert sign_at(previous)
for index,row in enumerate(intervals):
 a,b=fmpq(row['lo']),fmpq(row['hi']);assert previous<=a<b
 assert sign_at(a)*sign_at(b)==-1
 assert cover(previous,a,'gap')==0
 assert cover(a,b,'root')==1
 previous=b;print('root',index+1,'PASS',time.monotonic()-start,flush=True)
assert previous<bound;assert cover(previous,bound,'tail')==0
assert len(ones)==28
out={'status':'PASS','method':'integer homogeneous Mobius transforms, Descartes variations, rational bisection and exact power-of-two global root bound','root_count_above_minus_one':28,'global_absolute_root_bound':str(bound),'bound_M':M,'trace':trace,'isolated_root_intervals':ones,'seconds':time.monotonic()-start}
(p/'audit_real_roots.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS complete 28-root cover',time.monotonic()-start,flush=True)
