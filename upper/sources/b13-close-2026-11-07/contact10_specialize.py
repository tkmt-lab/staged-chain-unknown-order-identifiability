import resource
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
from flint import fmpz_mpoly_ctx,fmpz_mpoly,fmpq_poly
import json,time
from pathlib import Path
p=Path(__file__).parent;ctx=fmpz_mpoly_ctx.get(['h','k']);start=time.monotonic()
def parse(x):return fmpz_mpoly(x.replace('**','^'),ctx)
raw=json.loads((p/'lower_projection.json').read_text());new=json.loads((p/'upper_contact10.json').read_text())
G=parse(next(f['polynomial'] for f in raw['factors']['01']['factors'] if f['degrees'][1]==22));Hs={'02':parse(next(f['polynomial'] for f in raw['factors']['02']['factors'] if f['degrees'][1]>1)),'03':next(parse(f) for f,m in new['factors']['03'] if parse(f).degrees()[1]>1)}
out={}
for fs in ['3*h^2+8*h-4','14*h^2+40*h+29','3*h+5']:
 Hm=parse(fs);hd=Hm.degrees()[0];hdict=Hm.to_dict();H=fmpq_poly([int(hdict.get((i,0),0)) for i in range(hd+1)])
 def trim(P):
  while P and not P[-1]:P.pop()
  return P
 def poly(M):
  md=M.to_dict();P=[]
  for j in range(M.degrees()[1]+1):P.append(fmpq_poly([int(md.get((i,j),0)) for i in range(M.degrees()[0]+1)])%H)
  return trim(P)
 def monic(P):
  if not P:return []
  gcd,inv,_=P[-1].xgcd(H);assert gcd==1
  return [(x*inv)%H for x in P]
 def rem(A,B):
  A=A[:];B=monic(B)
  while len(A)>=len(B) and A:
   fac=A[-1];shift=len(A)-len(B)
   for i,x in enumerate(B):A[i+shift]=(A[i+shift]-fac*x)%H
   trim(A)
  return A
 def gcd(A,B):
  while B:A,B=B,rem(A,B)
  return monic(A)
 A=poly(G);res={}
 for tag,R in Hs.items():
  B=poly(R);Z=gcd(A,B);da=int(G.degrees()[1])-len(A)+1;db=int(R.degrees()[1])-len(B)+1;nullity=len(Z)-1+min(da,db)
  res[tag]={'G_original_degree':int(G.degrees()[1]),'H_original_degree':int(R.degrees()[1]),'G_specialized_degree':len(A)-1,'H_specialized_degree':len(B)-1,'gcd_degree':len(Z)-1,'gcd_coefficients_low_first':[str(x) for x in Z],'sylvester_nullity':nullity}
  print(fs,tag,'drops',da,db,'gcd',len(Z)-1,'nullity',nullity,time.monotonic()-start,flush=True)
 out[fs]=res
 (p/'contact10_specialize.json').write_text(json.dumps(out,indent=2)+'\n')
