import resource
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
from flint import fmpz_mpoly_ctx,fmpz_mpoly,fmpq
import json,time
from pathlib import Path
p=Path(__file__).parent;ctx=fmpz_mpoly_ctx.get(['h','k']);start=time.monotonic();h0=fmpq(-5,3);k0=fmpq(-3,2)
def parse(x):return fmpz_mpoly(x.replace('**','^'),ctx)
raw=json.loads((p/'lower_projection.json').read_text());new=json.loads((p/'upper_contact10.json').read_text())
G=parse(next(f['polynomial'] for f in raw['factors']['01']['factors'] if f['degrees'][1]==22));Hs={'02':parse(next(f['polynomial'] for f in raw['factors']['02']['factors'] if f['degrees'][1]>1)),'03':next(parse(f) for f,m in new['factors']['03'] if parse(f).degrees()[1]>1)}
def ev(M):return sum((fmpq(int(co))*h0**int(i)*k0**int(j) for (i,j),co in M.to_dict().items()),fmpq(0))
assert ev(G)==0;Gh=ev(G.derivative('h'));Gk=ev(G.derivative('k'));assert Gk!=0
out={'status':'PASS','point':['-5/3','-3/2'],'G_h':str(Gh),'G_k':str(Gk),'pairs':{}}
for tag,R in Hs.items():
 assert ev(R)==0;Rh=ev(R.derivative('h'));Rk=ev(R.derivative('k'));assert Rk!=0 and Gh*Rk-Gk*Rh==0
 out['pairs'][tag]={'H_h':str(Rh),'H_k':str(Rk),'jacobian_determinant':'0','resultant_valuation_lower_bound':2}
 print('PASS',tag,flush=True)
out['seconds']=time.monotonic()-start;(p/'contact10_tangent.json').write_text(json.dumps(out,indent=2)+'\n')
