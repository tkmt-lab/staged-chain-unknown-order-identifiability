import resource,sys
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));sys.set_int_max_str_digits(0)
import flint,json,time
from pathlib import Path
p=Path(__file__).parent;start=time.monotonic();C=flint.fmpz_mpoly_ctx.get(('h','k','t'));h,k,t=C.gens();J=json.loads((p/'upper_nine_subres.json').read_text());F=flint.fmpz_mpoly(next(f for f,d in zip(J['subresultants'],J['degrees']) if d==1),C)
A=C.from_dict({(hh,kk,0):v for (hh,kk,tt),v in F.to_dict().items() if tt==1});B=C.from_dict({(hh,kk,0):v for (hh,kk,tt),v in F.to_dict().items() if tt==0});content=A.gcd(B);co,fac=content.factor();a=A//content;b=B//content;print('content',[(f.degrees(),m) for f,m in fac],flush=True);print('a,b',a.degrees(),b.degrees(),len(a.to_dict()),len(b.to_dict()),time.monotonic()-start,flush=True)
out={'content':str(content),'content_factors':[[str(f),int(m)] for f,m in fac],'a':str(a),'b':str(b),'seconds':time.monotonic()-start};(p/'upper_nine_linear.json').write_text(json.dumps(out,indent=2)+'\n')
