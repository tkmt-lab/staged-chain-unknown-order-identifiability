"""Exact algebra and exterior-derivative checks supporting G22 closure."""
import resource,signal,time,json
from pathlib import Path
from collections import defaultdict
from flint import fmpz_mpoly_ctx,fmpz_mpoly
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));signal.alarm(1750)
P=Path(__file__).parent;start=time.monotonic();C=fmpz_mpoly_ctx.get(('h','k'));h,k=C.gens();C3=fmpz_mpoly_ctx.get(('h','k','t'));raw=json.loads((P/'lower_projection.json').read_text());F=fmpz_mpoly(raw['polynomials'][0],C3);d=F.degrees()[2];lead=C.from_dict({(i,j):co for (i,j,n),co in F.to_dict().items() if n==d});Q=5*h**4-4*h**3+3*h*h-2*h+1
quo,rem=divmod(lead,h**3*k*(1+h)*Q);assert not rem and len(quo)==1 and quo.degrees()==(0,0)
assert 5*Q==5*(h-1)**2+h*h*(5*h-2)**2+6*h*h
D={tuple(range(7)):1};expansions=[]
for step in range(3):
 E=defaultdict(int)
 for inds,co in D.items():
  for j in range(7):
   x=list(inds);x[j]+=1
   if len(set(x))<7:continue
   sign=(-1)**sum(x[a]>x[b] for a in range(7) for b in range(a+1,7));E[tuple(sorted(x))]+=co*sign
 D={key:v for key,v in E.items() if v};expansions.append({' '.join(map(str,key)):v for key,v in D.items()})
assert D=={(0,1,2,3,5,6,7):1,(0,1,2,3,4,6,8):2,(0,1,2,3,4,5,9):1}
lin=json.loads((P/'upper_nine_linear.json').read_text());assert sorted(lin['content_factors'])==sorted([['h',6],['k',18],['h + 1',40]])
out={'status':'PASS','F_t_degree':int(d),'leading_coefficient_nonzero_constant':str(quo),'leading_coefficient_factored':'constant*h^3*k*(1+h)*(5h^4-4h^3+3h^2-2h+1)','positive_quartic_SOS_identity':True,'D1_D2_D3_determinant_expansions':expansions,'S1_removed_content_is_unit_on_physical_domain':True,'seconds':time.monotonic()-start}
(P/'structure_G22_geometry.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS leading coefficient, SOS, exterior derivative identity, S1 units')
