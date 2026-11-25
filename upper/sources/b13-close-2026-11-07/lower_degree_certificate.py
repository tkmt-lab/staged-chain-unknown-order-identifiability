"""Exact Sylvester degree/valuation bounds plus finite-field gcd diagnostic.
Integer primal/dual certificates verify each assignment optimum independently.
"""
from pathlib import Path
import json,sys,time,resource,signal
from math import comb, isqrt
from flint import fmpz_mpoly_ctx,fmpz_mpoly,nmod_mpoly_ctx,nmod_mpoly
sys.set_int_max_str_digits(0);resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,)*2);signal.alarm(1750)
O=Path(__file__).resolve().parent;raw=json.loads((O/'lower_projection.json').read_text());old=json.loads((O/'lower_projection_modular.json').read_text());new=json.loads((O/'upper_contact10.json').read_text());newmod=json.loads((O/'upper_contact10_modular.json').read_text());Z=fmpz_mpoly_ctx.get(['h','k']);M=nmod_mpoly_ctx.get(['h','k'],32003)

assert all(32003 % d for d in range(2,isqrt(32003)+1)), 'modulus must be prime'

def hungarian(cost):
 n=len(cost);u=[0]*(n+1);v=[0]*(n+1);p=[0]*(n+1);way=[0]*(n+1)
 for i in range(1,n+1):
  p[0]=i;j0=0;minv=[10**9]*(n+1);used=[False]*(n+1)
  while True:
   used[j0]=True;i0=p[j0];delta=10**9;j1=0
   for j in range(1,n+1):
    if not used[j]:
     cur=cost[i0-1][j-1]-u[i0]-v[j]
     if cur<minv[j]:minv[j]=cur;way[j]=j0
     if minv[j]<delta:delta=minv[j];j1=j
   for j in range(n+1):
    if used[j]:u[p[j]]+=delta;v[j]-=delta
    else:minv[j]-=delta
   j0=j1
   if p[j0]==0:break
  while True:
   j1=way[j0];p[j0]=p[j1];j0=j1
   if j0==0:break
 assignment=[0]*n
 for j in range(1,n+1):assignment[p[j]-1]=j-1
 uu=u[1:];vv=v[1:];value=sum(cost[i][assignment[i]] for i in range(n))
 assert value==sum(uu)+sum(vv)
 assert len(set(assignment))==n
 assert all(uu[i]+vv[j]<=cost[i][j] for i in range(n) for j in range(n))
 assert max(cost[i][assignment[i]] for i in range(n))<10**9
 return {'value':value,'assignment':assignment,'row_dual':uu,'column_dual':vv}

def coeffs(poly):
 d=poly.to_dict();degree=int(poly.degrees()[1]);out=[]
 for j in range(degree,-1,-1):
  a={int(mon[0]):int(c) for mon,c in d.items() if mon[1]==j}
  if not a:out.append(None);continue
  deg=max(a);q=[sum(c*comb(i,t)*(-1)**(i-t) for i,c in a.items() if i>=t) for t in range(deg+1)]
  out.append({'degree':deg,'h':min(a),'q':next(i for i,c in enumerate(q) if c)})
 return out

def bound(g,h):
 a=coeffs(g);b=coeffs(h);m=len(a)-1;n=len(b)-1;N=m+n;mat=[]
 for i in range(n):mat.append([None]*i+a+[None]*(n-1-i))
 for i in range(m):mat.append([None]*i+b+[None]*(m-1-i))
 ans={'size':N,'degrees_k':[m,n]}
 for key in ['degree','h','q']:
  sg=-1 if key=='degree' else 1
  cost=[[10**9 if entry is None else sg*entry[key] for entry in row] for row in mat]
  result=hungarian(cost);result['bound']=sg*result['value'];ans[key]=result
 return ans

out={'status':'COMPLETED_BOUNDS','prime':32003,'prime_verified_by_integer_trial_division':True,'scope':'resultant upper degree and lower h/(h+1) valuations from exact Sylvester entries; modular gcd is promoted only if all stated criteria pass','branches':{}}
for tag in ['27','22']:
 gi=next(f for f in raw['factors']['01']['factors'] if f['degrees'][1]==int(tag));g=fmpz_mpoly(gi['polynomial'],Z);entry={};rr=[]
 for key in ['02','03']:
  hi=next(f['polynomial'] for f in raw['factors'][key]['factors'] if f['degrees'][1]>1) if key=='02' else next(f for f,m in new['factors'][key] if fmpz_mpoly(f,Z).degrees()[1]>1)
  h=fmpz_mpoly(hi,Z);bounds=bound(g,h)
  gm=nmod_mpoly(str(g),M);hm=nmod_mpoly(str(h),M)
  assert gm.degrees()[1]==g.degrees()[1] and hm.degrees()[1]==h.degrees()[1], 'input k degree dropped modulo prime'
  bounds['input_k_degrees_mod_p']=list(map(int,[gm.degrees()[1],hm.degrees()[1]]));bounds['input_k_degrees_preserved']=True
  modstr=old['branches'][tag]['res02']['polynomial'] if key=='02' else newmod['branches'][tag][key]['resultant']
  r=nmod_mpoly(modstr,M);assert r==nmod_mpoly(str(g),M).resultant(nmod_mpoly(str(h),M),'k');rr.append(r);bounds['mod_degree']=int(r.degrees()[0]);bounds['degree_preserved']=bounds['degree']['bound']==bounds['mod_degree'];entry[key]=bounds
  print(tag,key,'degree',bounds['degree']['bound'],bounds['mod_degree'],'h',bounds['h']['bound'],'q',bounds['q']['bound'],flush=True)
 gg=rr[0].gcd(rr[1]);ct,fac=gg.factor();entry['pair_mod_gcd']={'polynomial':str(gg),'factors':[[str(f),int(e)] for f,e in fac],'degree':int(gg.degrees()[0])};out['branches'][tag]=entry
 print('pair gcd',tag,[(str(f),e) for f,e in fac],flush=True)
 (O/'lower_degree_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('DONE',flush=True)
