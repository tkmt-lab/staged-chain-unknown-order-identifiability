import json
from pathlib import Path
D={tuple(range(7)):1};out={}
for j in range(4):
 assert all(max(t)<=6+j and len(t)==7 and len(set(t))==7 for t in D)
 out[str(j)]={','.join(map(str,t)):v for t,v in D.items()}
 N={}
 for t,v in D.items():
  for a in range(7):
   z=list(t);z[a]+=1
   if len(set(z))<7:continue
   inv=sum(z[i]>z[q] for i in range(7) for q in range(i+1,7));z=tuple(sorted(z));N[z]=N.get(z,0)+v*(-1)**inv
 D={t:v for t,v in N.items() if v}
print(json.dumps(out,indent=2));Path(__file__).with_suffix('.json').write_text(json.dumps({'status':'PASS','expansions':out},indent=2)+'\n')
