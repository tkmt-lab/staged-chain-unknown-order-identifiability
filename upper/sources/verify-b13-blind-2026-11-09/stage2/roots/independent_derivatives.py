"""Stage2 new implementation: derivative jets, fixed dyadic interval arithmetic.
Input RUR and isolation from source; no source classification or coefficient solver imported.
"""
import json,math,time,resource,signal,hashlib
from pathlib import Path
from fractions import Fraction as Q
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));signal.alarm(1750)
P=Path(__file__).parent;START=time.monotonic();BITS=6144;S=1<<BITS
class Box:
 def __init__(self,l,u=None):
  if u is None:
   q=Q(l);self.l=q.numerator*S//q.denominator;self.u=-((-q.numerator*S)//q.denominator)
  else:self.l,self.u=l,u
  assert self.l<=self.u
 def __add__(a,b):
  if not isinstance(b,Box):b=Box(b)
  return Box(a.l+b.l,a.u+b.u)
 __radd__=__add__
 def __neg__(a):return Box(-a.u,-a.l)
 def __sub__(a,b):return a+-b if isinstance(b,Box) else a+Box(-b)
 def __rsub__(a,b):return -a+b
 def __mul__(a,b):
  if not isinstance(b,Box):b=Box(b)
  v=[x*y for x in [a.l,a.u] for y in [b.l,b.u]]
  return Box(min(v)//S,-((-max(v))//S))
 __rmul__=__mul__
 def __truediv__(a,b):
  if not isinstance(b,Box):b=Box(b)
  assert b.l*b.u>0
  inv=Box(min(S*S//b.l,S*S//b.u),max(-((-S*S)//b.l),-((-S*S)//b.u)))
  return a*inv
 def sgn(a):return 1 if a.l>0 else -1 if a.u<0 else 0
 def data(a):return {'lo_scaled':str(a.l),'hi_scaled':str(a.u),'integer_enclosure':[a.l//S,-((-a.u)//S)],'sign':a.sgn()}
def poly(cs,x):
 z=Box(0)
 for c in reversed(cs):z=z*x+Box(c)
 return z
def exact_poly(cs,x):
 z=Q(0)
 for c in reversed(cs):z=z*x+Q(c)
 return z
C=json.loads((P/'lower_RUR_CANDIDATE.json').read_text());H=list(map(int,json.loads((P/'structure_H289_FACTORS.json').read_text())['primitive_coefficients_low_first']));HP=[i*H[i] for i in range(1,len(H))]
K=list(map(Q,C['beta_k_coefficients_low_first']));T=list(map(Q,C['beta_t_coefficients_low_first']));roots=json.loads((P/'structure_H289_SIGNS.json').read_text())['roots']
# Nontrivial unit tests for enclosure division, especially negative denominator.
for a,b,c,d in [(1,2,3,4),(-5,-2,2,3),(1,4,-3,-1),(-7,3,-9,-2)]:
 z=Box(a*S,b*S)/Box(c*S,d*S)
 for x in [a,b]:
  for y in [c,d]:assert Q(z.l,S)<=Q(x,y)<=Q(z.u,S)
results=[]
for row in roots:
 lo,hi=Q(row['h_lo']),Q(row['h_hi']);assert exact_poly(H,lo)*exact_poly(H,hi)<0
 h=Box(lo.numerator*S//lo.denominator,-((-hi.numerator*S)//hi.denominator));hp=poly(HP,h);assert hp.sgn()
 k=poly(K,h)*H[-1]/hp;t=poly(T,h)*H[-1]/hp;q=1+h;r=q*t
 assert q.l>0
 out={'id':row['index'],'k':k.data(),'r':r.data()}
 if k.l>0:out['reason']='k>0'
 elif r.u<0:out['reason']='r<0'
 else:
  assert k.u<0 and r.l>0
  # Direct derivatives of S'= -SX, X'=kX+SX, y=SX.
  sd=[Box(1)];xd=[r];yd=[]
  for n in range(9):
   yn=sum((math.comb(n,j)*sd[j]*xd[n-j] for j in range(n+1)),Box(0));yd.append(yn)
   sd.append(-yn);xd.append(k*xd[n]+yn)
  zd=[]
  for n in range(9):zd.append((yd[n]-sum((math.comb(n,j)*sd[j]*zd[n-j] for j in range(1,n+1)),Box(0)))/q)
  M=[[zd[i+2],zd[i+1],zd[i],-yd[i+2],-yd[i+1],-yd[i]] for i in range(6)];b=[-zd[i+3] for i in range(6)];piv=[]
  # Forward elimination followed by back substitution (source uses Gauss-Jordan).
  for col in range(6):
   pivot=next(j for j in range(col,6) if M[j][col].sgn());M[col],M[pivot]=M[pivot],M[col];b[col],b[pivot]=b[pivot],b[col];piv.append(M[col][col].data())
   for j in range(col+1,6):
    mul=M[j][col]/M[col][col]
    for c in range(col+1,6):M[j][c]=M[j][c]-mul*M[col][c]
    b[j]=b[j]-mul*b[col]
  sol=[None]*6
  for j in reversed(range(6)):sol[j]=(b[j]-sum((M[j][c]*sol[c] for c in range(j+1,6)),Box(0)))/M[j][j]
  d2,d1,d0,g,v,w=sol;gap=d2*d2-3*d1
  m0=(zd[2]+d2*zd[1]+d1*zd[0]-g*yd[1]-v*yd[0])/g
  tests={'g':g,'d2':d2,'trace_gap':gap,'m0':m0}
  reason=next((name for name,value in tests.items() if value.u<0),None);assert reason is not None
  out.update(reason=reason,tests={name:value.data() for name,value in tests.items()},pivots=piv)
 results.append(out);print(row['index'],out['reason'],flush=True)
assert len(results)==28
result={'status':'PASS','method':'new derivative recurrence and forward elimination/back substitution, dyadic outward intervals','bits':BITS,'all_28_excluded':True,'source_input_sha256':{n:hashlib.sha256((P/n).read_bytes()).hexdigest() for n in ['lower_RUR_CANDIDATE.json','structure_H289_FACTORS.json','structure_H289_SIGNS.json']},'roots':results,'seconds':time.monotonic()-START}
(P/'independent_derivatives.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS',result['seconds'],flush=True)
