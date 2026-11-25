"""Independent ordinary-Taylor enclosure of the ten physical H289 candidates.
Does not import the previous contact implementation or its derivative recursion.
"""
import resource,signal,time,json,sys,math,hashlib
from pathlib import Path
from fractions import Fraction
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));signal.alarm(1750);sys.set_int_max_str_digits(0)
P=Path(__file__).parent;start=time.monotonic();BITS=4608;SCALE=1<<BITS
class Iv:
 def __init__(self,lo,hi):self.lo=int(lo);self.hi=int(hi);assert self.lo<=self.hi
 @staticmethod
 def exact(v):
  if isinstance(v,Iv):return v
  q=Fraction(v);return Iv(q.numerator*SCALE//q.denominator,-((-q.numerator*SCALE)//q.denominator))
 def __add__(self,other):
  z=Iv.exact(other);return Iv(self.lo+z.lo,self.hi+z.hi)
 __radd__=__add__
 def __neg__(self):return Iv(-self.hi,-self.lo)
 def __sub__(self,other):return self+-Iv.exact(other)
 def __rsub__(self,other):return Iv.exact(other)+-self
 def __mul__(self,other):
  z=Iv.exact(other);v=[self.lo*z.lo,self.lo*z.hi,self.hi*z.lo,self.hi*z.hi];return Iv(min(v)//SCALE,-((-max(v))//SCALE))
 __rmul__=__mul__
 def __truediv__(self,other):
  z=Iv.exact(other);assert z.lo>0 or z.hi<0,'interval denominator includes zero'
  lows=[];highs=[]
  for a in (self.lo,self.hi):
   for b in (z.lo,z.hi):lows.append(a*SCALE//b);highs.append(-((-a*SCALE)//b))
  return Iv(min(lows),max(highs))
 def nonzero(self):return self.lo>0 or self.hi<0
 def sign(self):return 1 if self.lo>0 else -1 if self.hi<0 else 0
 def data(self):return {'lo_scaled':str(self.lo),'hi_scaled':str(self.hi),'strict_sign':self.sign(),'integer_enclosure':[self.lo//SCALE,-((-self.hi)//SCALE)]}
def horner(coeffs,x):
 y=Iv.exact(0)
 for c in reversed(coeffs):y=y*x+c
 return y
def integer_coeffs(vals):
 vals=list(map(Fraction,vals));den=math.lcm(*(v.denominator for v in vals));return [int(v*den) for v in vals],den
J=json.loads((P/'lower_RUR_CANDIDATE.json').read_text());Knum,Kden=integer_coeffs(J['beta_k_coefficients_low_first']);Tnum,Tden=integer_coeffs(J['beta_t_coefficients_low_first']);H=list(map(int,json.loads((P/'structure_H289_FACTORS.json').read_text())['primitive_coefficients_low_first']));HP=[(i+1)*H[i+1] for i in range(len(H)-1)];LC=H[-1]
roots=json.loads((P/'structure_H289_SIGNS.json').read_text());assert roots['status']=='PASS';selected=[r for r in roots['roots'] if r['passes_k_negative_t_positive']];assert [r['index'] for r in selected]==[0,3,4,6,10,14,15,16,17,18]
out={'status':'RUNNING','method':'ordinary Taylor coefficients of original left ODE; interval division for y/q; independently generated derivatives; interval Gaussian elimination','dyadic_bits':BITS,'ordinary_taylor_order':10,'input_RUR_sha256':hashlib.sha256((P/'lower_RUR_CANDIDATE.json').read_bytes()).hexdigest(),'input_root_brackets_sha256':hashlib.sha256((P/'structure_H289_SIGNS.json').read_bytes()).hexdigest(),'scope':'Conditional on exact RUR coverage/identities proved separately; excludes all ten physical RUR candidates by model-independent necessary signs','roots':[]}
for root in selected:
 a=Fraction(root['h_lo']);b=Fraction(root['h_hi']);h=Iv(a.numerator*SCALE//a.denominator,-((-b.numerator*SCALE)//b.denominator));hp=horner(HP,h);assert hp.nonzero();k=horner(Knum,h)*LC/Kden/hp;t=horner(Tnum,h)*LC/Tden/hp;assert k.hi<0 and t.lo>0
 q0=1+h;assert q0.lo>0;r=q0*t;SS=[Iv.exact(1)];XX=[r];yy=[]
 # Ordinary Taylor coefficients: no binomial derivative recurrence.
 for n in range(11):
  yy.append(sum((SS[i]*XX[n-i] for i in range(n+1)),Iv.exact(0)))
  SS.append(-yy[n]/(n+1));XX.append((k*XX[n]+yy[n])/(n+1))
 qq=SS[:];qq[0]=q0;zz=[]
 for n in range(11):zz.append((yy[n]-sum((qq[i]*zz[n-i] for i in range(1,n+1)),Iv.exact(0)))/q0)
 yd=[math.factorial(i)*yy[i] for i in range(11)];zd=[math.factorial(i)*zz[i] for i in range(11)]
 M=[[zd[i+2],zd[i+1],zd[i],-yd[i+2],-yd[i+1],-yd[i],-zd[i+3]] for i in range(6)];pivots=[]
 for col in range(6):
  row=next(j for j in range(col,6) if M[j][col].nonzero());M[col],M[row]=M[row],M[col];pivot=M[col][col];pivots.append({'column':col,'chosen_row':row,'interval':pivot.data()});M[col]=[v/pivot for v in M[col]]
  for j in range(6):
   if j!=col:
    mult=M[j][col];M[j]=[a-mult*b for a,b in zip(M[j],M[col])]
 d2,d1,d0,g,v,w=[M[j][-1] for j in range(6)];trace=d2*d2-3*d1;assert g.nonzero();m0=(zd[2]+d2*zd[1]+d1*zd[0]-g*yd[1]-v*yd[0])/g
 if g.hi<0:reason='g<0'
 elif d2.hi<0:reason='d2<0'
 elif trace.hi<0:reason='d2^2-3*d1<0'
 elif m0.hi<0:reason='free_response_m0<0'
 else:raise AssertionError(('not excluded',root['index']))
 entry={'index':root['index'],'all_six_pivots_nonzero':True,'pivots':pivots,'coefficients':{name:z.data() for name,z in zip(['d2','d1','d0','g','v','w'],[d2,d1,d0,g,v,w])},'trace_gap':trace.data(),'free_response_m0':m0.data(),'exclusion':reason};out['roots'].append(entry);out['seconds']=time.monotonic()-start;(P/'structure_TAYLOR_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n');print('excluded',root['index'],reason,'sec',time.monotonic()-start,flush=True)
out['status']='PASS';out['excluded_count']=len(out['roots']);out['seconds']=time.monotonic()-start;(P/'structure_TAYLOR_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS all ten excluded by independent ordinary Taylor intervals',flush=True)
