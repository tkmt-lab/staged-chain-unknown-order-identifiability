"""Independent rational-bracket refinement and dyadic integer sign checks.
RUR coverage/identities are certified separately by the upper担当.
"""
import resource,signal,time,json,sys,math,hashlib
from pathlib import Path
from fractions import Fraction
from flint import fmpz_poly,fmpq_poly,fmpq
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));signal.alarm(1750);sys.set_int_max_str_digits(0)
P=Path(__file__).parent;start=time.monotonic();ISO=json.loads((P/'structure_H289_ISOLATION.json').read_text());assert ISO['status']=='PASS' and ISO['integer_only_root_bounds'];J=json.loads((P/'lower_RUR_CANDIDATE.json').read_text());FAC=json.loads((P/'structure_H289_FACTORS.json').read_text());H=fmpz_poly(list(map(int,FAC['primitive_coefficients_low_first'])));HP=H.derivative();hf=fmpq_poly(H);hpf=fmpq_poly(HP)
def integer_coefficients(vals):
 vals=list(map(Fraction,vals));den=math.lcm(*(x.denominator for x in vals));return [int(x*den) for x in vals],den
BK,DK=integer_coefficients(J['beta_k_coefficients_low_first']);BT,DT=integer_coefficients(J['beta_t_coefficients_low_first']);HC=[int(HP[i]) for i in range(HP.degree()+1)]
sign=lambda x:1 if x>0 else -1 if x<0 else 0
def floor_scaled(x,p):return int(x.numerator)*(1<<p)//int(x.denominator)
def refine(lo,hi,bits):
 original_lo,original_hi=lo,hi;fl=hf(lo);fh=hf(hi);assert fl*fh<0;x=(lo+hi)/2;iters=0
 for p in [64,128,256,512,1024,2048,bits]:
  if hi-lo<fmpq(4,1<<p) and p<bits:continue
  eps=fmpq(1,1<<p)
  for j in range(80):
   y=hf(x);dy=hpf(x);z=x-y/dy if dy else (lo+hi)/2;z=fmpq(floor_scaled(z,p),1<<p)
   if not lo<z<hi:z=(lo+hi)/2
   yz=hf(z);assert yz
   if sign(yz)==sign(fl):lo,fl=z,yz
   else:hi,fh=z,yz
   z=fmpq(floor_scaled(z,p),1<<p);a=z-eps;b=z+eps
   if original_lo<a<b<original_hi:
    fa,fb=hf(a),hf(b)
    if fa*fb<0:lo,hi,fl,fh,x=a,b,fa,fb,z;iters+=j+1;break
   x=z
  else:raise ArithmeticError(('rational bracket refinement stalled',p))
 assert original_lo<lo<hi<original_hi and hf(lo)*hf(hi)<0
 return lo,hi,iters

def poly_interval(coeffs,a,b,p):
 S=1<<p;lo=hi=0
 for c in reversed(coeffs):
  vv=[lo*a,lo*b,hi*a,hi*b];lo=min(vv)//S+c*S;hi=-((-max(vv))//S)+c*S
 return lo,hi

def interval_sign(v):return 1 if v[0]>0 else -1 if v[1]<0 else 0
out={'status':'RUNNING','method':'exact rational Newton proposals accepted only by sign-changing subintervals; independent fixed dyadic integer Horner enclosures','dyadic_bits':3072,'RUR_source_sha256':hashlib.sha256((P/'lower_RUR_CANDIDATE.json').read_bytes()).hexdigest(),'scope':'Sign classification of the supplied RUR; exact RUR identities and original-system coverage are separate obligations','roots':[]}
for index,iv in enumerate(ISO['factors'][0]['intervals']):
 lo,hi,nit=refine(fmpq(iv['lo']),fmpq(iv['hi']),3072);a=floor_scaled(lo,3072);b=floor_scaled(hi,3072);assert fmpq(a,1<<3072)==lo and fmpq(b,1<<3072)==hi
 si=[interval_sign(poly_interval(cs,a,b,3072)) for cs in [HC,BK,BT]];assert all(si),('uncertain sign',index,si)
 sk=si[0]*si[1];st=si[0]*si[2];entry={'index':index,'h_lo':str(lo),'h_hi':str(hi),'refinement_steps':nit,'Hprime_sign':si[0],'beta_k_sign':si[1],'beta_t_sign':si[2],'k_sign':sk,'t_sign':st,'passes_k_negative_t_positive':sk<0 and st>0};out['roots'].append(entry)
 print('root',index,'k',sk,'t',st,'passes',entry['passes_k_negative_t_positive'],'sec',time.monotonic()-start,flush=True)
 out['seconds']=time.monotonic()-start;(P/'structure_H289_SIGNS.json').write_text(json.dumps(out,indent=2)+'\n')
out['status']='PASS';out['passes_count']=sum(v['passes_k_negative_t_positive'] for v in out['roots']);out['seconds']=time.monotonic()-start;(P/'structure_H289_SIGNS.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS',out['passes_count'],'of',len(out['roots']),'pass k/t signs',flush=True)
