"""Exact q-adic determinant valuation of Sylvester matrices, q=h+1.
Truncated rational series + unit Schur complements. No floating-point decisions.
"""
from pathlib import Path
import json,sys,time,resource,signal,hashlib
from math import comb
from flint import fmpz_mpoly_ctx,fmpz_mpoly,fmpq_poly
sys.set_int_max_str_digits(0);resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,)*2);signal.alarm(1750)
O=Path(__file__).resolve().parent;raw=json.loads((O/'lower_projection.json').read_text());new=json.loads((O/'upper_contact10.json').read_text());Z=fmpz_mpoly_ctx.get(['h','k']);start=time.time()

def coeffs(poly,precision):
 d=poly.to_dict();degree=int(poly.degrees()[1]);out=[]
 for j in range(degree,-1,-1):
  a={int(mon[0]):int(c) for mon,c in d.items() if mon[1]==j}
  values=[sum(c*comb(i,t)*(-1)**(i-t) for i,c in a.items() if i>=t) for t in range(precision)]
  out.append(fmpq_poly(values))
 return out

def inverse(p,N):
 assert p[0]
 u=fmpq_poly([1/p[0]]);size=1
 while size<N:
  size=min(2*size,N);u=(u*(2-p*u)).truncate(size)
 assert (p*u).truncate(N)==fmpq_poly([1])
 return u

def certificate(g,h,N,target):
 a=coeffs(g,N);b=coeffs(h,N);m=len(a)-1;n=len(b)-1;zero=fmpq_poly([]);matrix=[]
 for i in range(n):matrix.append([zero]*i+a+[zero]*(n-1-i))
 for i in range(m):matrix.append([zero]*i+b+[zero]*(m-1-i))
 precision=N;total=0;trace=[];step=0
 while matrix:
  size=len(matrix)
  # Exact truncation implies every true omitted entry is divisible by q^precision.
  nonzero=[(i,j,p) for i,row in enumerate(matrix) for j,p in enumerate(row) if p]
  if not nonzero:
   total+=size*precision;trace.append({'dimension':size,'all_zero_mod_q_power':precision,'total_lower_bound':total});break
  valuation=min(next(t for t in range(precision) if p[t]) for i,j,p in nonzero)
  if valuation:
   total+=size*valuation;precision-=valuation
   matrix=[[p.right_shift(valuation) for p in row] for row in matrix]
   trace.append({'dimension':size,'extracted_from_each_row':valuation,'total_lower_bound':total,'remaining_precision':precision})
   print('factor',size,valuation,'total',total,'precision',precision,flush=True)
   if total>=target:break
  choices=[(i,j,p[0]) for i,row in enumerate(matrix) for j,p in enumerate(row) if p[0]]
  # Smaller exact rational constants usually reduce coefficient growth.
  i,j,p0=min(choices,key=lambda z:abs(z[2].numerator).bit_length()+z[2].denominator.bit_length())
  matrix[0],matrix[i]=matrix[i],matrix[0]
  for row in matrix:row[0],row[j]=row[j],row[0]
  pivot=matrix[0][0];inv=inverse(pivot,precision)
  row=matrix[0][1:];column=[(r[0]*inv).truncate(precision) for r in matrix[1:]]
  rest=[]
  for oldrow,multiplier in zip(matrix[1:],column):
   rest.append([(oldrow[l+1]-multiplier*row[l]).truncate(precision) for l in range(size-1)])
  trace.append({'dimension':size,'unit_pivot_row':i,'unit_pivot_column':j,'unit_pivot_constant':str(p0),'total_lower_bound':total,'precision':precision})
  matrix=rest;step+=1
  if step%10==0:print('Schur',step,'remaining',len(matrix),'total',total,'seconds',time.time()-start,flush=True)
 return {'initial_precision':N,'target':target,'determinant_valuation_lower_bound':total,'passed':total>=target,'trace':trace}

out={'scope':'exact lower bound on v_(h+1)(Res_k(G,H)); formal series arithmetic over QQ','branches':{}}
tags=sys.argv[1:] or ['27','22']
for tag in tags:
 gi=next(f for f in raw['factors']['01']['factors'] if f['degrees'][1]==int(tag));g=fmpz_mpoly(gi['polynomial'],Z);target=98 if tag=='27' else 92
 for key in ['02','03']:
  hi=next(f['polynomial'] for f in raw['factors'][key]['factors'] if f['degrees'][1]>1) if key=='02' else next(f for f,m in new['factors'][key] if fmpz_mpoly(f,Z).degrees()[1]>1)
  print('START',tag,key,flush=True);result=certificate(g,fmpz_mpoly(hi,Z),12,target);out['branches'].setdefault(tag,{})[key]=result;out['seconds']=time.time()-start
  (O/'lower_qadic_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print('RESULT',tag,key,result['passed'],result['determinant_valuation_lower_bound'],flush=True)
  assert result['passed'],'increase series precision'
out['status']='PASS';(O/'lower_qadic_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS',time.time()-start,flush=True)
