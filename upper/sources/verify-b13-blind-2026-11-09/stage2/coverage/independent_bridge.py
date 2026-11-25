"""Direct exact checks for input edges omitted by the supplied replay sequence.
No new resultant factorization or projection search is needed.
"""
from pathlib import Path
import json,time,resource,sys,signal,hashlib
import flint
sys.set_int_max_str_digits(0);resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,)*2);signal.alarm(1750)
O=Path(__file__).resolve().parent;P=O.parents[2]/'b13-close-2026-11-07';start=time.monotonic()
read=lambda n:json.loads((P/n).read_text())
raw=read('lower_projection.json');new=read('upper_contact10.json');rad=read('upper_nine_radial.json');piv=read('structure_nine_pivots.json');base=json.loads((P.parent/'b13-lower3-2026-11-07/elim_general.json').read_text())
C=flint.fmpz_mpoly_ctx.get(('h','k','t'));h,k,t=C.gens();parse=lambda x:flint.fmpz_mpoly(x.replace('**','^'),C)
F=[parse(x) for x in raw['polynomials']];W=lambda x:-h*(1+h)*x.derivative('h')-k*(1+h)*x.derivative('k')-(t+1+k)*x.derivative('t')
assert W(F[0])==F[1] and W(F[1])==F[2];W3=W(F[2]);assert W3==parse(new['W3F'])
actual=F[0].resultant(W3,'t');assert actual==parse(new['resultants']['03'])
prod=C.constant(1)
for expr,power in new['factors']['03']:prod*=parse(expr)**power
quo,rem=divmod(actual,prod);assert not rem and quo.is_constant() and quo
# Check that factor selection drops only h,k,h+1, all nonzero physically.
removed={}
for key,allowed in [('01',[22,27]),('02',[49])]:
 selected=[];units=[]
 for f in raw['factors'][key]['factors']:
  z=parse(f['polynomial'])
  if int(z.degrees()[1]) in allowed:selected.append(int(z.degrees()[1]))
  else:assert z in [h,k,h+1];units.append([str(z),int(f['power'])])
 assert sorted(selected)==allowed;removed[key]=units
selected=[];units=[]
for expr,power in new['factors']['03']:
 z=parse(expr)
 if int(z.degrees()[1])==49:selected.append(z)
 else:assert z in [h,k,h+1];units.append([str(z),int(power)])
assert len(selected)==1;removed['03']=units
G=parse(next(f['polynomial'] for f in raw['factors']['01']['factors'] if f['degrees'][1]==27));assert h*G.derivative('h')+k*G.derivative('k')==parse(rad['E_G27'])
A=parse(base['A'].replace('r','t'));Asub=sum((int(v)*h**i*k**j*t**ell*(1+h)**ell for (i,j,ell),v in A.to_dict().items()),C.constant(0));assert Asub==parse(piv['transformed']['A'])*parse(piv['removed']['A'])
# All removed factors of A are units on h!=0,k<0,h+1>0,t>0.
ct,factors=parse(piv['removed']['A']).factor();assert ct and all(f in [h,k,t,h+1] for f,power in factors)
result={'status':'PASS','projection_omitted_factors_are_units':removed,'W_recurrence_to_W3':True,'Res_F_W3_exact_match':True,'H03_factorization_complete_up_to_nonzero_constant':str(quo),'radial_E_G27_exact_derivative':True,'pivot_A_substitution_exact':True,'pivot_removed_factors_are_physical_units':[[str(f),int(e)] for f,e in factors],'seconds':time.monotonic()-start}
(O/'independent_bridge.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
