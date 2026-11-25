"""Check the last modular-to-QQ gcd step and its exact certificate prerequisites."""
import flint,json
from pathlib import Path
p=Path(__file__).parent;prime=32003;c=flint.nmod_mpoly_ctx.get(('h','k'),prime);h,k=c.gens();DC=json.loads((p/'lower_degree_certificate.json').read_text());QC=json.loads((p/'lower_qadic_certificate.json').read_text());old=json.loads((p/'lower_projection_modular.json').read_text());new=json.loads((p/'upper_contact10_modular.json').read_text());out={'status':'PASS','prime':prime,'branches':{}}
for tag,a,b in [('27',90,98),('22',62,92)]:
 extra=c.constant(1) if tag=='27' else (3*h+5)**2*(3*h*h+8*h-4)*(14*h*h+40*h+29)
 D=h**a*(h+1)**b*extra;divs=[]
 for key in ['02','03']:
  r=flint.nmod_mpoly(old['branches'][tag]['res02']['polynomial'] if key=='02' else new['branches'][tag][key]['resultant'],c)
  B=DC['branches'][tag][key];assert B['degree_preserved'] and B['degree']['bound']==int(r.degrees()[0]);assert B['h']['bound']>=a
  cert=QC['branches'][tag][key];assert cert['passed'] and cert['determinant_valuation_lower_bound']>=b
  quo,rem=divmod(r,D);assert rem==0;assert quo.degrees()[0]==r.degrees()[0]-D.degrees()[0];divs.append(quo)
 assert divs[0].gcd(divs[1]).is_constant()
 out['branches'][tag]={'known_divisor_mod_p':str(D),'quotient_gcd':'1','quotient_degrees':[int(f.degrees()[0]) for f in divs]}
 print(tag,'quotient gcd 1, leading degrees preserved PASS',flush=True)
# Other exact divisor certificates are included by named files; their individual
# computation is not silently replaced by the modular check.
S=json.loads((p/'contact10_specialize.json').read_text());assert all(S[H][key]['gcd_degree']>=1 for H in ['3*h^2+8*h-4','14*h^2+40*h+29'] for key in ['02','03'])
T=json.loads((p/'contact10_tangent.json').read_text());assert T['status']=='PASS' and all(T['pairs'][key]['resultant_valuation_lower_bound']>=2 and T['pairs'][key]['jacobian_determinant']=='0' for key in ['02','03']);out['additional_exact_divisor_certificates']=['contact10_specialize.json','contact10_tangent.json'];out['note']='The tangent certificate supplies the second power of 3h+5; see upper_contact10_PROOF.md.'
(p/'upper_contact10_lift.json').write_text(json.dumps(out,indent=2)+'\n')
