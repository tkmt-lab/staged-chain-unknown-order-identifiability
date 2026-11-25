"""Exact certificate that the t-linear subresultant LC is nonzero on physical
G27 nine-contact candidates, via degree-preserved modular coprimality.
"""
from pathlib import Path
import json,hashlib
from flint import fmpz_mpoly,nmod_mpoly
O=Path(__file__).resolve().parent
D={'__file__':str(O/'lower_degree_certificate.py')};source=(O/'lower_degree_certificate.py').read_text();exec(source[:source.index("out={'status':")],D)
Q={'__file__':str(O/'lower_qadic_certificate.py')};source=(O/'lower_qadic_certificate.py').read_text();exec(source[:source.index("out={'scope':")],Q)
raw=json.loads((O/'lower_projection.json').read_text());lin=json.loads((O/'upper_nine_linear.json').read_text());old=json.loads((O/'lower_projection_modular.json').read_text());base=json.loads((O/'lower_degree_certificate.json').read_text());qbase=json.loads((O/'lower_qadic_certificate.json').read_text())
Z=D['Z'];M=D['M'];g=fmpz_mpoly(next(f['polynomial'] for f in raw['factors']['01']['factors'] if f['degrees'][1]==27),Z);a=fmpz_mpoly(lin['a'],Z);h,k=M.gens()
bounds=D['bound'](g,a);gm=nmod_mpoly(str(g),M);am=nmod_mpoly(str(a),M);assert gm.degrees()[1]==g.degrees()[1] and am.degrees()[1]==a.degrees()[1]
r=gm.resultant(am,'k');assert bounds['h']['bound']>=90 # Degree preservation of R02 alone forces every common primitive factor to survive modulo p.
print('DEGREE/H',bounds['degree']['bound'],bounds['h']['bound'],flush=True)
qcert=Q['certificate'](g,a,12,98);assert qcert['passed']
r02=nmod_mpoly(old['branches']['27']['res02']['polynomial'],M);common=r.gcd(r02);unit=h**90*(h+1)**98;assert common==unit
assert (r//unit).gcd(r02//unit)==M.constant(1)
assert base['branches']['27']['02']['degree_preserved'] and base['branches']['27']['02']['input_k_degrees_preserved']
assert base['branches']['27']['02']['h']['bound']>=90 and qbase['branches']['27']['02']['determinant_valuation_lower_bound']>=98
paths=['lower_projection.json','upper_nine_linear.json','lower_projection_modular.json','lower_degree_certificate.json','lower_qadic_certificate.json']
out={'status':'PASS','claim':'For h*k*(h+1)!=0, a=0 is impossible together with G27=0 and H02=0; hence a!=0 on every G27 nine-contact candidate.','prime':32003,'input_k_degrees':list(map(int,[g.degrees()[1],a.degrees()[1]])),'input_k_degrees_preserved':True,'new_resultant_mod_degree':int(r.degrees()[0]),'new_resultant_exact_degree_not_asserted':True,'degree_preserved_reference':'R02, sufficient to preserve every potential common factor','degree_and_h_certificate':bounds,'q_certificate':qcert,'mod_gcd':'h^90*(h+1)^98','quotient_mod_gcd':'1','dependencies':{p:hashlib.sha256((O/p).read_bytes()).hexdigest() for p in paths}}
(O/'lower_LINEAR_CERTIFICATE.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS a is nonzero on G27 nine-contact candidates')
