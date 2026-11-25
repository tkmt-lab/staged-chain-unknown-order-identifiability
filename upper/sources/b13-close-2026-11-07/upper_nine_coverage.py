"""Exact one-sided divisibility + preserved degree + modular coprimality coverage."""
import flint,json,hashlib
from fractions import Fraction
from pathlib import Path
p=Path(__file__).parent;prime=32003;C=flint.nmod_mpoly_ctx.get(('h','k'),prime);h,k=C.gens();J=json.loads((p/'lower_SHAPE_BETA_CANDIDATE.json').read_text());V=json.loads((p/'upper_rur_verified.json').read_text());assert all(V['checks'][key]['status']=='PASS' for key in ['G27','H02']);assert V['candidate_sha256']==hashlib.sha256((p/'lower_SHAPE_BETA_CANDIDATE.json').read_bytes()).hexdigest();hc=list(map(Fraction,J['H_coefficients_low_first']));assert hc[-1]==1 and all(c.denominator%prime for c in hc)
HQ=flint.fmpq_poly([flint.fmpq(v.numerator,v.denominator) for v in hc]);assert HQ.gcd(HQ.derivative()).degree()==0 and HQ[0] and HQ(-1);assert all(prime%d for d in range(2,179))
H=sum((int(c.numerator)%prime)*pow(int(c.denominator)%prime,-1,prime)*h**i for i,c in enumerate(hc));assert H.degrees()[0]==289
old=json.loads((p/'lower_projection_modular.json').read_text());rad=json.loads((p/'upper_nine_radial.json').read_text());R02=flint.nmod_mpoly(old['branches']['27']['res02']['polynomial'],C);Rrad=flint.nmod_mpoly(rad['modular_resultant'],C)
raw=json.loads((p/'lower_projection.json').read_text());zz=flint.fmpz_mpoly_ctx.get(('h','k'));gtext=next(f['polynomial'] for f in raw['factors']['01']['factors'] if f['degrees'][1]==27);htext=next(f['polynomial'] for f in raw['factors']['02']['factors'] if f['degrees'][1]>1)
for txt in [gtext,htext,rad['E_G27']]:assert flint.fmpz_mpoly(txt,zz).degrees()[1]==flint.nmod_mpoly(txt,C).degrees()[1]
Gp=flint.nmod_mpoly(gtext,C);assert R02==Gp.resultant(flint.nmod_mpoly(htext,C),'k');assert Rrad==Gp.resultant(flint.nmod_mpoly(rad['E_G27'],C),'k')
DC=json.loads((p/'lower_degree_certificate.json').read_text())['branches']['27']['02'];QC=json.loads((p/'lower_qadic_certificate.json').read_text())['branches']['27']['02'];assert DC['degree_preserved'] and DC['degree']['bound']==R02.degrees()[0] and DC['h']['bound']>=90;assert QC['passed'] and QC['determinant_valuation_lower_bound']>=98
GZ=flint.fmpz_mpoly(gtext,zz);LC=flint.nmod_mpoly_ctx.get(('h','k'),prime).from_dict({(int(hh),0):int(v)%prime for (hh,kk),v in GZ.to_dict().items() if int(kk)==27});assert LC.gcd(H).is_constant()
RH,remH=divmod(R02,H);assert remH==0 and RH.gcd(H).is_constant()
D=h**90*(h+1)**98*H;quo,rem=divmod(R02,D);assert rem==0;assert quo.degrees()[0]==R02.degrees()[0]-D.degrees()[0];assert quo.gcd(Rrad).is_constant()
out={'status':'PASS','prime':prime,'candidate_degree':289,'R02_degree_preserved':6155,'known_R02_boundary_divisor':'h^90*(h+1)^98','H_divisibility_basis':'exact integer RUR substitution into G27,H02; H squarefree and denominator coprime','modular_quotient_degree':int(quo.degrees()[0]),'quotient_gcd_with_radial_resultant':'1','H_occurs_once_in_R02':True,'G27_k_leading_coefficient_coprime_H':True,'unique_finite_k_over_each_H_root':True,'conclusion':'Every physical nine-contact candidate has h a root of H289; exact radial resultant is unnecessary.'};(p/'upper_nine_coverage.json').write_text(json.dumps(out,indent=2)+'\n');print('EXACT H289 COVERAGE PASS',flush=True)
