import sympy as s,json,argparse
parser=argparse.ArgumentParser();parser.add_argument("--explore",action="store_true",help="Display non-certifying numerical root approximations");args=parser.parse_args()
from pathlib import Path
O=Path(__file__).resolve().parent
a,b=s.symbols('a b'); data=json.loads((O/'derived.json').read_text());P=s.sympify(data['P']);Q=s.sympify(data['Q']);res=s.sympify(data['resultant_b'])
sr=s.subresultants(P,Q,b);linear=sr[-2];print('LINEAR',s.factor(linear),flush=True)
Bn=s.Poly(linear,b).coeff_monomial(1);Bd=s.Poly(linear,b).coeff_monomial(b);bb=s.factor(-Bn/Bd);print('b=',bb,flush=True)
for fac,mult in s.factor_list(res)[1]:
 print('FACTOR',fac,'mult',mult,'gcd denom',s.gcd(fac,Bd),flush=True)
 if s.degree(fac,a)==1:continue
 if s.gcd(fac,Bd)!=1:
  print('EXCEPTION: linear subresultant vanishes; classify separately in certify.py.',flush=True)
  continue
 if not args.explore:continue
 print('roots',s.nroots(fac,maxsteps=200),flush=True)
 for root in s.nroots(fac,maxsteps=200):
  if abs(s.im(root))>1e-10:continue
  bv=bb.subs(a,root);vals={key:s.N(s.sympify(data[key]).subs(b,bv).subs(a,root),9) for key in ['h','k','g','r','d1','d0','v']}
  print('ROOT',root,'b',bv,'vals',vals,flush=True)
data.update(b_formula=str(bb),linear_subresultant=str(s.factor(linear)))
(O/'classified.json').write_text(json.dumps(data,indent=2)+'\n')
