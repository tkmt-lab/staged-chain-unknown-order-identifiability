"""Exceptional cubic factor: reject finite b points before division."""
import sympy as s,json
from pathlib import Path
O=Path(__file__).resolve().parent;a,b=s.symbols('a b');D=json.loads((O/'derived.json').read_text());P=s.sympify(D['P']);Q=s.sympify(D['Q']);C=432*a**3-252*a**2+24*a-1
pr=s.rem(P,C,a);qr=s.rem(Q,C,a)
assert s.degree(pr,b)==s.degree(qr,b)==1
R=s.resultant(pr,qr,b);U,V,one=s.gcdex(R,C,a)
assert one==1 and s.expand(U*R+V*C)==1
print('P mod C =',s.factor(pr));print('Q mod C =',s.factor(qr))
print('linear resultant =',s.factor(R));print('gcd with C =',one)
print('PASS: no finite common b on cubic branch.')
