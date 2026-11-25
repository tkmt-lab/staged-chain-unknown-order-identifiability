import sympy as s,json
from pathlib import Path
a,b=s.symbols('a b');j=json.loads(Path(__file__).with_name('general.json').read_text()); E=s.sympify(j['E']);R=[s.sympify(e) for e in j['residuals']]
f=s.Poly(s.cancel(R[0]*a**3*E**2/(8*b*b)),b).as_expr();g=s.Poly(s.cancel(R[1]*6*a**4*E**2/(b*b)),b).as_expr()
res=s.factor(s.resultant(f,g,b)); print('resultant',res,flush=True)
sub=s.subresultants(f,g,b); print('linear',s.factor(sub[-2]),flush=True)
L=s.Poly(sub[-2],b); bb=s.factor(-L.nth(0)/L.nth(1)); print('b',bb,flush=True)
fac=s.factor_list(res)[1];polys=[p for p,n in fac if s.degree(p,a)>1]
sol={k:s.sympify(e).subs(b,bb) for k,e in {**j['solutions'],**j['basic']}.items()}
out={'status':'EXPLORATION_ONLY_NOT_A_CERTIFICATE','warning':'C3 factor shown below is an extraneous resultant factor; substituting its root into the generic b formula does not produce a candidate. See classify_exact.json for exact exclusion.','resultant':str(res),'linear':str(s.factor(sub[-2])),'b':str(bb),'f':str(f),'g':str(g),'roots':[]}
for p in polys:
 print('degree',s.degree(p,a),flush=True)
 intervals=s.polys.polytools.intervals(p,eps=s.Rational(1,10**20))
 for iv,m in intervals:
  av=sum(iv)/2; bv=bb.subs(a,av); values={k:float(e.subs(a,av)) for k,e in sol.items()}
  print('root',tuple(map(str,iv)),'a',float(av),'b',float(bv),values,flush=True)
  out['roots'].append({'interval':list(map(str,iv)),'multiplicity':m,'a_approx':float(av),'b_approx':float(bv),'values_approx':values})
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
