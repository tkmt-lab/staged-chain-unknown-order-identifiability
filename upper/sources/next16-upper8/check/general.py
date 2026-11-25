"""Independent physical-coordinate differentiation and jet interpolation, QQ(a,b)."""
import sympy as s,json
from pathlib import Path
S,X,L,k,h,r,d0,d1,g,v,a,b=s.symbols('S X L k h r d0 d1 g v a b')
Dt=lambda f:-S*X*s.diff(f,S)+(k+S)*X*s.diff(f,X)
y=S*X;z=y/(S+h)
F=s.Poly(s.expand(s.cancel((S+h)**3/y*(Dt(Dt(z))+d1*Dt(z)+d0*z-g*Dt(y)-v*y)).subs(X,r+1-S-k*L)),S,L)
assert set(F.monoms())=={(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1),(3,0),(3,1),(4,0)}
assert set(s.Poly(F.as_expr().subs(h,0),S,L).monoms())=={(2,0),(3,0),(3,1),(4,0)}
Dx=lambda f:s.expand(S*s.diff(f,S)+s.diff(f,L))
mons=[S**i*L**j for i,j in F.monoms()]
rows=[];vs=mons[:]
for n in range(11):
 rows.append([e.subs({S:1,L:0}) for e in vs]);vs=[Dx(e) for e in vs]
M=s.Matrix(rows);assert M.det()!=0 and s.Matrix(rows[:7]).rank()==7
wv=M.inv()*s.Matrix([0]*10+[1]);W=sum(e*c for e,c in zip(mons,wv))
t,u,w=s.symbols('t u w'); pars=s.solve([1+t+u+w-a,1+2*t+4*u+8*w-b,1+3*t+9*u+27*w-54*a*b],[t,u,w])
E=(216*a-6)*b+4*a-1
vals={h:1/(12*a),k:-3456*a*b/E,g:1-1728*b*b/E}
lam=s.factor((288*k*k*h*h).subs(vals));V=(W+t*Dx(W)+u*Dx(Dx(W))+w*Dx(Dx(Dx(W)))).subs(pars)
# Reduce each coefficient before expanding the whole large rational function.
G={mon:s.factor(F.coeff_monomial(mon).subs(vals)-lam*s.Poly(V,S,L).coeff_monomial(mon)) for mon in mons}
assert all(G[m]==0 for m in [L*L,S*L*L,S**4,S**3*L,S**2*L])
mat,rhs=s.linear_eq_to_matrix([G[m] for m in [S*L,L]],[r,d1]);det=s.factor(mat.det());rd=[s.factor(e) for e in mat.inv()*rhs]
assert s.factor(det-27648*b*b/(a*E*E))==0
assert s.cancel(s.diff(G[S**3],v)+1)==0 and s.cancel(s.diff(G[S*S],d0)-1)==0
vv=s.factor(-G[S**3].subs(dict(zip([r,d1],rd))).subs(v,0)/s.diff(G[S**3],v))
dd=s.factor(-G[S*S].subs(dict(zip([r,d1],rd))).subs(v,vv).subs(d0,0)/s.diff(G[S*S],d0))
res=[s.factor(G[m].subs(dict(zip([r,d1],rd))).subs({v:vv,d0:dd})) for m in [S,s.Integer(1)]]
for mon in mons:
 if mon not in [S,s.Integer(1)]:
  assert s.cancel(G[mon].subs(dict(zip([r,d1],rd))).subs({v:vv,d0:dd}))==0
out={'status':'PASS','verified_support_size':11,'verified_h0_support_size':4,'zero_coefficients_after_substitution':9,'parameters':{str(k):str(e) for k,e in pars.items()},'E':str(E),'det':str(det),'solutions':{str(k):str(e) for k,e in zip([r,d1,v,d0],rd+[vv,dd])},'residuals':list(map(str,res)),'basic':{str(k):str(e) for k,e in vals.items()}}
print(json.dumps(out,indent=2),flush=True)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
