import sympy as s,json
from pathlib import Path
S,X,L,k,h,r,d0,d1,g,v=s.symbols('S X L k h r d0 d1 g v')
Dt=lambda f:-S*X*s.diff(f,S)+(k+S)*X*s.diff(f,X)
y=S*X;z=y/(S+h)
F=s.Poly(s.expand(s.cancel((S+h)**3/y*(Dt(Dt(z))+d1*Dt(z)+d0*z-g*Dt(y)-v*y)).subs(X,r+1-S-k*L)),S,L)
Dx=lambda f:s.expand(S*s.diff(f,S)+s.diff(f,L))
mons=[S**i*L**j for i,j in F.monoms()]
rows=[];vs=mons[:]
for n in range(11):
 rows.append([e.subs({S:1,L:0}) for e in vs]);vs=[Dx(e) for e in vs]
M=s.Matrix(rows);wv=M.inv()*s.Matrix([0]*10+[1]);W=sum(e*c for e,c in zip(mons,wv))
# P(z)=-(z-2)(z-3)(z-4)/24=1-13z/12+3z²/8-z³/24.
V=W-s.Rational(13,12)*Dx(W)+s.Rational(3,8)*Dx(Dx(W))-s.Rational(1,24)*Dx(Dx(Dx(W)))
G=s.Poly(s.expand(F.as_expr().subs({h:s.Rational(1,3),g:1})-32*k*k*V),S,L)
assert all(G.coeff_monomial(m)==0 for m in [L*L,S*L*L,S**4,S**3*L,S**2*L])
mat,rhs=s.linear_eq_to_matrix([G.coeff_monomial(m) for m in [S*L,L]],[r,d1]);det=s.factor(mat.det());rd=[s.factor(e) for e in mat.inv()*rhs]
vv=s.solve(G.coeff_monomial(S**3).subs(dict(zip([r,d1],rd))),v)[0]
dd=s.solve(G.coeff_monomial(S*S).subs(dict(zip([r,d1],rd))).subs(v,vv),d0)[0]
res=[s.factor(G.coeff_monomial(m).subs(dict(zip([r,d1],rd))).subs({v:vv,d0:dd})) for m in [S,s.Integer(1)]]
print('det',det,'r,d1',rd,'v',vv,'d0',dd,'residuals',res)
assert det==4*k*k/27
assert res==[3791*k*k/972,-10013*k*k/8748]
nums=[s.fraction(e)[0] for e in res]
print('gcd',s.factor(s.gcd(*nums)))
Path(__file__).with_suffix('.json').write_text(json.dumps({'det':str(det),'r':str(rd[0]),'d1':str(rd[1]),'v':str(vv),'d0':str(dd),'residuals':list(map(str,res)),'gcd':str(s.factor(s.gcd(*nums)))},indent=2)+'\n')
