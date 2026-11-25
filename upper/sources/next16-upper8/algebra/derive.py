import sympy as s, json, time
from pathlib import Path
O=Path(__file__).resolve().parent
S,L,k,h,r,d0,d1,g,v,t,u,w,a,b=s.symbols('S L k h r d0 d1 g v t u w a b')
X=r+1-S-k*L
Dx=lambda f:s.diff(f,L)+S*s.diff(f,S)
Dt=lambda f:-X*Dx(f)
y=S*X;z=y/(S+h)
F=s.Poly(s.cancel((S+h)**3/y*(Dt(Dt(z))+d1*Dt(z)+d0*z-g*Dt(y)-v*y)),S,L)
W=L**2/s.Integer(288)+s.Rational(59,1728)*L+s.Rational(671,6912)+S*(-L**2/s.Integer(24)-L/s.Integer(36)-s.Rational(103,432))+S**2*(-L/s.Integer(16)+s.Rational(1,8))+S**3*(-L/s.Integer(216)+s.Rational(7,432))+S**4/s.Integer(6912)
sol=s.solve([1+t+u+w-a,1+2*t+4*u+8*w-b,1+3*t+9*u+27*w-54*a*b],[t,u,w]);print('sol',sol,flush=True)
E=216*a*b-6*b+4*a-1
hh=1/(12*a);kk=-3456*a*b/E;ll=288*kk**2*hh**2;gg=s.factor(1-ll*E/13824)
G=s.Poly(F.as_expr()-288*k*k*h*h*(W+t*Dx(W)+u*Dx(Dx(W))+w*Dx(Dx(Dx(W)))),S,L)
subs={h:hh,k:kk,g:gg,**sol}
coef={m:s.factor(q.subs(subs)) for m,q in G.terms()}
for m in [(0,2),(1,2),(4,0),(3,1),(2,1)]:assert coef.get(m,0)==0
M,B=s.linear_eq_to_matrix([coef[(1,1)],coef[(0,1)]],[r,d1]);det=s.factor(M.det());print('det',det,flush=True)
rr,dd=map(s.factor,M.inv()*B);print('r',rr,'d1',dd,flush=True)
after={m:s.factor(q.subs({r:rr,d1:dd})) for m,q in coef.items()}
vv=s.factor(s.solve(after[(3,0)],v)[0]);ddd=s.factor(s.solve(after[(2,0)].subs(v,vv),d0)[0]);print('v',vv,'d0',ddd,flush=True)
assert s.factor(s.diff(after[(3,0)],v)+1)==0
assert s.factor(s.diff(after[(2,0)],d0)-1)==0
end={m:s.factor(q.subs({v:vv,d0:ddd})) for m,q in after.items()};print('END',end,flush=True)
assert all(z==0 for m,z in end.items() if m not in [(1,0),(0,0)])
out={'tuw':{str(x):str(y) for x,y in sol.items()},'h':str(hh),'k':str(kk),'g':str(gg),'det':str(det),'r':str(rr),'d1':str(dd),'v':str(vv),'d0':str(ddd),'end':{str(x):str(y) for x,y in end.items()}}
(O/'derived.json').write_text(json.dumps(out,indent=2)+'\n')
p,q=[s.Poly(s.fraction(end[m])[0],a,b).primitive()[1].as_expr() for m in [(1,0),(0,0)]]
p=s.cancel(p/b**2);q=s.cancel(q/b**2)
print('P',p,'Q',q,flush=True)
res=s.factor(s.resultant(p,q,b));print('RESULTANT',res,flush=True)
out.update(P=str(p),Q=str(q),resultant_b=str(res));(O/'derived.json').write_text(json.dumps(out,indent=2)+'\n')
