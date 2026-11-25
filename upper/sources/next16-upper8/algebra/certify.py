"""Exact QQ elimination and Sturm certificate; no float enters any verdict."""
import sympy as s,json,platform,time,resource
from pathlib import Path
O=Path(__file__).resolve().parent;a,b=s.symbols('a b');D=json.loads((O/'classified.json').read_text());P=s.sympify(D['P']);Q=s.sympify(D['Q']);res=s.sympify(D['resultant_b'])
C=432*a**3-252*a**2+24*a-1
H=s.Poly(s.div(res,2916*(12*a+1)**7*C,a)[0],a)
assert s.expand(res-2916*(12*a+1)**7*C*H.as_expr())==0
assert s.factor(s.resultant(P,Q,b)-res)==0
cp,cq=[s.rem(z,C,a) for z in (P,Q)]
assert s.degree(cp,b)==s.degree(cq,b)==1
cr=s.resultant(cp,cq,b);U,V,one=s.gcdex(cr,C,a);assert one==1 and s.expand(U*cr+V*C)==1
bb=s.sympify(D['b_formula']);lin=s.sympify(D['linear_subresultant']);assert s.factor(lin-s.subresultants(P,Q,b)[-2])==0
bn,bd=s.fraction(bb);assert s.gcd(bd,H.as_expr())==1
# The stated b is forced by the subresultant on all H roots.
assert s.factor(lin.subs(b,bb))==0
assert s.gcd(C,H.as_expr())==1
assert s.rem(s.fraction(s.cancel(P.subs(b,bb)))[0],H.as_expr(),a)==0
assert s.rem(s.fraction(s.cancel(Q.subs(b,bb)))[0],H.as_expr(),a)==0
assert H.count_roots(-s.oo,s.oo)==4
intervals=[(-777954024,-777954023,'g'),(-118906418,-118906417,'d1'),(-114192897,-114192896,'d1'),(25568747,25568748,'r')]
rows=[]
for lo,hi,key in intervals:
 lo=s.Rational(lo,10**9);hi=s.Rational(hi,10**9)
 assert H.count_roots(lo,hi)==1
 expr=s.cancel(s.sympify(D[key]).subs(b,bb));num,den=s.fraction(expr)
 # Signs are certified on the original canceled numerator and denominator.
 num=s.Poly(num,a);den=s.Poly(den,a)
 assert num.count_roots(lo,hi)==0 and den.count_roots(lo,hi)==0
 ns=s.sign(num.eval(lo));ds=s.sign(den.eval(lo));assert ns*ds==-1
 rows.append({'a_interval':[str(lo),str(hi)],'H_root_count':1,'negative_parameter':key,'numerator':str(num.as_expr()),'denominator':str(den.as_expr()),'num_interval_root_count':0,'den_interval_root_count':0,'num_left_sign':int(ns),'den_left_sign':int(ds)})
# Negative controls: a false nonnegative verdict and a damaged resultant fail.
assert all(row['num_left_sign']*row['den_left_sign']!=1 for row in rows)
assert s.expand((res+1)-s.resultant(P,Q,b))!=0
out={'status':'PASS','input_commit':'4b3486ad1f4750833c97411f21d391b80ab11d5f','python':platform.python_version(),'sympy':s.__version__,'H':str(H.as_expr()),'C':str(C),'resultant':str(res),'C_P_remainder':str(s.factor(cp)),'C_Q_remainder':str(s.factor(cq)),'C_remainder_resultant':str(s.factor(cr)),'C_bezout_U':str(U),'C_bezout_V':str(V),'C_bezout_identity':'U*Res_b(P mod C,Q mod C)+V*C=1','b_formula':str(bb),'H_total_real_root_count':4,'roots':rows,'max_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'conclusion':'All general-branch real contact candidates fail g>0, d1>0, or r>0.'}
(O/'CERTIFICATE.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':out['status'],'roots':[{k:v for k,v in z.items() if k not in ('numerator','denominator')} for z in rows]},indent=2))
