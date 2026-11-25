"""Exact Sturm counts and rational interval certificates; no floating decisions."""
import sympy as s,json
from fractions import Fraction as Q
from pathlib import Path
a,b=s.symbols('a b'); j=json.loads(Path(__file__).with_name('general.json').read_text()); E=s.sympify(j['E']);R=list(map(s.sympify,j['residuals']))
f=s.cancel(R[0]*a**3*E**2/(8*b*b));g=s.cancel(R[1]*6*a**4*E**2/(b*b))
res=s.factor(s.resultant(f,g,b));factors=s.factor_list(res)[1];C=next(p for p,n in factors if s.degree(p,a)==3);H=next(p for p,n in factors if s.degree(p,a)==12)
# No division is used on the leading-coefficient exceptional factor.
fc=s.rem(f,C,a);gc=s.rem(g,C,a);assert s.degree(fc,b)==s.degree(gc,b)==1
Dc=s.factor(s.resultant(fc,gc,b));u,v,d=s.gcdex(C,Dc,a);assert d==1 and s.expand(u*C+v*Dc)==1
# Unique b on H=0 from the penultimate subresultant.
lin=s.Poly(s.subresultants(f,g,b)[-2],b)
bval=s.factor(-lin.nth(0)/lin.nth(1)); bn,bd=s.fraction(bval)
assert s.gcd(H,bd)==1 and s.gcd(H,C)==1
assert s.rem(s.fraction(s.cancel(f.subs(b,bval)))[0],H,a)==0
assert s.rem(s.fraction(s.cancel(g.subs(b,bval)))[0],H,a)==0
# Count all real roots using signs at +/- infinity of a rational Sturm chain.
st=s.sturm(H,a)
def variations(signs):
 return sum(x!=y for x,y in zip([x for x in signs if x],[x for x in signs if x][1:]))
def signs_inf(right):
 return [int(s.sign(s.LC(s.Poly(p,a))))*int((-1)**s.degree(p,a) if not right else 1) for p in st]
def var_at(x): return variations([int(s.sign(p.subs(a,x))) for p in st])
realcount=variations(signs_inf(False))-variations(signs_inf(True));assert realcount==4
intervals=[('-0.777954023820','-0.777954023819','g'),('-0.118906417694','-0.118906417693','d1'),('-0.114192896216','-0.114192896215','d1'),('0.025568747572','0.025568747573','r')]
def add(I,J):return (I[0]+J[0],I[1]+J[1])
def mul(I,J):
 p=[x*y for x in I for y in J];return min(p),max(p)
def poly_interval(p,I):
 ans=(Q(0),Q(0))
 for c in s.Poly(p,a).all_coeffs():
  cq=Q(int(s.numer(c)),int(s.denom(c)));ans=add(mul(ans,I),(cq,cq))
 return ans
out={'status':'PASS','resultant':str(res),'C':str(C),'H':str(H),'C_remainders':[str(fc),str(gc)],'C_linear_determinant':str(Dc),'C_bezout':[str(u),str(v)],'b':str(bval),'real_root_count':realcount,'sturm_signs_minus_infinity':signs_inf(False),'sturm_signs_plus_infinity':signs_inf(True),'root_certificates':[]}
sol={**j['solutions'],**j['basic']}
for lo,hi,key in intervals:
 l,u=s.Rational(lo),s.Rational(hi);assert var_at(l)-var_at(u)==1
 l,u=s.refine_root(H,l,u,eps=s.Rational(1,10**35));assert var_at(l)-var_at(u)==1
 interval=(Q(int(s.numer(l)),int(s.denom(l))),Q(int(s.numer(u)),int(s.denom(u))))
 val=s.factor(s.sympify(sol[key]).subs(b,bval));num,den=s.fraction(val)
 num=s.rem(num,H,a);den=s.rem(den,H,a)
 assert s.gcd(H,den)==1
 ni=poly_interval(num,interval);di=poly_interval(den,interval)
 assert ni[0]*ni[1]>0 and di[0]*di[1]>0 and ni[0]*di[0]<0
 bounds=[x/y for x in ni for y in di];bounds=(min(bounds),max(bounds));assert bounds[1]<0
 out['root_certificates'].append({'a_interval':[str(l),str(u)],'sturm_variations':[var_at(l),var_at(u)],'violated_positive_quantity':key,'reduced_numerator':str(num),'reduced_denominator':str(den),'numerator_interval':list(map(str,ni)),'denominator_interval':list(map(str,di)),'value_interval':list(map(str,bounds))})
 print('root',lo,hi,'reject',key,'negative; exact interval certificate saved',flush=True)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS: all four H real candidates fail necessary positivity; C has no finite candidate.')
