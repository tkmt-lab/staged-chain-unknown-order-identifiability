#!/usr/bin/env python3
"""Exact dyadic interval checks for realization of certified root8."""
from fractions import Fraction
import json
from certify_interval import I,Q
BASE='./'
def interval(pair):
    return I.raw(I(Fraction(pair[0])).a,I(Fraction(pair[1])).b)
def endpoint(x,side):return I.raw(getattr(x,side),getattr(x,side))
def matmul(A,B):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)]for row in A]
def solve(A,b):
    M=[row[:]+[x] for row,x in zip(A,b)];pivots=[]
    for k in range(len(b)):
        pivot=max(range(k,len(b)),key=lambda j:abs((M[j][k].a+M[j][k].b)//2))
        M[k],M[pivot]=M[pivot],M[k];p=M[k][k];assert not p.a<=0<=p.b;pivots.append(p)
        M[k]=[x/p for x in M[k]]
        for j in range(len(b)):
            if j!=k:
                p=M[j][k];M[j]=[x-p*y for x,y in zip(M[j],M[k])]
    return [row[-1]for row in M],pivots

data=json.load(open(BASE+'lower_candidate.json'));cert=json.load(open(BASE+'lower_root_certificate.json'))
z=[interval(b)for b in cert['root_box_exact']];i,g,p2,p1,p0,nb,nq1,nq0=z
# Every root is defined as the unique root in this rational sign bracket.
radius=I('1e-28').b;rates=[];signs=[]
for j in range(3):
    center=I(-Fraction(data['right']['K'][j][j])/25).mid()
    r=I.raw(center.a-radius,center.b+radius)
    def f(t):return t*t*t-p2*t*t+p1*t-p0
    fl,fr=f(endpoint(r,'a')),f(endpoint(r,'b'))
    assert (fl.b<0 and fr.a>0) or (fl.a>0 and fr.b<0)
    assert r.a>0
    rates.append(r);signs.append([fl.bounds(),fr.bounds()])
assert all(rates[j].b<rates[j+1].a for j in range(2))
a,b,c=rates;u=a/2;v=b/2;q1=nq1/nb;q0=nq0/nb;beta=-25*nb
C=[I(1),(q1-b-c)/u,(c*c-q1*c+q0)/(u*v)]
K=[[-25*a,I(0),I(0)],[25*u,-25*b,I(0)],[I(0),25*v,-25*c]]
CK=matmul([C],K)[0];CK2=matmul([CK],K)[0]
O=[C,CK,CK2]
y0=i;y1=i*(1-g-i);y2=y1*(1-g-i)+i*(-i-(1-g)*i)
sright=I(Fraction(141,200));v0=i/sright;v1=(y1+v0*i)/sright;v2=(y2+2*v1*i+v0*y1)/sright
w0=v0/beta;w1=v1/beta;w2=v2/beta
X,pivots=solve(O,[w0,w1-y0,w2-y1-CK[0]*y0])
Rleft=100-1-i;Rright=100-sright-sum(X)
controllability=(25*u)**2*(25*v)
for value in [i,g,beta,*C,*X,Rleft,Rright,controllability]:assert value.a>0
column_sums=[sum(K[j][k] for j in range(3))for k in range(3)]
assert all(x.b<0 for x in column_sums)
F8=interval(cert['residual8_exact']);gap=-sright*(25**11)*40320*F8;assert gap.a>0
result={'method':'exact integer dyadic intervals, 384 bits; cubic sign brackets and interval Gaussian elimination','rates_fixed_T_exact':[r.bounds()for r in rates],'cubic_endpoint_signs_exact':signs,'betaR_exact':beta.bounds(),'C_exact':[x.bounds()for x in C],'X0_exact':[x.bounds()for x in X],'observable_elimination_pivots_exact':[x.bounds()for x in pivots],'controllability_determinant_exact':controllability.bounds(),'column_sums_exact':[x.bounds()for x in column_sums],'Rleft_exact':Rleft.bounds(),'Rright_exact':Rright.bounds(),'gap_order11_exact':gap.bounds(),'approx':{'rates_physical':[(25*r).approx()for r in rates],'betaR':beta.approx(),'C':[x.approx()for x in C],'X0':[x.approx()for x in X],'Rleft':Rleft.approx(),'Rright':Rright.approx(),'gap_order11':gap.approx()},'all_checks_passed':True}
with open(BASE+'lower_realization_certificate.json','w')as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result['approx'],indent=2));print('All exact checks passed.')
