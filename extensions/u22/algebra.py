#!/usr/bin/env python3
"""Rebuild the degree-eight constraints by Newton interpolation.
Use only our universal.json for model coefficients. Producer certificates
are read only for coefficient comparison, never used to construct constraints.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import json,math,time,resource,hashlib,argparse
from pathlib import Path
resource.setrlimit(resource.RLIMIT_AS,(7*1024**3,7*1024**3))
resource.setrlimit(resource.RLIMIT_CPU,(7000,7000))
import sympy as s
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'b22data'
x,h,u,lam=s.symbols('x h u lam')
a,b,c,d,e,g,v,A=s.symbols('a b c d e g v A')
def eq(f,z=0):assert s.cancel(f-z)==0, str(f-z)[:300]
def primitive(f,*vv):
    poly=s.Poly(s.cancel(f),*vv).clear_denoms()[1].primitive()[1]
    return poly if poly.LC()>0 else -poly
def write(name,obj):
    obj.update(status='PASS',seconds=time.monotonic()-START,rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (HERE/(name+'.json')).write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n')
    print(name,'PASS',obj['seconds'],'seconds',flush=True)
def interpolate():
    maxima={int(k):z for k,z in json.loads((HERE/'universal.json').read_text())['support_bound'].items()}
    mult=[sum(i<=mx for mx in maxima.values()) for i in range(max(maxima.values())+1)]
    assert sum(mult)==82
    # Recover full partial fractions by local truncated series of 1/prod_{j!=i}(x-j)^m_j.
    chi=s.Poly(s.prod((x-i)**mi for i,mi in enumerate(mult)),x,domain=s.QQ)
    partial_numerator=s.Poly(0,x,domain=s.QQ)
    wi=[];second=[]
    for i,mi in enumerate(mult):
        coeff=[s.S.One]+[s.S.Zero]*(mi-1)
        for j,mj in enumerate(mult):
            if i==j:continue
            fac=[(-1)**n*s.binomial(mj+n-1,n)/s.Integer(i-j)**(mj+n) for n in range(mi)]
            coeff=[sum(coeff[k]*fac[n-k] for k in range(n+1)) for n in range(mi)]
        for n,cn in enumerate(coeff):
            divisor=s.Poly((x-i)**(mi-n),x)
            partial_numerator+=chi.exquo(divisor).mul_ground(cn)
        wi.append(coeff[0]/s.factorial(mi-1))
        second.append(coeff[1]/s.factorial(mi-2) if mi>1 else 0)
    assert partial_numerator==s.Poly(1,x,domain=s.QQ)
    rr=[w/wi[0] for w in wi]
    for i in range(5):
        sigma=sum(s.Rational(mult[j],i-j) for j in range(len(mult)) if i!=j)
        eq(second[i]/wi[0],-5*rr[i]*sigma)
    nodes=[0,1,2,3,4,8,12,16]
    values=[h**4,2*h**3,3*h*h,2*h,1,8*u,s.Rational(163,27)*u*u,-s.Rational(2,27)*u**3]
    # Newton divided differences, rather than Vandermonde inversion or Lagrange sum.
    work=[val/rr[node] for node,val in zip(nodes,values)]
    coeffs=[]
    while work:
        coeffs.append(work[0]); order=len(coeffs)
        work=[s.expand((work[k+1]-work[k])/(nodes[k+order]-nodes[k])) for k in range(len(work)-1)]
    p7=s.expand(sum(val*s.prod(x-nodes[j] for j in range(k)) for k,val in enumerate(coeffs)))
    kernel=s.prod(x-i for i in nodes)
    kernel=s.expand(kernel/(rr[15]*kernel.subs(x,15)))
    p=p7+lam*kernel
    for node,value in zip(nodes,values):eq(rr[node]*p.subs(x,node),value)
    eq(rr[15]*kernel.subs(x,15),1)
    vals={i:s.expand(rr[i]*p.subs(x,i)) for i in [6,7,10,11,14,15]}
    qs=[]
    for i in range(9):
        qs.append(s.expand(rr[i]*p.subs(x,i)) if mult[i]==5 else s.expand(5*rr[i]*s.diff(p,x).subs(x,i)+second[i]/wi[0]*p.subs(x,i)))
    return mult,rr,p7,p,vals,qs

def constraints():
    mult,rr,p7,p,vals,qs=interpolate()
    source=json.loads((HERE/'universal.json').read_text())
    edge=lambda i,j:s.sympify(source['edges'][f'{i},{j}'])
    t=g-1;delta=a-d;C0=108*h**7*b**5;u0=t*t/(h*h*b)
    N=[s.cancel(edge(i,4)/(36*b**4*h**(11-i))) for i in range(9)]
    eq(-N[1]+N[2]-N[4]+N[5],4*delta**2+(47*t-4)*delta*h+(1-36*t-21*t*t)*h*h)
    cc=[1,-2,1,1,-2,1,1,-14,s.Rational(215,8)]
    eq(sum(k*n for k,n in zip(cc,N)),150*h*h*t)
    actual={i:s.cancel(edge(i,j)/C0) for i,j in [(6,4),(7,4),(10,3),(11,3),(14,2),(15,2)]}
    eq(actual[7],u0/3*(7*delta/t+50*h-16*h/t))
    eq(actual[11],u0**2/27*(259*delta/t+804*h+76*h/t))
    eq(actual[15],-8*u0**3/27*(3*h-h/t))
    eq(54*u0*actual[11]-222*u0*u0*actual[7]-4509*actual[15]-1916*h*u0**3)
    D,Z,eta,xi=s.symbols('D Z eta xi')
    tri=[D*D-10*eta+D*(15*h-4*Z)-3*xi-26*h*h-18*h*Z,
      387*D*D-372*eta+D*(4136*h-528*Z)+612*xi+4080*h*h+704*h*Z-348*Z*Z,
      3*D*D-54*eta-D*(28*h+16*Z)+56*xi-252*h*h+92*h*Z+152*Z*Z]
    subs={D:delta/t,Z:h/t,eta:a*h/t,xi:h*(v-c)/t**2}
    for form,left in zip(tri,[3*actual[6]/u0,108*actual[10]/u0**2,54*actual[14]/u0**3]):eq(form.subs(subs),left)
    combo=s.expand(6108*tri[0]+361*tri[1]-3618*tri[2])
    assert not combo.has(eta,xi)
    expected=134961*D*D+(1686020*h-157152*Z)*D+2225808*h*h-188656*h*Z-675564*Z*Z
    eq(combo,expected)
    assert s.cancel(sum(k*n for k,n in zip(cc,N))+1-150*h*h*t)!=0
    # u=0 is tested before any division by u or g-1.
    zero_lam=-s.expand(vals[15].subs({u:0,lam:0}))
    eq(vals[15].subs({u:0,lam:zero_lam}),0)
    zero7=primitive(vals[7].subs({u:0,lam:zero_lam}),h)
    positive=91417638312000*h**4-612612000*h**3+157080*h*h-231*h+2
    eq(zero7.as_expr(),positive)
    qq1=91417638312000*h*h-612612000*h+78540;qq2=78540*h*h-231*h+2
    eq(positive,h*h*qq1+qq2)
    discriminants=[int(s.discriminant(z,h)) for z in [qq1,qq2]]
    assert all(z<0 for z in discriminants)
    assert all(s.Poly(z,h).LC()>0 for z in [qq1,qq2])
    write('interpolation_u0',dict(multiplicities=mult,partial_fractions_exact=True,newton_interpolation_exact=True,
        all_universal_identities=True,u0_positive_quartic=str(positive),u0_discriminants=discriminants,
        negative_control='L4 identity +1 rejected'))
    kval=sum(cc[i]*h**i*qs[i] for i in range(9));kone=-h*qs[1]+h*h*qs[2]-h**4*qs[4]+h**5*qs[5]
    zz=3*h+s.Rational(27,8)*vals[15]/u**3
    dd=(3*vals[7]/u-50*h+16*zz)/7
    raw={
        'F':54*u*vals[11]-222*u*u*vals[7]-4509*vals[15]-1916*h*u**3,
        'E':3*u*u*kval-450*h**8*u**3-s.Rational(2025,4)*h**7*vals[15],
        'Q':u**5*(3*kone-h**6*u*(4*dd*dd+(47*h-4*zz)*dd+zz*zz-36*h*zz-21*h*h)),
        'R':u**6*(18324*vals[6]/u+38988*vals[10]/u**2-195372*vals[14]/u**3-combo.subs({D:dd,Z:zz}))}
    eq(s.diff(raw['F'],lam),-27*(760293059706180000*u*u-343854447300*u+167))
    raw={k:primitive(z,lam,h,u).as_expr() for k,z in raw.items()}
    prod=json.loads((SOURCE/'upper_degree8.json').read_text())
    for k,z in raw.items():eq(z,primitive(s.sympify(prod['raw'][k]).subs(s.Symbol('ell'),lam),lam,h,u).as_expr())
    f=s.Poly(raw['F'],lam);f1,f0=f.all_coeffs()
    reduced={};factors={}
    for name in ['E','Q','R']:
        z=s.Poly(raw[name],lam);degree=z.degree()
        det=sum(z.nth(j)*(-f0)**j*f1**(degree-j) for j in range(degree+1))
        # Clear only powers known to be nonzero on this branch.
        poly=primitive(det,h,u)
        ih=min(m[0] for m,co in poly.terms());iu=min(m[1] for m,co in poly.terms())
        factors[name]=[ih,iu]
        reduced[name]=primitive(poly.as_expr()/(h**ih*u**iu),h,u).as_expr()
        eq(poly.as_expr(),h**ih*u**iu*reduced[name])
    assert factors=={'E':[0,1],'Q':[1,2],'R':[0,2]}
    certificates=json.loads((SOURCE/'upper_degree8_resultant.json').read_text())
    for k,z in reduced.items():eq(z,primitive(s.sympify(certificates['inputs'][k]),h,u).as_expr())
    write('constraints',dict(raw={k:str(z) for k,z in raw.items()},reduced={k:str(z) for k,z in reduced.items()},
      divided_factors=factors,producer_coefficient_comparison=True,lambda_degrees={k:s.degree(z,lam).__int__() for k,z in raw.items()},
      reduced_degrees={k:[int(s.degree(z,h)),int(s.degree(z,u))] for k,z in reduced.items()},
      lambda_pivot_never_divided=True))

if __name__=='__main__':
    START=time.monotonic();constraints()
