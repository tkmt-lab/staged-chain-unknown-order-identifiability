#!/usr/bin/env python3
"""Independent contact verifier: Python integer interval arithmetic, no research imports.
All endpoints are integers divided by 2**BITS. Round every primitive outwards.
"""
import json, math, hashlib, time, resource, signal, sys
from pathlib import Path
from fractions import Fraction as F
BITS=1024
UNIT=1<<BITS
BASE=Path(__file__).resolve().parent
ROOT=BASE
SOURCE=BASE

def ceildiv(a,b): return -((-a)//b)
class I:
    __slots__=('l','u')
    def __init__(self,a=0,b=None,raw=False):
        if raw: self.l,self.u=a,b; return
        if isinstance(a,I): self.l,self.u=a.l,a.u; return
        a=F(a); b=a if b is None else F(b)
        assert a<=b
        self.l=(a.numerator*UNIT)//a.denominator
        self.u=ceildiv(b.numerator*UNIT,b.denominator)
    def __add__(self,b):
        b=I(b); return I(self.l+b.l,self.u+b.u,True)
    __radd__=__add__
    def __neg__(self): return I(-self.u,-self.l,True)
    def __sub__(self,b): return self+-I(b)
    def __rsub__(self,b): return I(b)+-self
    def __mul__(self,b):
        b=I(b); p=[self.l*b.l,self.l*b.u,self.u*b.l,self.u*b.u]
        return I(min(p)//UNIT,ceildiv(max(p),UNIT),True)
    __rmul__=__mul__
    def inv(self):
        assert self.l>0 or self.u<0, 'division by zero interval'
        return I((UNIT*UNIT)//self.u,ceildiv(UNIT*UNIT,self.l),True)
    def __truediv__(self,b): return self*I(b).inv()
    def __rtruediv__(self,b): return I(b)*self.inv()
    def zero(self): return self.l<=0<=self.u
    def mid(self): return F(self.l+self.u,2*UNIT)
    def bound(self): return F(max(abs(self.l),abs(self.u)),UNIT)
    def out(self): return [str(F(self.l,UNIT)),str(F(self.u,UNIT))]
    def display(self): return [float(F(self.l,UNIT)),float(F(self.u,UNIT))]

class AD:
    __slots__=('v','d')
    def __init__(self,v,d=None):
        if isinstance(v,AD): self.v,self.d=v.v,v.d; return
        self.v=I(v); self.d=[I(0) for _ in range(3)] if d is None else d
    def __add__(self,b):
        b=AD(b); return AD(self.v+b.v,[x+y for x,y in zip(self.d,b.d)])
    __radd__=__add__
    def __neg__(self): return AD(-self.v,[-x for x in self.d])
    def __sub__(self,b): return self+-AD(b)
    def __rsub__(self,b): return AD(b)+-self
    def __mul__(self,b):
        b=AD(b);return AD(self.v*b.v,[x*b.v+self.v*y for x,y in zip(self.d,b.d)])
    __rmul__=__mul__
    def __truediv__(self,b):
        b=AD(b); v=self.v/b.v
        return AD(v,[(x-v*y)/b.v for x,y in zip(self.d,b.d)])
    def __rtruediv__(self,b): return AD(b)/self

def value(x): return x.v if isinstance(x,AD) else x

def solve(A,b):
    """Gaussian elimination without pivot swaps; every interval pivot excludes zero."""
    A=[list(row)+[v] for row,v in zip(A,b)]; piv=[]; n=len(b)
    for j in range(n):
        p=A[j][j]; assert not value(p).zero();piv.append(value(p))
        for k in range(j+1,n):
            m=A[k][j]/p
            for l in range(j+1,n+1): A[k][l]=A[k][l]-m*A[j][l]
            A[k][j]=p*0
    x=[None]*n
    for j in reversed(range(n)):
        x[j]=(A[j][n]-sum((A[j][k]*x[k] for k in range(j+1,n)),0))/A[j][j]
    return x,piv

def contact(x):
    h,k,r=x; S=[r*0+1]; X=[r]; Y=[]
    for n in range(16):
        y=sum((math.comb(n,j)*S[j]*X[n-j] for j in range(n+1)),0)
        Y.append(y); S.append(-y); X.append(k*X[n]+y)
    q=S[:];q[0]=q[0]+h; Z=[]
    for n in range(16):
        Z.append((Y[n]-sum((math.comb(n,j)*q[j]*Z[n-j] for j in range(1,n+1)),0))/q[0])
    rows=[[Z[n+3],Z[n+2],Z[n+1],Z[n],-Y[n+3],-Y[n+2],-Y[n+1],-Y[n]] for n in range(12)]
    u,piv=solve(rows[:8],[-Z[n+4] for n in range(8)])
    res=[Z[n+4]+sum((rows[n][j]*u[j] for j in range(8)),0) for n in range(12)]
    return res,u,Y,Z,piv

def inverse_fraction(A):
    n=len(A); a=[list(row)+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for j in range(n):
        p=a[j][j]; assert p
        a[j]=[v/p for v in a[j]]
        for i in range(n):
            if i!=j:
                m=a[i][j];a[i]=[v-m*w for v,w in zip(a[i],a[j])]
    inv=[row[n:] for row in a]
    assert all(sum(A[i][k]*inv[k][j] for k in range(n))==F(i==j) for i in range(n) for j in range(n))
    return inv

def root_certificate(box,centers):
    ad=[AD(v,[I(i==j) for j in range(3)]) for i,v in enumerate(box)]
    R,u,*_=contact(ad); J=[a.d for a in R[8:11]]
    C=inverse_fraction([[a.mid() for a in row] for row in J])
    center=[I(v) for v in centers]; f=contact(center)[0][8:11]
    E=[[I(i==j)-sum((I(C[i][k])*J[k][j] for k in range(3)),I(0)) for j in range(3)] for i in range(3)]
    L=max(sum(a.bound() for a in row) for row in E)
    image=[center[i]-sum((I(C[i][j])*f[j] for j in range(3)),I(0))+sum((E[i][j]*(box[j]-center[j]) for j in range(3)),I(0)) for i in range(3)]
    assert L<1
    assert all(a.l<b.l and b.u<a.u for a,b in zip(box,image))
    return image,dict(box=box,image=image,preconditioner=[[str(v) for v in row] for row in C],jacobian=J,error=E,contraction=str(L),contraction_display=float(L),pivots=[value(p) for p in contact(ad)[4]])

def mv(A,x): return [sum((v*w for v,w in zip(row,x)),I(0)) for row in A]
def rowmul(v,A): return [sum((v[i]*A[i][j] for i in range(len(v))),I(0)) for j in range(len(v))]
def poly(d,x):
    out=I(1)
    for a in d:out=out*x+a
    return out

def reconstruct(x,pole_boxes):
    R,u,Y,z,piv=contact(x); d=u[:4];g=u[4];assert g.l>0
    poles=[]; pole_proofs=[]
    for ends in pole_boxes:
        a,b=map(F,ends); a0,b0=a,b
        pa,pb=poly(d,I(a)),poly(d,I(b))
        assert (pa.l>0 and pb.u<0) or (pa.u<0 and pb.l>0)
        sign=1 if pa.l>0 else -1
        der=poly([3*d[0]/4,d[1]/2,d[2]/4],I(a,b))*4
        assert not der.zero()
        # Guaranteed interval bisection. Stop when coefficient uncertainty dominates.
        for _ in range(300):
            m=(a+b)/2; pm=poly(d,I(m))
            if pm.zero():break
            if (pm.l>0)==(sign>0):a=m
            else:b=m
        poles.append(I(a,b));pole_proofs.append(dict(input=ends,endpoint_values=[pa,pb],derivative=der,refined=I(a,b)))
    assert all(poles[i].u<poles[i+1].l for i in range(3))
    assert all(p.u<0 for p in poles)
    K=[[poles[i] if i==j else I(i==j+1) for j in range(4)] for i in range(4)]
    sums=[sum((K[i][j] for i in range(4)),I(0)) for j in range(4)]
    assert all(v.u<0 for v in sums)
    # Numerator P(s)=s^3+(b2/g)s^2+(b1/g)s+b0/g; expand P(s)/D(s) at infinity.
    n=[I(1)]+[v/g for v in u[5:]]; H=[]
    for j in range(4):H.append(n[j]-sum((d[i-1]*H[j-i] for i in range(1,j+1)),I(0)))
    e=[I(1),I(0),I(0),I(0)]; cols=[e]
    for j in range(3):cols.append(mv(K,cols[-1]))
    Q=[[cols[j][i] for j in range(4)] for i in range(4)]
    assert all(Q[i][j].l==Q[i][j].u==(UNIT if i==j else 0) for i in range(4) for j in range(i+1))
    C,cp=solve(cols,H);assert C[0].l==C[0].u==UNIT
    O=[C]
    for j in range(3):O.append(rowmul(O[-1],K))
    M=[z[j]/g-sum((H[j-1-l]*Y[l] for l in range(j)),I(0)) for j in range(4)]
    state,op=solve(O,M)
    assert all(v.l>0 for v in C+state)
    h,k,r=x
    recoveries=[I(100)-1-r,I(100)-1-h-sum(state,I(0))]
    assert all(v.l>0 for v in recoveries)
    assert h.l>-UNIT and k.u<0 and r.l>0
    return dict(u=u,residuals=R,Y=Y,z=z,linear_pivots=piv,poles=pole_proofs,K=K,C=C,X=state,g=g,Q=Q,H=H,O=O,observability_pivots=op,column_sums=sums,recoveries=recoveries)

def ode(beta,K,C,S0,X0,N=100,order=15):
    S=[I(S0)];X=[[I(v) for v in X0]];Rec=[I(N)-S[0]-sum(X[0],I(0))];Y=[]
    for n in range(order+1):
        y=I(beta)*sum((math.comb(n,j)*S[j]*sum((C[i]*X[n-j][i] for i in range(len(C))),I(0)) for j in range(n+1)),I(0));Y.append(y)
        v=mv(K,X[n]);S.append(-y);X.append([v[i]+(y if i==0 else 0) for i in range(len(C))]);Rec.append(-sum(v,I(0)))
    return dict(Y=Y,S=S,X=X,R=Rec)

def arithmetic_tests():
    count=0
    samples=[F(-7,3),F(-1,10),F(0),F(2,7),F(11,3)]
    intervals=[(a,b) for a in samples for b in samples if a<=b]
    def contains(c,v):return F(c.l,UNIT)<=v<=F(c.u,UNIT)
    for a,b in intervals:
      for c,d in intervals:
        x,y=I(a,b),I(c,d)
        for z,values in [(x+y,[a+c,b+d]),(x-y,[a-d,b-c]),(x*y,[a*c,a*d,b*c,b*d])]:
            assert all(contains(z,v) for v in values);count+=1
        if c*d>0:
            z=x/y;assert all(contains(z,v) for v in [a/c,a/d,b/c,b/d]);count+=1
    try:I(-1,1).inv()
    except AssertionError:count+=1
    else:raise AssertionError('zero divisor accepted')
    a,b=F(2,7),F(-7,3);x=AD(a,[I(1),I(0),I(0)]);y=AD(b,[I(0),I(1),I(0)])
    z=(x*x+y)/(x-y)
    refs=[((2*a)*(a-b)-(a*a+b))/(a-b)**2,((a-b)+(a*a+b))/(a-b)**2,F(0)]
    assert all(contains(v,r) for v,r in zip(z.d,refs));count+=3
    return count

def exact_scaling_test():
    # Auxiliary rational models test the scaling code, not existence of contact.
    def series(beta,K,C,s,x,N):
        S=[s];X=[x];R=[N-s-sum(x)];Y=[]
        for n in range(16):
            y=beta*sum(math.comb(n,j)*S[j]*sum(C[i]*X[n-j][i] for i in range(len(C))) for j in range(n+1))
            Y.append(y);v=[sum(K[i][j]*X[n][j] for j in range(len(C))) for i in range(len(C))]
            S.append(-y);X.append([v[i]+(y if i==0 else 0) for i in range(len(C))]);R.append(-sum(v))
        return S,X,R,Y
    a=F(7,300); count=0
    for dim in [1,4]:
        K=[[F(-i-2) if i==j else F(i==j+1) for j in range(dim)] for i in range(dim)]
        C=[F(i+1) for i in range(dim)];x=[F(1,i+3) for i in range(dim)];s=F(1);beta=F(2,3)
        original=series(beta,K,C,s,x,F(100));scaled=series(beta/a,K,C,a*s,[a*v for v in x],F(7,3))
        for group in [0,2,3]:
            for v,w in zip(original[group],scaled[group]):assert a*v==w;count+=1
        for row,row2 in zip(original[1],scaled[1]):
            for v,w in zip(row,row2):assert a*v==w;count+=1
        for model,pop in [(original,F(100)),(scaled,F(7,3))]:
            for n in range(17):assert model[0][n]+sum(model[1][n])+model[2][n]==(pop if n==0 else 0);count+=1
        bad=series(beta,K,C,a*s,[a*v for v in x],F(7,3))
        assert bad[3][0]!=scaled[3][0]
        assert scaled[0][0]+sum(scaled[1][0])+original[2][0]!=F(7,3)
    return count

def serialize(o):
    if isinstance(o,I):return o.out()
    raise TypeError(type(o))

def main():
    signal.alarm(590);resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
    start=time.monotonic(); tests=arithmetic_tests(); scaling_tests=exact_scaling_test();print('arithmetic checks',tests,'exact scaling checks',scaling_tests,flush=True)
    manifest=json.loads((BASE/'READS.json').read_text());hashes={r['path']:r['sha256'] for r in manifest}
    def read(name):
        p=SOURCE/name;data=p.read_bytes();assert hashlib.sha256(data).hexdigest()==hashes[str(p.relative_to(ROOT))];return json.loads(data)
    src=read('CONTACT11_BOX.json'); chain=read('POSITIVE_CHAIN11.json')
    centers=list(map(F,src['center']));rad=F(src['radius'])
    original=[list(map(F,b)) for b in src['certificate']['box']]
    box=[I(*b) for b in original]
    x,cert=root_certificate(box,centers)
    # Existence inside both the original exact rational box and nominal center +/- radius.
    assert all(F(xi.l,UNIT)>max(a,c-rad) and F(xi.u,UNIT)<min(b,c+rad) for xi,(a,b),c in zip(x,original,centers))
    print('root contraction',cert['contraction_display'],flush=True)
    model=reconstruct(x,[p['box'] for p in chain['pole_isolation']]);h,k,r=x
    L=ode(1,[[k]],[I(1)],1,[r]);R=ode(model['g'],model['K'],model['C'],1+h,model['X'])
    diff=[a-b for a,b in zip(L['Y'],R['Y'])]
    assert all(v.zero() for v in diff[:15]);assert diff[15].l>0
    bridge=(1+h)*model['residuals'][11];assert bridge.l>0
    assert max(bridge.l,diff[15].l)<=min(bridge.u,diff[15].u)
    # A coarse simple rational enclosure certified directly, no float decisions.
    assert F(diff[15].l,UNIT)>128995100000000000 and F(diff[15].u,UNIT)<128995200000000000
    for states in [L,R]:
        assert (states['S'][0]+sum(states['X'][0],I(0))+states['R'][0]-100).zero()
        assert all((states['S'][j]+sum(states['X'][j],I(0))+states['R'][j]).zero() for j in range(1,17))
    print('15th difference',diff[15].display(),'bridge',bridge.display(),flush=True)
    # Deliberately move h; proving residual nonzero is stronger than fixed-point-test failure.
    shifted=[I(c-rad,c+rad) for c in centers];shifted[0]=shifted[0]+I(F(1,10**25))
    bad=contact(shifted)[0][8:11];assert any(not v.zero() for v in bad)
    badstate=model['X'][:];badstate[0]=-badstate[0]
    assert badstate[0].u<0
    badY=model['g']*(1+h)*sum((c*z for c,z in zip(model['C'],badstate)),I(0))
    assert not (L['Y'][0]-badY).zero()
    # Change a transfer numerator coefficient: its first affected Markov coefficient differs.
    alteredH1=(model['u'][5]/model['g']+1)-model['u'][0]
    assert not (alteredH1-model['H'][1]).zero()
    controls=dict(shifted_h_residuals=bad,flipped_state_output_difference=L['Y'][0]-badY,altered_markov_difference=alteredH1-model['H'][1])
    # A second population, computed from the actual reconstructed model.
    a=F(7,300)
    L2=ode(I(1)/a,[[k]],[I(1)],a,[a*r],F(7,3))
    R2=ode(model['g']/a,model['K'],model['C'],a*(1+h),[a*v for v in model['X']],F(7,3))
    diff2=[v-w for v,w in zip(L2['Y'],R2['Y'])]
    assert all(v.zero() for v in diff2[:15]) and diff2[15].l>0
    for j in range(16):
        assert (L2['Y'][j]-a*L['Y'][j]).zero() and (R2['Y'][j]-a*R['Y'][j]).zero()
    result=dict(schema='finite-jet.contact11.v1',bits=BITS,input_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [SOURCE/'CONTACT11_BOX.json',SOURCE/'POSITIVE_CHAIN11.json']},code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),arithmetic_checks=tests,exact_scaling_checks=scaling_tests,root=cert,root_enclosure=x,model=model,left=L,right=R,differences=diff,bridge=bridge,negative_controls=controls,scaled_population='7/3',scaled_differences=diff2,seconds=time.monotonic()-start,rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (BASE/'RESULT.json').write_text(json.dumps(result,default=serialize,indent=2)+'\n')
    print('PASS root, model, jet separation, controls, scaling',result['seconds'],result['rss_KiB'],flush=True)
if __name__=='__main__':main()
