from independent import *
from flint import fmpz_mat
C=fmpz_mpoly_ctx.get(('S','X','k'));S,X,k=C.gens()
D=lambda f:-S*X*f.derivative(0)+(S+k)*X*f.derivative(1)
# At h=0, z=X. Time derivatives X and y, divided by X before multiplying S^4.
Z=[C.constant(1)];Y=[S]
for j in range(4):
 Z.append(D(Z[-1])+(S+k)*Z[-1]);Y.append(D(Y[-1])+(S+k)*Y[-1])
B=[S**4*z for z in Z]+[S**4*y for y in Y[:4]]
T=fmpz_mpoly_ctx.get(('S','L','r','k'));s,l,r,kk=T.gens()
support=set()
for b in B:
 for e in b.compose(s,r+1-s-kk*l,kk,ctx=T).monoms():support.add((e[0],e[1]))
basis=[(4,0)]+[(i,j) for i in range(5,9) for j in range(9-i)]
assert support<=set(basis)
A=fmpz_mat([[math.factorial(n)//math.factorial(n-j)*i**(n-j) if n>=j else 0 for i,j in basis] for n in range(11)])
assert A.det()==1146617856
T3=fmpz_mpoly_ctx.get(('h','k','t'));h,k,t=T3.gens();F=poly('exact/F4.json',T3)
lc=T3.from_dict({(e[0],e[1],0):int(c) for e,c in F.to_dict().items() if e[2]==12})
L,rem=divmod(lc,-3359232000*h**6*k*(1+h)**2);assert not rem
lk=L.subs({0:2});assert lk==6*(19399*k*k-2860*k-414120)
assert lk(0,-5,0)>0 and lk(0,-4,0)<0
save('DEGENERATE.json',dict(status='PASS',h0_basis=basis,h0_support=sorted(support),h0_jet_determinant=int(A.det()),lc_has_admissible_zero=True))
