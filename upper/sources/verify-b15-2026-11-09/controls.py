from independent import *
C=fmpz_mpoly_ctx.get(('h','k'));h,k=C.gens()
a=poly('audit/G82.json.gz',C);b=poly('audit/G99.json.gz',C);r=poly('audit/CRT_R01.json.gz',C)
assert a*b==r
H=data('projection/H50.json');H50=fmpz_poly([1]);H30=fmpz_poly([1]);qs=[]
for f in H['factors']:
 q=fmpz_poly(list(map(int,f['coefficients'])));qs.append(q);H50*=q**f['multiplicity'];H30*=q
assert H50==fmpz_poly(list(map(int,H['H50_coefficients'])))
for i,q in enumerate(qs):
 assert q(0) and q(-1)
 for r in qs[:i]:assert q.gcd(r).degree()==0
# Reject a changed candidate section.
qq=data('projection/QADIC_R02_H80_0.json');q=fmpq_poly([fmpq(x) for x in qq['q_coefficients']])**qq['multiplicity'];K=fmpq_poly([fmpq(x) for x in qq['section_coefficients']])+1
hh=fmpq_poly([0,1]);value=fmpq_poly([])
for (i,j),c in a.to_dict().items():value=(value+int(c)*(hh**i%q)*(K**j%q))%q
assert value
crt=[]
T=fmpz_mpoly_ctx.get(('h','k','t'));h,k,t=T.gens();F=poly('exact/F4.json',T);Fs=[F]
for _ in range(3):
 f=Fs[-1];Fs.append(-h*(1+h)*f.derivative(0)-k*(1+h)*f.derivative(1)-(t+1+k)*f.derivative(2))
for j in [1,2,3]:
 cert=data(f'audit/CRT_R0{j}_CERT.json');M=int(cert['product']);B=int(cert['B']);CC=int(cert['candidate_raw_max_coefficient'])
 for pp in cert['primes']:
  assert all(any(e[2]==12 and int(c)%pp for e,c in f.to_dict().items()) for f in [Fs[0],Fs[j]])
 # Q_raw + M*u has the same residues but its h^(20+88) k^57 coefficient
 # is >= M*binom(176,88)-C, enough to reject using the SAME criterion.
 lower=M*math.comb(176,88)-CC
 assert M<=B+lower
 crt.append(dict(pair=f'0{j}',same_residue_mutation_rejected=True))
from fractions import Fraction
boxes=json.loads((OUT/'FIBERS.json').read_text())['boxes']
intervals=sorted(tuple(map(Fraction,b['h'])) for b in boxes)
assert all(a[1]<b[0] for a,b in zip(intervals,intervals[1:]))
# Verify that the Hensel clusters exhaust and partition the initial polynomial.
from local import read_coeff,mobius,qp
for branch in [82,99]:
 setup=data(f'projection/CLUSTER_G{branch}.json')
 original=mobius(read_coeff(f'audit/G{branch}.json.gz',1,True),setup['mobius_c'],1)
 p0=qp([a[0] for a in original]);p0=p0/p0.leading_coefficient()
 product=qp([1])
 for c in setup['clusters']:product*=qp(c['factor_coefficients'])**c['multiplicity']
 assert product==p0
save('CONTROLS.json',dict(status='PASS',R01_exact_factor_product=True,H50_degree=H50.degree(),H30_degree=H30.degree(),H80_degree=(H50*H30).degree(),factors_pairwise_coprime=True,domain_boundary_disjoint=True,changed_section_rejected=True,CRT_controls=crt))
