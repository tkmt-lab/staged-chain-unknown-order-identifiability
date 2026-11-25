import os,sys,json,gzip,time,math,resource
from pathlib import Path
resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
sys.path.insert(0,str(Path(__file__).parent/'.deps'))
from flint import fmpz_mpoly_ctx,fmpq_poly,fmpz_poly,nmod_poly,fmpq
from read_input import read,OUT
BASE='sources/b15-close-2026-11-09/'
def data(p):
 b=read(BASE+p)
 return json.loads(gzip.decompress(b) if p.endswith('.gz') else b)
def poly(p,ctx):
 d=data(p); return ctx.from_dict({tuple(e):int(c) for e,c in d['coefficients']})
def save(name,d):
 d.update(seconds=time.monotonic()-START,rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 (OUT/name).write_text(json.dumps(d,indent=2,default=int)+'\n')
 print(name,d.get('status','done'),d['seconds'],flush=True)
START=time.monotonic()
def determinant():
 C=fmpz_mpoly_ctx.get(('S','X','h','k'));S,X,h,k=C.gens();q=S+h
 def Dt(f):return -S*X*f.derivative(0)+(S+k)*X*f.derivative(1)
 Y=[C.constant(1)];Z=[C.constant(1)]
 for i in range(4):
  Y.append(Dt(Y[-1])+(S+k-X)*Y[-1])
  Z.append(q*(Dt(Z[-1])+(S+k-X)*Z[-1])+(i+1)*S*X*Z[-1])
 B=[q*Z[3],q*q*Z[2],q**3*Z[1],q**4,-q**5*Y[3],-q**5*Y[2],-q**5*Y[1],-q**5,Z[4]]
 # Derive jets before restricting S, not by using the supplied W generator.
 C3=fmpz_mpoly_ctx.get(('h','k','r'));hh,kk,rr=C3.gens()
 def restrict(f):return f.compose(C3.constant(1),rr,hh,kk,ctx=C3)
 M=[]
 for j in range(9):
  M.append([restrict(f) for f in B]); B=[S*f.derivative(0)-(S+k)*f.derivative(1) for f in B]
 # Fraction-free elimination, all divisions verified. Source generation scripts not read.
 prev=C3.constant(1);sign=1
 for a in range(8):
  if not M[a][a]:
   b=next(b for b in range(a+1,9) if M[b][a]);M[a],M[b]=M[b],M[a];sign=-sign
  pivot=M[a][a]
  for i in range(a+1,9):
   for j in range(a+1,9):
    numerator=pivot*M[i][j]-M[i][a]*M[a][j]
    quotient,remainder=divmod(numerator,prev);assert not remainder
    M[i][j]=quotient
  prev=pivot;print('pivot',a,'seconds',time.monotonic()-START,flush=True)
 det=sign*M[8][8]
 T=fmpz_mpoly_ctx.get(('h','k','t'));h,k,t=T.gens()
 det=det.compose(h,k,(1+h)*t,ctx=T)
 f,rem=divmod(det,-12288*h**4*(1+h)**10);assert not rem
 expected=poly('exact/F4.json',T);assert f==expected
 rows=[]
 for j in range(1,5):
  f=-h*(1+h)*f.derivative(0)-k*(1+h)*f.derivative(1)-(t+1+k)*f.derivative(2)
  assert f==poly(f'exact/W{j}.json',T);rows.append(list(f.degrees()))
 save('DETERMINANT.json',dict(status='PASS',F4_terms=len(expected),F4_degrees=list(expected.degrees()),W_degrees=rows,method='time derivative numerator recurrence; x derivative before specialization; exact fraction-free determinant'))

def crt(j):
 from flint import nmod_mpoly_ctx
 from sympy import isprime
 C=fmpz_mpoly_ctx.get(('h','k','t'));h,k,t=C.gens()
 F=poly('exact/F4.json',C);Fs=[F]
 for i in range(3):Fs.append(-h*(1+h)*Fs[-1].derivative(0)-k*(1+h)*Fs[-1].derivative(1)-(t+1+k)*Fs[-1].derivative(2))
 def norm(f):
  arr=[0]*13
  for e,c in f.to_dict().items():arr[e[2]]+=abs(int(c))
  return sum(x*x for x in arr)
 B=norm(Fs[0])**6*norm(Fs[j])**6
 cert=data(f'audit/CRT_R0{j}_CERT.json');primes=cert['primes'];assert len(set(primes))==len(primes) and all(isprime(p) for p in primes)
 rawdata=data(f'audit/CRT_R0{j}.json.gz')
 print('candidate metadata', {k:v for k,v in rawdata.items() if k!='coefficients'},flush=True)
 R=C.from_dict({tuple(e)+(0,):int(c) for e,c in rawdata['coefficients']})
 content=int(cert['reconstructed_content']); raw=content*R*h**20*k**57*(h+1)**176
 actualC=max(abs(int(x)) for x in raw.coeffs());M=math.prod(primes)
 assert B==int(cert['B']); assert actualC==int(cert['candidate_raw_max_coefficient']);assert M>B+actualC
 checks=[]
 for p in primes:
  Cp=nmod_mpoly_ctx.get(('h','k','t'),p)
  ff=[Cp.from_dict({e:int(c)%p for e,c in f.to_dict().items()}) for f in [Fs[0],Fs[j]]]
  assert all(int(f.degrees()[2])==12 for f in ff)
  result=ff[0].resultant(ff[1],2)
  expected=Cp.from_dict({e:int(c)%p for e,c in raw.to_dict().items()})
  assert result==expected,(j,p)
  checks.append(p);print('crt',j,len(checks),len(primes),'elapsed',time.monotonic()-START,flush=True)
  (OUT/f'CRT_0{j}_progress.json').write_text(json.dumps(checks))
 save(f'CRT_0{j}.json',dict(status='PASS',primes=primes,B_bits=B.bit_length(),C_bits=actualC.bit_length(),M_bits=M.bit_length(),criterion=True,all_raw_coefficients_recomputed=True))
if __name__=='__main__':
 if len(sys.argv)==1:determinant()
 elif sys.argv[1]=='crt':crt(int(sys.argv[2]))
