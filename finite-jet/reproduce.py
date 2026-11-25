#!/usr/bin/env python3
"""Replay bounded certificates in a disposable output copy, never in source."""
import argparse,json,os,resource,shutil,subprocess,sys,time
from pathlib import Path

def limit():
    resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
    resource.setrlimit(resource.RLIMIT_CPU,(590,590))

def main():
    assert not sys.flags.optimize,'do not use -O'
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);a=ap.parse_args()
    out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False)
    src=Path(__file__).resolve().parent
    for name in ['one','two','three','four']:shutil.copytree(src/name,out/name)
    jobs=[('one',['verify.py']),('two',['verify.py','pair.json','--out','RESULT.json']),
          ('three',['certify_interval.py','lower_candidate.json','lower_root_certificate.json']),
          ('three',['realize.py']),('four',['verify.py'])]
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
    rows=[]
    for i,(folder,args) in enumerate(jobs):
        start=time.monotonic();cmd=[sys.executable,'-B',*args]
        p=subprocess.run(cmd,cwd=out/folder,env=env,capture_output=True,text=True,timeout=595,preexec_fn=limit)
        (out/f'{i}.stdout').write_text(p.stdout);(out/f'{i}.stderr').write_text(p.stderr)
        rows.append(dict(cwd=folder,argv=['python3','-B',*args],exit=p.returncode,seconds=time.monotonic()-start))
        (out/'RUN.json').write_text(json.dumps(rows,indent=2)+'\n')
        if p.returncode:raise RuntimeError(f'{folder}: see {out}/{i}.stderr')
    print(json.dumps({'status':'PASS','checks':5,'scope':'upper 4; lower 4,8,11,15; upper 8,11,15 and analytic bridges not re-proved by this runner','output':str(out)}))
if __name__=='__main__':main()
