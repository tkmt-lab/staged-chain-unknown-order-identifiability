#!/usr/bin/env python3
"""Sequential bounded certificate replays in a disposable copy."""
import argparse,json,os,resource,shutil,subprocess,sys,time
from pathlib import Path

def limits():
    resource.setrlimit(resource.RLIMIT_AS,(7*1024**3,7*1024**3))
    resource.setrlimit(resource.RLIMIT_CPU,(7000,7000))

def main():
    if sys.flags.optimize:raise ValueError('assertions must be enabled')
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);a=ap.parse_args()
    src=Path(__file__).resolve().parent;out=Path(a.output).resolve()
    if out==src or src in out.parents:raise ValueError('output must be outside the extension source tree')
    out.mkdir(parents=True,exist_ok=False)
    for folder in ['l18','u42','b22data','l22','u22','u23']:shutil.copytree(src/folder,out/folder)
    jobs=[('lower18','l18',['independent.py']),
          ('upper42-algebra','u42',['check42.py','certificate42.json','--out','CHECK.json']),
          ('upper42-projection','u42',['check_projection42.py','--out','PROJECTION_CHECK.json']),
          ('lower22','l22',['check.py','--output','RESULT.json']),
          ('upper22-bridge','u22',['check.py','bridge']),
          ('upper22-universal','u22',['check.py','universal']),
          ('upper22-interpolation','u22',['algebra.py']),
          ('upper22-elimination','u22',['eliminate.py']),
          ('upper23','u23',['checker.py'])]
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1');env.pop('PYTHONOPTIMIZE',None)
    rows=[];start=time.monotonic()
    for name,folder,cmd in jobs:
        t=time.monotonic()
        try:
            p=subprocess.run([sys.executable,'-B',*cmd],cwd=out/folder,env=env,capture_output=True,text=True,preexec_fn=limits,timeout=7100)
            (out/(name+'.stdout')).write_text(p.stdout);(out/(name+'.stderr')).write_text(p.stderr)
            row=dict(job=name,argv=cmd,exit=p.returncode,seconds=time.monotonic()-t)
        except subprocess.TimeoutExpired as e:
            row=dict(job=name,argv=cmd,exit='timeout',seconds=time.monotonic()-t)
        rows.append(row);(out/'RUN.json').write_text(json.dumps(rows,indent=2)+'\n')
        print(json.dumps(row),flush=True)
        if row['exit']!=0:raise RuntimeError('certificate job failed: '+name)
    result=dict(status='PASS',jobs=len(rows),seconds=time.monotonic()-start,max_child_rss_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,processes='one runner plus one child',threads_per_child=1,child_AS_GiB=7,child_CPU_seconds=7000,child_wall_seconds=7100,scope='Finite algebra and certified witness checks only; analytic implications require the bundled prose proofs. Earlier large CRT computations are not rerun.')
    (out/'SUMMARY.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
