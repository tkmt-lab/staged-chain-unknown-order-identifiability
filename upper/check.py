#!/usr/bin/env python3
"""Hash verification or bounded sequential replay in a fresh copy (two processes)."""
import argparse, hashlib, json, os, resource, shutil, subprocess, sys, time
from pathlib import Path

def limits():
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3))
    resource.setrlimit(resource.RLIMIT_CPU,(590,590))

def main():
    assert not sys.flags.optimize, 'do not use -O'
    root=Path(__file__).resolve().parent
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--group',choices=['u8','u11','u15','u15-without-crt','smoke','all']);ap.add_argument('--job');ap.add_argument('--output');a=ap.parse_args()
    rows=json.loads((root/'SOURCE-MAP.json').read_text())['files']
    for r in rows:
        p=root.parent/r['path'];assert p.is_file(),r['path'];assert p.stat().st_size==r['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],r['path']
    if a.verify:print(json.dumps({'status':'PASS','files':len(rows),'scope':'integrity, not mathematical recomputation'}));return
    jobs=json.loads((root/'JOBS.json').read_text())
    smoke=['u8-general','u8-exception','u8-classify_exact','u11-recompute_reduction','u11-audit_real_roots','u15-degenerate','u15-controls']
    if a.job:jobs=[j for j in jobs if j['id']==a.job]
    elif a.group=='u15-without-crt':jobs=[j for j in jobs if j['id'].startswith('u15-') and not j['id'].startswith('u15-crt')]
    elif a.group=='smoke':jobs=[j for j in jobs if j['id'] in smoke]
    elif a.group and a.group!='all':jobs=[j for j in jobs if j['id'].startswith(a.group+'-')]
    else:assert a.group=='all','choose --verify, --group or --job'
    assert jobs and a.output
    out=Path(a.output).resolve();assert not out.exists(),'output must be new'
    shutil.copytree(root,out/'work');env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    # Keep numeric input access inside the extracted copy, even when another source tree is present.
    launcher="""import os,sys,runpy,pathlib,json
from importlib.metadata import distribution,PackageNotFoundError
allowed_meta=set()
try:allowed_meta.add(pathlib.Path(distribution('numpy')._path).resolve()/'direct_url.json')
except PackageNotFoundError:pass
root=pathlib.Path.cwd().resolve();script=pathlib.Path(sys.argv[1]).resolve();sys.path.insert(0,str(script.parent));sys.argv=sys.argv[1:];seen=set()
def audit(event,args):
 if event=='open' and isinstance(args[0],(str,bytes)):
  p=pathlib.Path(os.fsdecode(args[0])).resolve()
  if p.suffix in ('.json','.gz'):
   if p in allowed_meta:return
   if not p.is_relative_to(root):raise RuntimeError('numeric input outside isolated copy: '+str(p))
   seen.add(str(p.relative_to(root)))
sys.addaudithook(audit)
try:runpy.run_path(str(script),run_name='__main__')
finally:print('INPUT_PATHS '+json.dumps(sorted(seen)))
"""
    results=[]
    for j in jobs:
        start=time.monotonic();print('RUN',j['id'],flush=True)
        with (out/(j['id']+'.log')).open('w') as log:
            try:r=subprocess.run([sys.executable,'-B','-c',launcher,*j['argv']],cwd=out/'work',env=env,stdout=log,stderr=subprocess.STDOUT,timeout=595,preexec_fn=limits);code=r.returncode
            except subprocess.TimeoutExpired:code='timeout'
        results.append(dict(id=j['id'],exit=code,seconds=time.monotonic()-start))
        (out/'RUN.json').write_text(json.dumps({'jobs':results,'status':'PASS' if len(results)==len(jobs) and all(r['exit']==0 for r in results) else 'INCOMPLETE','child_peak_rss_KiB':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss},indent=2)+'\n')
        print(results[-1],flush=True)
        if code:raise SystemExit(1)
    print('PASS',len(jobs),'bounded jobs; prose obligations remain as documented')
if __name__=='__main__':main()
