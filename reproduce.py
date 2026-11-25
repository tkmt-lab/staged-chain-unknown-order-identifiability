#!/usr/bin/env python3
"""Offline, standard-library reproduction of bounded algebra checks.

Run from any directory. Copied source programs run in private output folders;
source, input snapshots and baselines are never modified. No proof certification.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
NAMES = ('darboux', 'equilibrium', 'boundary', 'weighted', 'recovery')
VOLATILE = {'elapsed_seconds', 'maxrss_bytes', 'peak_rss_bytes', 'ru_maxrss_platform_units'}


def scientific(value):
    if isinstance(value, dict):
        return {k: scientific(v) for k, v in value.items() if k not in VOLATILE}
    if isinstance(value, list):
        return [scientific(v) for v in value]
    return value


def integrity(root):
    count = 0
    for line in (root/'MANIFEST.sha256').read_text().splitlines():
        expected, name = line.split('  ', 1)
        p = root/name
        if p.is_symlink() or not p.is_file():
            raise ValueError('missing or linked manifest entry: '+name)
        if hashlib.sha256(p.read_bytes()).hexdigest() != expected:
            raise ValueError('hash mismatch: '+name)
        count += 1
    return count


def limits():
    resource.setrlimit(resource.RLIMIT_CPU, (60, 60))


def run(output):
    if sys.flags.optimize:
        raise ValueError('assertions must be enabled; do not use -O')
    count = integrity(ROOT)
    output.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    records = []
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONNOUSERSITE='1')
    env.pop('PYTHONPATH', None)
    env.pop('PYTHONOPTIMIZE', None)
    for name in NAMES:
        work = output/name
        work.mkdir(exist_ok=True)
        copied = work/(name+'.py')
        shutil.copyfile(ROOT/'checks'/(name+'.py'), copied)
        r = subprocess.run([sys.executable, '-I', '-B', str(copied)],
                           cwd=work, env=env, capture_output=True, text=True,
                           timeout=60, preexec_fn=limits)
        (work/'stdout.txt').write_text(r.stdout)
        (work/'stderr.txt').write_text(r.stderr)
        if r.returncode:
            raise RuntimeError(name+' failed: '+r.stderr[-2000:])
        if name == 'darboux':
            actual = json.loads((work/'independent-results.json').read_text())
        else:
            actual = json.loads(r.stdout)
        expected = json.loads((ROOT/'baseline'/(name+'.json')).read_text())
        if scientific(actual) != scientific(expected):
            raise AssertionError(name+': scientific output differs from baseline')
        (work/'result.json').write_text(json.dumps(actual, indent=2)+'\n')
        rss = next((actual[k] for k in ('peak_rss_bytes','maxrss_bytes','ru_maxrss_platform_units') if k in actual), None)
        if rss is not None and 'peak_rss_bytes' not in actual and sys.platform != 'darwin':
            rss *= 1024
        records.append({'check':name, 'baseline_match':True,
                        'elapsed_seconds':actual.get('elapsed_seconds'),
                        'peak_rss_bytes':rss})
    result = {'status':'PASS', 'manifest_entries':count, 'checks':records,
              'elapsed_seconds':time.perf_counter()-start,
              'cpu_limit_seconds_per_child':60, 'wall_limit_seconds_per_child':60,
              'rss_limit_enforced':False, 'cost_usd':0,
              'network_required':False, 'formal_proof_certification':False}
    (output/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--extensions', action='store_true')
    args = parser.parse_args()
    dest = args.output.resolve()
    # Avoid accidental replacement of shipped files or protected directories.
    if dest == ROOT or any(dest == ROOT/name or ROOT/name in dest.parents
                           for name in ('checks','baseline','inputs','validation')):
        parser.error('choose a separate output directory, e.g. build/checks')
    try:
        print(json.dumps(run(dest), indent=2))
        if args.extensions:
            sys.stdout.flush(); os.execv(sys.executable, [sys.executable, '-B', str(ROOT/'extensions/reproduce.py'), '--output', str(dest/'extensions')])
    except (ValueError, AssertionError, RuntimeError, subprocess.TimeoutExpired) as exc:
        print('FAIL:', exc, file=sys.stderr)
        raise SystemExit(1)
