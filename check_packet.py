#!/usr/bin/env python3
"""Bounded integrity, source-reference, and deliberate-error controls."""
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import tempfile
import reproduce

ROOT = Path(__file__).resolve().parent


def main():
    checks = []
    n = reproduce.integrity(ROOT)
    checks.append('manifest hashes match')
    provenance = json.loads((ROOT/'PROVENANCE.json').read_text())
    for entry in provenance['inputs']:
        b = (ROOT/entry['copy']).read_bytes()
        assert len(b) == entry['bytes']
        assert hashlib.sha256(b).hexdigest() == entry['sha256']
    checks.append('all pinned input and checker copies match provenance')
    tex = (ROOT/'paper.tex').read_text()
    labels = re.findall(r'\\label\{([^}]+)\}', tex)
    assert len(labels) == len(set(labels))
    refs = re.findall(r'\\(?:eqref|ref)\{([^}]+)\}', tex)
    assert set(refs) <= set(labels)
    cites = re.findall(r'\\cite\{([^}]+)\}', tex)
    bib = re.findall(r'\\bibitem\{([^}]+)\}', tex)
    assert set(cites) <= set(bib)
    assert not re.search(r'R04[234]|mailbox|next[1-9]|commit|check-[1-4]',tex)
    checks.append('TeX labels, citations and separation from internal workflow')
    checks.append('all computational assertions enabled')
    assert not __import__('sys').flags.optimize
    spec = importlib.util.spec_from_file_location('boundary_control', ROOT/'checks/boundary.py')
    mod = importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    u=mod.u;x=mod.x;J=sum(u)
    velocities={'x':-x*J,'u1':x*J-u[0]}
    velocities.update({'u'+str(i+1):u[i-1]-u[i] for i in range(1,5)})
    ratios=[mod.Poly.cast(1)]+[v*u[0]**-1 for v in u[1:]]
    c=mod.invariants(ratios);t=ratios[1]
    actual=mod.derive(c[2],velocities)
    correct=x*(J*u[0]**-1)*(t*c[1]-c[2]+t**2*mod.F(1,2))
    assert not (actual-correct).terms
    wrong=x*(J*u[0]**-1)*(t*c[1]-c[2]-t**2*mod.F(1,2))
    assert (actual-wrong).terms
    assert actual.terms  # freezing the boundary invariant in the full system fails
    checks.append('wrong sign and frozen-invariant controls rejected')
    build=ROOT/'build';build.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='tamper-',dir=build) as d:
        p=Path(d);(p/'input').write_bytes(b'original')
        (p/'MANIFEST.sha256').write_text(hashlib.sha256(b'original').hexdigest()+'  input\n')
        assert reproduce.integrity(p)==1
        (p/'input').write_bytes(b'tampered')
        try: reproduce.integrity(p)
        except ValueError: pass
        else: raise AssertionError('tampering was not rejected')
    checks.append('altered input rejected by manifest')
    print(json.dumps({'status':'PASS','manifest_entries':n,'checks':checks,
                      'independent_proof_review':False},indent=2))

if __name__=='__main__': main()
