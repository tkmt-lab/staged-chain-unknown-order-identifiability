"""Distribution privacy: text, JSON keys and values, and decoded gzip JSON."""
from pathlib import Path, PurePosixPath
import gzip, json, re, sys
ROOT = Path(__file__).resolve().parent
def check(root):
    lines = (root/'MANIFEST.sha256').read_text().splitlines()
    failures = []
    decoded = 0
    # Construct forbidden prefixes so this scanner can scan its own source.
    prefixes = ['/' + s + '/' for s in ('workspace', 'home', 'Users', 'Volumes', 'tmp')]
    markers = ['mailbox' + '.sqlite', 'INTERNAL' + '_GOAL', 'CANDIDATES' + '.md', 'goal' + '-tree']
    envelopes = {'limits', 'writes', 'external_contact'}
    for line in lines:
        _, name = line.split('  ', 1)
        path = PurePosixPath(name)
        if path.is_absolute() or '..' in path.parts:
            failures.append([name, 'unsafe manifest path']); continue
        p = root/name
        if p.suffix == '.gz':
            text = gzip.decompress(p.read_bytes()).decode('utf-8')
            json.loads(text); decoded += 1
        else:
            try: text = p.read_text(encoding='utf-8')
            except UnicodeDecodeError:
                if p.suffix != '.pdf': failures.append([name, 'unreviewed binary'])
                continue
        if any(s in text for s in prefixes + markers) or re.search(r'[A-Za-z]:\\(?:Users|workspace)\\', text):
            failures.append([name, 'private absolute path or marker'])
        if p.suffix == '.json':
            obj = json.loads(text)
            if p.name == 'INPUT.json' and isinstance(obj, dict) and envelopes.intersection(obj):
                failures.append([name, 'operational envelope in mathematical input'])
    if failures: raise ValueError(json.dumps(failures))
    return {'status':'PASS', 'manifest_entries':len(lines), 'gzip_json_decoded':decoded,
            'scope':'all manifested UTF-8 text including JSON keys; decoded gzip JSON; PDF requires separate typesetting/text review; no Git-history audit'}
if __name__ == '__main__':
    assert not sys.flags.optimize
    print(json.dumps(check(ROOT)))
