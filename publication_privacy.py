from pathlib import Path
import json
root=Path(__file__).resolve().parent
for line in (root/'MANIFEST.sha256').read_text().splitlines():
    name=line.split('  ',1)[1]
    assert not any(p in name for p in ('inputs/','instruction','review-inputs','CANDIDATES','INTERNAL_GOAL','goal-tree'))
    p=root/name
    if p.suffix in ('.md','.tex','.json','.py','.yml','.cff') and p.name != 'publication_privacy.py':
        text=p.read_text()
        assert not any(s in text for s in ('/Users/','/Volumes/','mailbox.sqlite','INTERNAL_GOAL','CANDIDATES.md','goal-tree')),(name,'private marker')
print(json.dumps({'status':'PASS','scope':'manifested paths and text; PDF separately inspected'}))
