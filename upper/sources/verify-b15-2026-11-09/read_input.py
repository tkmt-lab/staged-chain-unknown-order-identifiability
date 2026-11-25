import sys,json,hashlib,datetime
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/'sources/verify-b15-2026-11-09'
def read(path,display=False):
 p=ROOT/path; b=p.read_bytes()
 row=dict(path=path,sha256=hashlib.sha256(b).hexdigest(),utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),mode='display' if display else 'computation')
 with (OUT/'READS.jsonl').open('a') as f:f.write(json.dumps(row,ensure_ascii=False)+'\n')
 if display:print('\nFILE '+path+'\n'+b.decode())
 return b
if __name__=='__main__':
 for p in sys.argv[1:]:read(p,True)
