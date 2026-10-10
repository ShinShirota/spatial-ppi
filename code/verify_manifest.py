"""Verify the package's recorded files without executing scientific code."""
from pathlib import Path
import hashlib,json,sys
root=Path(__file__).resolve().parents[1]
rows=json.loads((root/'manifest.json').read_text())
fail=[]
for row in rows:
    p=root/row['path']
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:
        fail.append(row['path'])
print(json.dumps({'checked':len(rows),'failures':fail},indent=2))
sys.exit(bool(fail))
