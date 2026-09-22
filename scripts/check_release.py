#!/usr/bin/env python3
"""Check packaged routes, scripts and B5 lineage; never rerun research."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import unquote, urlsplit
import hashlib, json, re, subprocess, tempfile
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'dist'
class Links(HTMLParser):
    def __init__(self):super().__init__();self.links=[];self.ids=set();self.code=[];self.inscript=False;self.active=''
    def handle_starttag(self,tag,attrs):
        a=dict(attrs)
        if 'id' in a:self.ids.add(a['id'])
        for k in ('href','src','poster'):
            if k in a:self.links.append(a[k])
        if tag=='script':self.inscript=('src' not in a and a.get('type')!='application/json');self.active=''
    def handle_data(self,data):
        if self.inscript:self.active+=data
    def handle_endtag(self,tag):
        if tag=='script' and self.inscript:self.code.append(self.active);self.inscript=False
errors=[];parsed={};links=0;checked_js=0
for p in D.rglob('*.html'):
    obj=Links();obj.feed(p.read_text());parsed[p.resolve()]=obj
for p,o in parsed.items():
    for url in o.links:
        q=urlsplit(url)
        if q.scheme or q.netloc or url.startswith('data:'):continue
        target=(p.parent/unquote(q.path)).resolve() if q.path else p
        links+=1
        if not target.exists():errors.append(f'{p.relative_to(D)} -> missing {url}')
        elif q.fragment and target in parsed and q.fragment not in parsed[target].ids:errors.append(f'{p.relative_to(D)} -> missing anchor {url}')
    for code in o.code:
        if not code.strip():continue
        run=subprocess.run(['node','--check'],input=code,text=True,capture_output=True)
        checked_js+=1
        if run.returncode:errors.append(f'{p.name} inline JS: {run.stderr[:400]}')
for p in D.rglob('*.js'):
    run=subprocess.run(['node','--check',str(p)],capture_output=True,text=True);checked_js+=1
    if run.returncode:errors.append(f'{p.relative_to(D)} JS: {run.stderr[:300]}')
b5=ROOT/'research/hf4_b5/input/BH_HF4_B4_fields_geometry_masks.npz'
digest=hashlib.sha256(b5.read_bytes()).hexdigest()
assert digest=='a5e1c19802ded7f79da5475eb267f1228c6e26208f69c88337346faf617ee8e8'
sizes=sorted([(p.stat().st_size,str(p.relative_to(D))) for p in D.rglob('*') if p.is_file()],reverse=True)
result={'html_pages':len(parsed),'local_references_checked':links,'javascript_syntax_checks':checked_js,
        'errors':errors,'largest_static_asset':sizes[0],'static_bytes':sum(s for s,_ in sizes),
        'B5_input_sha256':digest,'previous_scientific_analysis_recomputed':False,
        'browser_WebGL_render_verified_this_release':False,
        'browser_limitation':'Static preview is not supported by the available Sites preview; native publishing services return an access error.'}
(ROOT/'docs/release-validation.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result,ensure_ascii=False,indent=2));raise SystemExit(bool(errors))
