#!/usr/bin/env python3
"""Validate the reading layer without running scientific calculations."""
from pathlib import Path
from html.parser import HTMLParser
import hashlib,json,re,subprocess
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'dist'
def tracked(*paths):return subprocess.check_output(['git','ls-files','-z',*paths],cwd=ROOT).decode().rstrip('\0').split('\0')
def baseline(name):return subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)
files=[p for p in D.rglob('*') if p.is_file()];before={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
subprocess.run(['python',str(ROOT/'scripts/localize_site.py')],check=True,capture_output=True)
assert all(hashlib.sha256(p.read_bytes()).hexdigest()==before[str(p)] for p in files),'Non-idempotent build'
class Visible(HTMLParser):
 def __init__(self):super().__init__();self.skip=0;self.text=[]
 def handle_starttag(self,t,a):
  if t in ['style','script']:self.skip+=1
 def handle_endtag(self,t):
  if t in ['style','script']:self.skip-=1
 def handle_data(self,s):
  if not self.skip and re.search('[А-Яа-яЁё]',s) and s.strip()!='Русский' and not s.strip().startswith('prior/Пространства/'):self.text.append(s.strip())
left={}
for p in (D/'en').rglob('*.html'):
 q=Visible();q.feed(p.read_text())
 if q.text:left[str(p.relative_to(D))]=q.text
assert not left,left
protected=tracked('research','dist/downloads')
for name in protected:assert (ROOT/name).read_bytes()==baseline(name),name
pat=r'<script\b[^>]*type=["\']application/json["\'][^>]*>(.*?)</script>'
for name in ['hf1','qc0','yitian','b5']:
 path='dist/viewers/'+name+'.html';base=baseline(path).decode()
 for p in [ROOT/path,D/'en/viewers'/f'{name}.html']:assert re.findall(pat,base,re.S)==re.findall(pat,p.read_text(),re.S),str(p)
media=[]
for name in tracked('dist/assets','dist/viewers'):
 if name.endswith(('.html','.css','reader.js','b5-controller.js')):continue
 assert (ROOT/name).read_bytes()==baseline(name),name
 media.append(name)
def lum(h):
 rgb=[int(h[i:i+2],16)/255 for i in (1,3,5)]
 return sum(v*w for v,w in zip([v/12.92 if v<=.04045 else ((v+.055)/1.055)**2.4 for v in rgb],[.2126,.7152,.0722]))
colors={'light':['#f8fafc','#ffffff','#eaf0f5','#203247','#4e6378','#195c92'],'dark':['#151b24','#1d2631','#253141','#dee6ef','#b4c2d1','#95c7ee'],'warm':['#ede5d7','#f5eddf','#e2d7c5','#39352e','#655c4f','#27576c']};contrast={}
for theme,vs in colors.items():
 ratios=[(max(lum(fg),lum(bg))+.05)/(min(lum(fg),lum(bg))+.05) for bg in vs[:3] for fg in vs[3:]];contrast[theme]=round(min(ratios),3);assert min(ratios)>=4.5
p=ROOT/'docs/reading-validation.json';old=json.loads(p.read_text()) if p.exists() else {}
r={'date':'2026-09-23','html_pages':len(list(D.rglob('*.html'))),'english_reports':18,'builder_idempotent':True,'untranslated_visible_prose':left,'scientific_files_byte_identical':len(protected),'standalone_data_and_media_byte_identical':len(media),'embedded_scientific_JSON_byte_identical':True,'minimum_normal_text_contrast':contrast,'scientific_computations_repeated':False,'browser_QA':old.get('browser_QA','Pending deployed-site check')};p.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n');print(json.dumps(r,indent=2))
