#!/usr/bin/env python3
"""Build the bilingual reading layer; never run scientific calculations."""
from pathlib import Path
import json,re,os,subprocess
from urllib.parse import urlsplit
ROOT=Path(__file__).resolve().parents[1];DIST=ROOT/'dist'
mapping=json.loads((ROOT/'scripts/reading-ui.json').read_text());mapping.update(json.loads((ROOT/'scripts/reading-dynamic.json').read_text()))
pattern=re.compile('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)))
def translate(s):
 s=pattern.sub(lambda m:mapping[m[0]],s)
 s=re.sub(r'Таблицы, численные результаты и программы · (\d+) файлов',r'Tables, numerical results and programs · \1 files',s)
 return re.sub(r'(\d+) КБ',r'\1 KiB',s)
def strip(s):return re.sub(r'<!-- reader:start -->.*?<!-- reader:end -->','',s,flags=re.S)
def relative(target,page):return os.path.relpath(DIST/target,(DIST/page).parent).replace(os.sep,'/')
pages=['index.html','materials.html']+[f'viewers/{x}.html' for x in ['hf1','qc0','yitian','b5']]+[f'reports/{p.name}' for p in sorted((DIST/'reports').glob('*.html'))]
def decorate(s,page):
 en=page.startswith('en/');original=page[3:] if en else page;ru=relative(original,page);eng=relative('en/'+original,page)
 bar=f'''<!-- reader:start --><div class="reader-bar" aria-label="{'Reading preferences' if en else 'Настройки чтения'}"><nav aria-label="{'Language' if en else 'Язык'}"><a href="{ru}" data-reader-language="ru" lang="ru" {'aria-current="page"' if not en else ''}>Русский</a><a href="{eng}" data-reader-language="en" lang="en" {'aria-current="page"' if en else ''}>English</a></nav><label for="reader-theme">{'Theme' if en else 'Тема'} <select id="reader-theme"><option value="system">{'System' if en else 'Как в системе'}</option><option value="light">{'Light' if en else 'Светлая'}</option><option value="dark">{'Dark' if en else 'Тёмная'}</option><option value="warm">{'Warm' if en else 'Тёплая'}</option></select></label></div><!-- reader:end -->'''
 head=f'<!-- reader:start --><link rel="alternate" hreflang="ru" href="{ru}"><link rel="alternate" hreflang="en" href="{eng}"><link rel="stylesheet" href="{relative("assets/reader.css",page)}"><script src="{relative("assets/reader.js",page)}"></script><!-- reader:end -->'
 if '</head>' in s:s=s.replace('</head>',head+'</head>',1)
 else:s=s.replace('<main>',head+'<main>',1)
 if re.search(r'<body[^>]*>',s):s=re.sub(r'(<body[^>]*>)',lambda m:m[0]+bar,s,count=1)
 else:s=s.replace('<main>',bar+'<main>',1)
 return s

def rebase(s,old,new):
 def url(m):
  val=m[2];u=urlsplit(val)
  if u.scheme or u.netloc or not u.path or u.path.startswith('/'):return m[0]
  target=os.path.normpath(str(Path(old).parent/u.path)).replace(os.sep,'/')
  if target in pages or target=='viewers/b5-controller.js':target='en/'+target
  return m[1]+relative(target,new)+('?' +u.query if u.query else '')+('#'+u.fragment if u.fragment else '')+m[3]
 return re.sub(r'((?:href|src|poster)=["\'])([^"\']+)(["\'])',url,s)
for page in pages:
 p=DIST/page;s=strip(p.read_text())
 if page=='viewers/hf1.html' and '/* reader-layout */' not in s:
  s=s.replace('function layout(){const axis=',"function layout(){/* reader-layout */const rp=window.HorizonReader.palette();const axis=")
  for a,b in [("color:'#9fb0c2'","color:rp.ink"),("backgroundcolor:'#071421'","backgroundcolor:rp.bg"),("gridcolor:'#24374a'","gridcolor:rp.line"),("paper_bgcolor:'#071421'","paper_bgcolor:rp.bg"),("plot_bgcolor:'#071421'","plot_bgcolor:rp.bg"),("font:{color:'#e8eff6'}","font:{color:rp.ink}")]:s=s.replace(a,b)
 if page=='viewers/qc0.html' and '/* reader-layout */' not in s:
  s=s.replace(' const current=traces();'," /* reader-layout */const rp=window.HorizonReader.palette();layout.paper_bgcolor=rp.bg;layout.font={color:rp.ink};layout.legend.bgcolor=rp.paper;for(const a of ['xaxis','yaxis','zaxis'])Object.assign(layout.scene[a],{color:rp.ink,backgroundcolor:rp.bg,gridcolor:rp.line});\n const current=traces();")
 if page=='viewers/yitian.html' and '/* reader-curve */' not in s:
  s=s.replace('function drawCurve(frame){', 'function drawCurve(frame){/* reader-curve */const rp=window.HorizonReader.palette();')
  s=s.replace("ctx.strokeStyle='#e4ebf0'",'ctx.strokeStyle=rp.line').replace("ctx.fillStyle='#657b8d'",'ctx.fillStyle=rp.ink').replace("ctx.strokeStyle='#263d50'",'ctx.strokeStyle=rp.ink')
  s=s.replace('frameCount:NF};render();',"frameCount:NF};document.addEventListener('reader-theme-change',()=>drawCurve(+$('time').value));render();")
 if not page.startswith('reports/'):
  chunks=re.split(r'(<script\b[^>]*type=["\']application/json["\'][^>]*>.*?</script>)',s,flags=re.S|re.I)
  e=''.join(c if re.match(r'<script\b',c,re.I) else translate(c) for c in chunks)
  e=e.replace(".replace('.',',')",'').replace('lang="ru"','lang="en"').replace("toLocaleString('ru-RU'","toLocaleString('en-GB'")
  parts=re.split(r'(<script\b[^>]*>.*?</script>|<style\b[^>]*>.*?</style>|<[^>]+>)',e,flags=re.S|re.I)
  e=''.join(c if c.startswith('<') else re.sub(r'(?<=\d),(?=\d)', '.', c) for c in parts)
  e=rebase(e,page,'en/'+page);target=DIST/'en'/page;target.parent.mkdir(parents=True,exist_ok=True);target.write_text(decorate(e,'en/'+page))
 p.write_text(decorate(s,page))
p=DIST/'en/viewers/b5-controller.js';p.write_text(translate((DIST/'viewers/b5-controller.js').read_text()))
subprocess.run(['node',str(ROOT/'scripts/render_english.mjs')],check=True)
for p in (DIST/'en/reports').glob('*.html'):p.write_text(decorate(strip(p.read_text()),str(p.relative_to(DIST))))
print('Bilingual reading layer built; no scientific computations run.')
