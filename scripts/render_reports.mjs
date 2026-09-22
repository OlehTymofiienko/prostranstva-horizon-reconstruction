import fs from 'node:fs';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {marked} from './vendor/marked.mjs';
const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const dist=path.join(root,'dist'),reports=path.join(dist,'reports');
const escape=s=>s.replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('>','&gt;').replaceAll('"','&quot;');
let catalog=JSON.parse(fs.readFileSync(path.join(root,'docs/report-catalog.json'),'utf8'));
for(const name of ['PROJECT_REPORT','REPRODUCIBILITY','SOURCES']){
  fs.copyFileSync(path.join(root,'docs',name+'.md'),path.join(reports,name+'.md'));
  catalog.push({stage:'Итоговый выпуск',title:name,markdown:name+'.md',html:name+'.html'});
}
const reportNames=new Set(fs.readdirSync(reports));
for(const item of catalog){
 let body=marked.parse(fs.readFileSync(path.join(reports,item.markdown),'utf8'));
 // Preserve source text, but remove expired local links from historical reports.
 body=body.replace(/<a href="([^"]+)"([^>]*)>([\s\S]*?)<\/a>/g,(all,url,extra,label)=>{
   if(/^(https?:|mailto:|#)/.test(url))return all;
   const base=path.basename(url.replace(/^sandbox:/,''));
   if(reportNames.has(base))return `<a href="${escape(base)}">${label}</a>`;
   return `<span>${label}</span>`;
 });
 body=body.replace(/<img src="([^"]+)"([^>]*)>/g,(all,url,extra)=>{
   const base=path.basename(url.replace(/^sandbox:/,''));
   return reportNames.has(base)?`<img src="${escape(base)}"${extra}>`:'';
 });
 body=body.replaceAll('<table>','<div class="table-wrap"><table>').replaceAll('</table>','</table></div>');
 const note=item.stage==='Итоговый выпуск'?'':`<p class="status-note">Отчёт сохранён как результат своего этапа. Исправления и границы всего выпуска сведены в <a href="PROJECT_REPORT.html">итоговом отчёте</a>.</p>`;
 const page=`<!doctype html><html lang="ru"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>${escape(item.title)} — Пространства</title><link rel="stylesheet" href="../assets/style.css"></head><body><article class="document"><a href="../materials.html">← Каталог проекта</a>${note}${body}<hr><p><a href="${escape(item.markdown)}" download>Исходный отчёт Markdown</a></p></article></body></html>`;
 fs.writeFileSync(path.join(reports,item.html),page);
}
const groupNames={horizon_audit:'Yitian A1 · исходный аудит',hf3_continuation:'Yitian A2 · продолжение',hf3_a3:'Yitian A3 · ранняя динамика',hf3_a4:'Yitian A4 · смешение компонент',hf3_r7b:'HF3-R7B · финальное кольцевание QC0',hf4_a1:'HF4-A1 · скалярные карты',hf4_a2:'HF4-A2 · динамический атлас',hf4_a3:'HF4-A3 · инвентаризация геометрии',hf4_b2:'HF4-B2 · поля на поверхностях',hf4_b3:'HF4-B3 · физическая мера',hf4_b4:'HF4-B4 · все общие кадры',hf4_b5:'HF4-B5 · угловое разложение',hf1_visual_v2:'HF1-V2 · настоящая поверхность'};
let groups='';
for(const [stage,name] of Object.entries(groupNames)){
 const rs=catalog.filter(x=>x.stage===stage);let files=[];
 const walk=dir=>{if(!fs.existsSync(dir))return;for(const e of fs.readdirSync(dir,{withFileTypes:true})){const p=path.join(dir,e.name);if(e.isDirectory())walk(p);else files.push(p);}};
 walk(path.join(dist,'downloads',stage));
 groups+=`<section><h2>${escape(name)}</h2><ul>${rs.map(r=>`<li><a href="reports/${r.html}">${escape(r.title)}</a></li>`).join('')}</ul><details><summary>Таблицы, численные результаты и программы · ${files.length} файлов</summary><ul>${files.map(p=>`<li><a href="${path.relative(dist,p).split(path.sep).join('/')}" download>${escape(path.relative(path.join(dist,'downloads',stage),p))}</a> <span class="status-note">${(fs.statSync(p).size/1024).toFixed(0)} КБ</span></li>`).join('')}</ul></details></section>`;
}
const material=`<!doctype html><html lang="ru"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Материалы исследования — Пространства</title><link rel="stylesheet" href="assets/style.css"></head><body><article class="document"><a href="index.html">← К реконструкциям</a><p class="eyebrow" style="margin-top:25px">ВЫПУСК 1.0 · HF1–HF4-B5</p><h1>Материалы исследования</h1><p>Готовые результаты собраны по этапам. Для понимания выводов начните с итогового отчёта; программы нужны для независимого воспроизведения.</p><div class="actions"><a class="button" href="reports/PROJECT_REPORT.html">Итоговый научный отчёт</a><a href="reports/REPRODUCIBILITY.html">Как воспроизвести</a></div><h2>Реконструкции</h2><ul><li><a href="viewers/hf1.html">HF1: выделенная область на 3D-поверхности</a></li><li><a href="viewers/qc0.html">QC0: отдельные и общий горизонты</a></li><li><a href="viewers/yitian.html">Yitian: динамический атлас мультипольных карт</a></li><li><a href="viewers/b5.html">B5: исходное поле, аппроксимация и остаток</a></li><li><a href="assets/hf1-evolution.mp4" download>Видео HF1</a> · <a href="assets/yitian-dynamics.mp4" download>Видео Yitian</a></li></ul>${groups}<section id="sources"><h2>Источники и авторство данных</h2><p>SXS:BBH:0305, собственные QC0-расчёты и данные Yitian для SXS:BBH:0389 имеют разное происхождение. Их параметры и определения не смешиваются.</p><a href="reports/SOURCES.html">Источники, соглашения и контрольные суммы</a><p class="status-note">Исходный авторский архив из переписки и большие сырые сетки не включены. Для B5 сохранён необходимый вход B4. Полнота каждого пути воспроизведения явно указана в документации.</p></section></article></body></html>`;
fs.writeFileSync(path.join(dist,'materials.html'),material);
fs.writeFileSync(path.join(root,'docs/report-catalog.json'),JSON.stringify(catalog,null,2)+'\n');
console.log(JSON.stringify({reports:catalog.length,materials:'dist/materials.html'}));
