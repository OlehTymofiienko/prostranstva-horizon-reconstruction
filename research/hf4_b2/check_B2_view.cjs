// Structural and interaction checks; this is not a WebGL/browser render test.
const fs = require('fs'), vm = require('vm'), assert = require('assert'), path = require('path');
const out = path.join(__dirname, 'results');
const D = JSON.parse(fs.readFileSync(path.join(out, 'BH_HF4_B2_view_data.json'), 'utf8'));
const elements = new Map(), intervals = new Map(); let counter=0, calls=0;
function element(id) {
  if (!elements.has(id)) elements.set(id, {id,value:'',checked:false,textContent:'',attrs:{},setAttribute(k,v){this.attrs[k]=v;}});
  return elements.get(id);
}
element('data').textContent=JSON.stringify(D); element('field').value='repsi2';element('mask').checked=true;
const context = {document:{getElementById:element},window:{},console,
  setInterval:fn=>{intervals.set(++counter,fn);return counter;},clearInterval:id=>intervals.delete(id),
  Plotly:{react:(node,traces,layout)=>{calls++;node.lastTraces=traces;return Promise.resolve();},
          relayout:()=>Promise.resolve()}};
vm.createContext(context);
vm.runInContext(fs.readFileSync(path.join(out,'BH_HF4_B2_controller.js'),'utf8'),context);
(async()=>{
 const api=context.window.HF4B2; await api.ready();
 let tested=0;
 for(const f of D.fields) for(let index=0;index<D.frames.length;index++)
 for(const mask of [false,true]) for(const inners of [false,true]) {
  element('field').value=f;element('mask').checked=mask;element('inners').checked=inners;
  await api.setFrame(index);const frame=D.frames[index], traces=element('plot').lastTraces;
  const common=frame.objects.some(o=>o.horizon===3);
  const bases=traces.filter(t=>t.name&&/^AH[123]$/.test(t.name));
  assert.equal(bases.length,common?(inners?3:1):2);
  assert.equal(bases.some(t=>t.name==='AH3'),common);
  for(const b of bases) {
   assert.equal(b.x.length,2520);assert.equal(b.i.length,4896);
   for(const a of [b.x,b.y,b.z]) assert(a.every(Number.isFinite));
   if(b.name==='AH3') {assert.deepEqual(Array.from(b.intensity),frame.fields[f]);assert.equal(b.cmin,D.ranges[f][0]);assert.equal(b.cmax,D.ranges[f][1]);}
  }
  for(const t of traces.filter(t=>t.type==='mesh3d')) {
   assert(t.i.length===t.j.length&&t.i.length===t.k.length);
   for(const a of [t.i,t.j,t.k])assert(a.every(k=>Number.isInteger(k)&&k>=0&&k<t.x.length));
  }
  const marks=traces.filter(t=>t.type==='scatter3d');assert.equal(marks.length,common&&mask?1:0);
  if(marks.length){const obj=frame.objects.find(o=>o.horizon===3);assert.equal(marks[0].x.length,frame.masks[f].length);
   for(let k=0;k<marks[0].x.length;k++){const id=frame.masks[f][k];assert.equal(marks[0].x[k],obj.xyz[3*id]);assert.equal(marks[0].y[k],obj.xyz[3*id+1]);assert.equal(marks[0].z[k],obj.xyz[3*id+2]);}}
  tested++;
 }
 await api.setFrame(9);element('play').onclick();assert(api.getState().playing);
 for(const fn of [...intervals.values()])fn();
 assert.equal(api.getState().index,9);assert(!api.getState().playing);assert(element('status').textContent.includes('1,3125'));
 element('next').onclick();await api.ready();assert.equal(api.getState().index,10);assert(element('status').textContent.includes('1,3125'));
 await api.setFrame(4);element('play').onclick();for(const fn of [...intervals.values()])fn();await api.ready();assert.equal(api.getState().index,5);element('play').onclick();
 element('orbit').onclick();assert(api.getState().orbiting);for(const fn of [...intervals.values()])fn();element('orbit').onclick();assert(!api.getState().orbiting);
 assert.equal(intervals.size,0);assert.equal(element('error').textContent,'');
 const result={controller_states_tested:tested,render_calls:calls,source_masks_match:true,gap_playback_stops:true,camera_control_dispatch:true,browser_WebGL_render_tested:false};
 fs.writeFileSync(path.join(out,'BH_HF4_B2_viewer_checks.json'),JSON.stringify(result,null,2)+'\n');console.log(JSON.stringify(result));
})().catch(e=>{console.error(e);process.exit(1);});
