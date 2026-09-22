// Checks real serialized data and controller transitions. No WebGL claim.
const fs=require('fs'),vm=require('vm'),assert=require('assert'),path=require('path');
const out=path.join(__dirname,'results');const D=JSON.parse(fs.readFileSync(path.join(out,'BH_HF4_B4_view_data.json'),'utf8'));
const elements=new Map(),timers=new Map();let counter=0,calls=0;
function el(id){if(!elements.has(id))elements.set(id,{id,value:'',checked:false,textContent:'',setAttribute(k,v){this[k]=v;}});return elements.get(id);}
el('data').textContent=JSON.stringify(D);el('field').value='repsi2';el('mask').checked=true;
const context={document:{getElementById:el},window:{},console,setInterval:fn=>{timers.set(++counter,fn);return counter;},clearInterval:id=>timers.delete(id),
Plotly:{react:(node,traces)=>{calls++;node.lastTraces=traces;return Promise.resolve();},relayout:()=>Promise.resolve()}};
vm.createContext(context);vm.runInContext(fs.readFileSync(path.join(out,'BH_HF4_B4_controller.js'),'utf8'),context);
(async()=>{
 const api=context.window.HF4B4;await api.ready();let checked=0;
 // All frames with the default field; both modes on every common frame.
 for(let n=0;n<D.frames.length;n++)for(const inners of (D.frames[n].fields?[false,true]:[false])){
  el('inners').checked=inners;await api.setFrame(n);const F=D.frames[n];const traces=el('plot').lastTraces;
  const expected=F.objects.filter(o=>!F.fields||inners||o.horizon===3);
  const bases=traces.filter(t=>/^AH[123]$/.test(t.name||''));assert.equal(bases.length,expected.length);
  for(const t of bases){const o=expected.find(o=>'AH'+o.horizon===t.name);assert(o);assert.equal(t.x.length,o.xyz.length/3);
   assert.equal(t.i.length,D.meshes[o.mesh].length);for(const a of [t.x,t.y,t.z])assert(a.every(Number.isFinite));}
  for(const t of traces.filter(t=>t.type==='mesh3d'))for(const a of [t.i,t.j,t.k])assert(a.every(k=>Number.isInteger(k)&&k>=0&&k<t.x.length));
  if(F.fields){const t=bases.find(t=>t.name==='AH3');assert.deepEqual(Array.from(t.intensity),F.fields.repsi2);}
  checked++;
 }
 // Every field's marker indices and values on all 55 common frames.
 for(const field of D.fields)for(let n=D.birth_index;n<D.frames.length;n++){
  el('field').value=field;el('inners').checked=false;el('mask').checked=true;await api.setFrame(n);
  const F=D.frames[n],o=F.objects.find(o=>o.horizon===3),tr=el('plot').lastTraces,marks=tr.find(t=>t.type==='scatter3d');
  assert.equal(marks.x.length,F.masks[field].length);
  F.masks[field].forEach((id,k)=>{assert.equal(marks.x[k],o.xyz[id*3]);assert.equal(marks.y[k],o.xyz[id*3+1]);assert.equal(marks.z[k],o.xyz[id*3+2]);});
  const base=tr.find(t=>t.name==='AH3');assert.deepEqual(Array.from(base.intensity),F.fields[field]);
  assert.equal(base.cmin,D.ranges[field][0]);assert.equal(base.cmax,D.ranges[field][1]);checked++;
 }
 el('mask').checked=false;await api.setFrame(D.birth_index);assert(!el('plot').lastTraces.some(t=>t.type==='scatter3d'));
 // A missing pre-common snapshot clears surfaces and pauses playback.
 await api.setFrame(5);el('play').onclick();for(const fn of [...timers.values()])fn();await api.ready();
 assert.equal(api.getState().index,6);assert(!api.getState().playing);assert.equal(el('plot').lastTraces.length,0);
 // Playback traverses the former B3 gap using actual new frames.
 await api.setFrame(9696/32);el('play').onclick();for(const fn of [...timers.values()])fn();await api.ready();
 assert.equal(D.frames[api.getState().index].iteration,9728);assert(api.getState().playing);el('play').onclick();
 // Missing individual horizons never remove the stored common horizon.
 for(const it of [10656,10752]){el('inners').checked=true;await api.setFrame(it/32);assert.equal(el('plot').lastTraces.filter(t=>/^AH[123]$/.test(t.name||'')).length,1);}
 el('begin').onclick();await api.ready();assert.equal(api.getState().index,0);
 el('birth').onclick();await api.ready();assert.equal(api.getState().index,D.birth_index);
 el('orbit').onclick();assert(api.getState().orbiting);for(const fn of [...timers.values()])fn();el('orbit').onclick();assert(!api.getState().orbiting);
 assert.equal(timers.size,0);assert.equal(el('error').textContent,'');
 const result={controller_states_checked:checked,render_calls:calls,all_common_fields_and_marker_indices_verified:true,
 missing_frames_clear_and_pause:true,former_B3_gap_traversed_with_real_frames:true,missing_individuals_do_not_hide_common:true,
 camera_control_dispatch:true,browser_WebGL_render_tested:false};
 fs.writeFileSync(path.join(out,'BH_HF4_B4_controller_checks.json'),JSON.stringify(result,null,2)+'\n');console.log(JSON.stringify(result));
})().catch(e=>{console.error(e);process.exit(1);});
