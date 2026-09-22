// Controller and geometry-data check without a browser or WebGL renderer.
const fs=require('fs'),vm=require('vm'),assert=require('assert'),path=require('path');
const root=__dirname, data=JSON.parse(fs.readFileSync(path.join(root,'results/BH_HF1_V2_viewer_data.json'),'utf8'));
const template=fs.readFileSync(path.join(root,'viewer_template.html'),'utf8');
const scripts=[...template.matchAll(/<script>([\s\S]*?)<\/script>/g)];const code=scripts.at(-1)[1];
const elements=new Map(),timers=new Map();let now=0,nextTimer=1,lastTraces,lastLayout;
function element(id){if(!elements.has(id))elements.set(id,{id,style:{},handlers:{},events:{},value:'0',textContent:'',
 setAttribute(k,v){this[k]=v},addEventListener(k,v){this.handlers[k]=v},on(k,v){this.events[k]=v}});return elements.get(id);}
element('sourceData').textContent=JSON.stringify(data);
const jump=[0,2,5,10,20,50].map(t=>({...element('jump'+t),dataset:{time:String(t)}}));
const context={document:{getElementById:element,querySelectorAll:()=>jump},window:{},structuredClone,console,
 performance:{now:()=>now},requestAnimationFrame:f=>{const id=nextTimer++;timers.set(id,f);return id},cancelAnimationFrame:id=>timers.delete(id),
 Plotly:{react:async(el,traces,layout)=>{lastTraces=traces;lastLayout=layout},relayout:async(el,change)=>{if(el.events.plotly_relayout)el.events.plotly_relayout(change)}}};
vm.createContext(context);new vm.Script(code).runInContext(context);
const drain=()=>new Promise(resolve=>setImmediate(resolve));
async function run(){
 await drain();const api=context.window.HF1V2;assert(api);
 assert.strictEqual(api.getState().current,api.nearest(2));
 for(let k=0;k<3;k++){
  element('field').handlers.change({target:{value:String(k)}});await drain();
  const st=api.getState(),f=data.frames[st.current];assert.strictEqual(st.field,k);
  assert.strictEqual(lastTraces[2].x.length,f.masks[k].length);
  assert.strictEqual(lastTraces[0].facecolor.length,lastTraces[0].i.length);
  assert(element('status').textContent.includes(String(f.masks[k].length)));
 }
 const st=api.getState();const eye={...st.camera.eye};element('opposite').handlers.click();
 assert.strictEqual(api.getState().current,st.current);assert.strictEqual(api.getState().camera.eye.x,-eye.x);
 element('time').handlers.input({target:{value:'0'}});await drain();
 element('play').handlers.click();assert(api.getState().playing&&!api.getState().orbiting);
 now=1000;let id=Math.max(...timers.keys()),cb=timers.get(id);timers.delete(id);cb(now);await drain();
 assert.strictEqual(api.getState().current,api.previous(.30));
 element('orbit').handlers.click();const fixed=api.getState().current;assert(!api.getState().playing&&api.getState().orbiting);
 now=1600;id=Math.max(...timers.keys());cb=timers.get(id);timers.delete(id);cb(now);await drain();
 assert.strictEqual(api.getState().current,fixed);
 jump.find(b=>b.dataset.time==='20').handlers.click();await drain();assert.strictEqual(element('late').style.display,'block');
 assert(!api.getState().playing&&!api.getState().orbiting);
 const geometry=[];
 for(const t of [0,2,5,10,20,50])for(let k=0;k<3;k++){
  const j=api.nearest(t),s=api.split(data.frames[j],k);geometry.push({row:data.frames[j].row,k,...s});
 }
 fs.writeFileSync(path.join(root,'results/viewer_geometry_check.json'),JSON.stringify(geometry));
 const result={javascript_syntax:true,field_switch_preserves_node_counts:true,opposite_view_preserves_time:true,
  play_uses_simulation_time:true,orbit_stops_time_playback:true,late_frames_are_explicitly_separate:true,
  browser_webgl_render_tested:false,viewer_snapshot_count:data.frames.length};
 fs.writeFileSync(path.join(root,'results/BH_HF1_V2_viewer_controller_check.json'),JSON.stringify(result,null,2));
 console.log(JSON.stringify(result));
}
run().catch(e=>{console.error(e);process.exit(1)});
