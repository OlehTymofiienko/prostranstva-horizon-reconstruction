// Run the shipped controller against its shipped data; no browser-render claim.
const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert');
const root=path.resolve(__dirname,'..'),dir=path.join(root,'dist/viewers'),elements=new Map(),timers=new Map();let next=0,draws=0;
function ctx(){return {createImageData:(w,h)=>({data:new Uint8ClampedArray(w*h*4)}),putImageData(){},drawImage(){draws++},beginPath(){},moveTo(){},lineTo(){},stroke(){}};}
function el(id){if(!elements.has(id))elements.set(id,{value:'0',checked:false,hidden:false,textContent:'',handlers:{},getContext:ctx,setAttribute(k,v){this[k]=v},addEventListener(k,v){this.handlers[k]=v}});return elements.get(id);}
el('degree').value='3';el('mask').checked=true;
const context={console,Uint8Array,Uint8ClampedArray,Float32Array,Map,atob:s=>Buffer.from(s,'base64').toString('binary'),
 document:{hidden:false,getElementById:el,createElement:()=>({getContext:ctx}),addEventListener(){}},
 setInterval:f=>{timers.set(++next,f);return next;},clearInterval:id=>timers.delete(id)};
context.window=context;vm.createContext(context);
for(const name of fs.readdirSync(path.join(dir,'b5-data')).sort())vm.runInContext(fs.readFileSync(path.join(dir,'b5-data',name),'utf8'),context);
vm.runInContext(fs.readFileSync(path.join(dir,'b5-controller.js'),'utf8'),context);
assert.equal(context.B5.frames.length,55);assert.equal(el('loading').hidden,true);let count=0;
for(let t=0;t<55;t++)for(let f=0;f<3;f++)for(let l=0;l<4;l++){
 el('time').value=String(t);el('field').value=String(f);el('degree').value=String(l);el('field').handlers.change();
 const m=context.B5.frames[t].fields[f].metrics[l];
 assert.equal(el('r2').textContent,(100*m.R2).toFixed(4)+'%');assert.equal(el('dice').textContent,m.dice.toFixed(4));
 assert(el('clock').textContent.includes(context.B5.frames[t].t.toFixed(4)));assert.equal(el('error').textContent,'');count++;
}
el('time').value='54';el('play').handlers.click();assert.equal(el('time').value,0);assert.equal(timers.size,1);
for(const fn of timers.values())fn();assert.equal(el('time').value,1);el('time').handlers.input();assert.equal(timers.size,0);
el('mask').checked=false;el('mask').handlers.change();assert.equal(el('error').textContent,'');
const result={states_checked:count,canvas_draw_dispatches:draws,fields_degrees_times_and_metric_readouts:true,playback_restarts_at_first_frame:true,scrubbing_pauses_playback:true,browser_canvas_render_tested:false};
fs.writeFileSync(path.join(root,'docs/b5-viewer-validation.json'),JSON.stringify(result,null,2)+'\n');console.log(JSON.stringify(result));
