/* Functional and numerical audit using Node and a real canvas engine.
   This is not a full browser/layout test. No network requests are made. */
const fs=require('fs'),vm=require('vm');
const {createCanvas}=require(process.env.CODEX_PRIMARY_RUNTIME_NODE_MODULES+'/@napi-rs/canvas');
const filename=process.argv[2],out=process.argv[3];
const html=fs.readFileSync(filename,'utf8');
const data=html.match(/<script id="atlas-data" type="application\/json">([\s\S]*?)<\/script>/)[1];
const source=html.match(/<script>\s*([\s\S]*?)<\/script>/)[1];
const elements={};
for(const m of html.matchAll(/id="([^"]+)"/g))elements[m[1]]={textContent:'',value:'',replaceChildren(...items){this.children=items}};
elements['atlas-data'].textContent=data;
for(const [id,value] of Object.entries({field:'0',frame:'original',scale:'linear',speed:'1',time:'0'}))elements[id].value=value;
for(const [id,w,h] of [['map4',582,294],['map6',582,294],['map8',582,294],['curve',630,180],['legend',768,18]]){
 elements[id]=createCanvas(w,h);elements[id].getBoundingClientRect=()=>({left:0,top:0,width:w,height:h});
}
const document={getElementById:id=>elements[id],createElement:tag=>tag==='canvas'?createCanvas(1,1):{textContent:''}};
let lastCallback=null;
const context={document,window:{},Uint8Array,Float64Array,Math,JSON,atob:s=>Buffer.from(s,'base64').toString('binary'),requestAnimationFrame:f=>{lastCallback=f;return 1},cancelAnimationFrame:()=>{}};
vm.runInNewContext(source,context,{timeout:30000});
const api=context.window.HF4A2;
const numeric=api.validate();
if(numeric.maxAbsError>1e-10)throw Error('Numerical atlas mismatch');
let modes=0;
for(let field=0;field<4;field++)for(let rotate of ['original','rotated'])for(let scale of ['linear','asinh']){
 elements.field.value=String(field);elements.frame.value=rotate;elements.scale.value=scale;
 api.setFrame(field===3?300:160);
 const settings=api.getSettings();if(settings.field!==field||settings.rotate!==(rotate==='rotated')||settings.scale!==scale)throw Error('Control mismatch');
 const values=api.synthesize(field,settings.frame,8,settings.rotate);
 if(!Array.from(values).every(Number.isFinite))throw Error('Non-finite field');
 modes++;
}
api.setFrame(0);elements.next.onclick();if(api.getSettings().frame!==1)throw Error('Next control');
elements.back.onclick();if(api.getSettings().frame!==0)throw Error('Previous control');
elements.play.onclick();lastCallback(1000);if(api.getSettings().frame!==1)throw Error('Playback control');elements.play.onclick();
elements.field.value='3';elements.frame.value='original';elements.scale.value='linear';api.setFrame(0);
for(const cut of [4,6,8])fs.writeFileSync(out.replace(/\.json$/,`_map${cut}.png`),elements['map'+cut].toBuffer('image/png'));
const result={method:'Node VM with real canvas and minimal DOM; full browser layout unavailable',numeric,control_modes_checked:modes,previous_next_play_checked:true,frames:api.frameCount,all_pass:true};
fs.writeFileSync(out,JSON.stringify(result,null,2)+'\n');
console.log(JSON.stringify({maxAbsError:numeric.maxAbsError,modes,frames:api.frameCount,all_pass:true}));
