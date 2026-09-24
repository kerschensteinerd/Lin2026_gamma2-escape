#target illustrator
(function(){
var root='/Users/kerschensteinerd/Documents/g2_loom_cg/';
var log=new File(root+'tmp/fig4_v9/native_export_log.txt');
var prior=app.userInteractionLevel;
function report(s){log.open('a');log.writeln(s);log.close();}
function quote(s){return '"'+String(s).replace(/\\/g,'\\\\').replace(/"/g,'\\"').replace(/\r/g,'\\r').replace(/\n/g,'\\n')+'"';}
function json(v){if(v===null)return 'null';if(typeof v==='string')return quote(v);if(typeof v!=='object')return String(v);var a=[],k;if(v instanceof Array){for(k=0;k<v.length;k++)a.push(json(v[k]));return '['+a.join(',')+']';}for(k in v)if(v.hasOwnProperty(k))a.push(quote(k)+':'+json(v[k]));return '{'+a.join(',')+'}';}
function inspect(d){var r=d.artboards[0].artboardRect;return {name:d.name,artboards:d.artboards.length,page_points:[r[2]-r[0],r[1]-r[3]],editable_text_frames:d.textFrames.length,vector_paths:d.pathItems.length,external_links:d.placedItems.length,embedded_rasters:d.rasterItems.length};}
var jobs=[{source:'figures/main/PDF/Fig4_v9.pdf',target:'figures/main/AI/Fig4_v9.ai',qa:'figures/main/assets/Fig4_v9/native_illustrator_qa.json',width:493.228,height:257.509,minpaths:90000,mintext:150},{source:'output/pdf/Figure7_v8.pdf',target:'figures/main/AI/Fig7_v8.ai',qa:'figures/main/assets/Fig7_v8/native_illustrator_qa.json',width:493.228,height:535,minpaths:1000,mintext:300}];
try{
app.userInteractionLevel=UserInteractionLevel.DONTDISPLAYALERTS;
for(var k=0;k<jobs.length;k++){
var job=jobs[k],src=new File(root+job.source),dst=new File(root+job.target);
if(!src.exists)throw Error('Missing source '+src.fsName);if(dst.exists)throw Error('Destination exists '+dst.fsName);
report('Opening '+src.fsName);var d=app.open(src);report('Opened '+d.name);var before=inspect(d);
if(before.artboards!==1||Math.abs(before.page_points[0]-job.width)>.02||Math.abs(before.page_points[1]-job.height)>.02)throw Error('Wrong artboard');
if(before.editable_text_frames<job.mintext||before.vector_paths<job.minpaths||before.external_links!==0)throw Error('Unexpected editability '+json(before));
var opt=new IllustratorSaveOptions();opt.pdfCompatible=true;opt.compressed=true;opt.embedICCProfile=true;opt.embedLinkedFiles=true;opt.fontSubsetThreshold=100;
if(job.source.indexOf('Fig4_v9.pdf')>=0){var allText='';for(var ti=0;ti<d.textFrames.length;ti++)allText+=d.textFrames[ti].contents+' ';if(allText.indexOf('25% retained')<0||allText.indexOf('GABA removed')>=0)throw Error('Wrong model labels imported');}
report('Saving native '+dst.fsName+' '+json(before));d.saveAs(dst,opt);if(!d.saved||!dst.exists)throw Error('Save incomplete');d.close(SaveOptions.DONOTSAVECHANGES);
report('Reopening '+dst.fsName);var reopened=app.open(dst);var after=inspect(reopened);
if(after.editable_text_frames!==before.editable_text_frames||after.vector_paths!==before.vector_paths||after.external_links!==0)throw Error('Editability changed after save');
var qa=new File(root+job.qa);qa.open('w');qa.write(json({source:src.fsName,output:dst.fsName,illustrator_version:app.version,native_save_reopened:true,pdf_compatible:true,before:before,after:after}));qa.close();reopened.close(SaveOptions.DONOTSAVECHANGES);report('SUCCESS '+dst.fsName);
}
return 'SUCCESS: both native artifacts saved and reopened';
}catch(e){report('ERROR '+e+' line '+e.line);throw e;}finally{app.userInteractionLevel=prior;}
})();
