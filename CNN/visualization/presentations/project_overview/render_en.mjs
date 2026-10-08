import fs from 'node:fs/promises';
import {FileBlob, PresentationFile} from '@oai/artifact-tool';
const root='/Users/fabioali/SND@LHC/nue_analysis/CNN';
await fs.mkdir(root+'/visualization/presentations/project_overview/en', {recursive:true});
const p=await PresentationFile.importPptx(await FileBlob.load(root+'/output/presentazione_progetto/CNN_SND_Project_Overview_EN.pptx'));
for(let i=0;i<p.slides.items.length;i++){
 const b=await p.export({slide:p.slides.items[i],format:'png',scale:1});
 await fs.writeFile(root+'/visualization/presentations/project_overview/en/final-'+String(i+1).padStart(2,'0')+'.png',new Uint8Array(await b.arrayBuffer()));
}
