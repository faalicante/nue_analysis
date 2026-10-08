import fs from 'node:fs/promises';
import {FileBlob, PresentationFile} from '@oai/artifact-tool';
const root='/Users/fabioali/SND@LHC/nue_analysis/CNN';
const p=await PresentationFile.importPptx(await FileBlob.load(root+'/output/presentazione_progetto/Sintesi_progetto_CNN_SND.pptx'));
for(let i=0;i<p.slides.items.length;i++){
 const b=await p.export({slide:p.slides.items[i],format:'png',scale:1});
 await fs.writeFile(root+'/visualization/presentations/project_overview/final-'+String(i+1).padStart(2,'0')+'.png',new Uint8Array(await b.arrayBuffer()));
}
