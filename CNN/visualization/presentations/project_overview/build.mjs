import fs from 'node:fs/promises';
import path from 'node:path';
import { Presentation, PresentationFile } from '@oai/artifact-tool';
import { resolvePresentationFont, applyPresentationChartFont, finalizePresentation } from '/Users/fabioali/.codex/plugins/cache/openai-primary-runtime/presentations/26.905.11957/skills/presentations/container_tools/artifact_tool_utils.mjs';
const W='/Users/fabioali/SND@LHC/nue_analysis/CNN';
const B=path.join(W,'visualization/presentations/project_overview'), O=path.join(W,'output/presentazione_progetto');
const S='/Users/fabioali/.codex/plugins/cache/openai-primary-runtime/presentations/26.905.11957/skills/presentations';
const D=JSON.parse(await fs.readFile(path.join(B,'data.json'),'utf8'));
const family=resolvePresentationFont({fontFamily:'Arial'}); console.log('font',family);
const p=Presentation.create({slideSize:{width:1280,height:720}});
const C={ink:'#102A43',blue:'#126782',muted:'#526777',white:'#FFFFFF',pale:'#EDF3F6',orange:'#C56730'};
function text(s,t,x,y,w,h,size=28,color=C.ink,bold=false){const a=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});a.text=t;a.text.style={typeface:family,fontSize:size,color,bold,autoFit:'none'};return a;}
function slide(title,notes){const s=p.slides.add();s.background.fill=C.white;text(s,title,64,42,1152,100,43,C.ink,true);text(s,String(p.slides.items.length).padStart(2,'0'),1160,674,55,26,17,C.muted);s.speakerNotes.textFrame.setText(notes);return s;}
function foot(s,t){text(s,t,64,624,1100,50,19,C.muted);}
function para(s,head,body,y,x=64,w=1090){text(s,head,x,y,w,43,28,C.blue,true);text(s,body,x,y+46,w,78,25);}
async function img(s,rel,x,y,w,h){s.images.add({blob:new Uint8Array(await fs.readFile(path.join(W,rel))),contentType:'image/png',alt:rel,fit:'contain',position:{left:x,top:y,width:w,height:h}});}
function table(s,values,widths,y=190,h=330){const t=s.tables.add({rows:values.length,columns:values[0].length,left:64,top:y,width:1152,height:h,columnWidths:widths,values}); t.borders.assign({fill:'#FFFFFF',width:1}); for(let r=0;r<values.length;r++){for(let c=0;c<values[r].length;c++){const a=t.getCell(r,c);a.fill=r===0?C.ink:(r%2?C.pale:C.white);a.text.style={typeface:family,fontSize:r===0?23:25,color:r===0?C.white:C.ink,bold:r===0};}}return t;}
const src=n=>`Fonti: runs/cnn21d_signal_p24_l5/fitlt5_${n}/full/summary.json e scan_t090/score_size_operating_points/summary.json.`;
const fmt=(v,n=1)=>v.toFixed(n).replace('.',',');
// 1
{
const s=p.slides.add();s.background.fill=C.ink;text(s,'SND@LHC',72,70,1000,45,27,'#8DCCDF',true);text(s,'Ricerca di sciami\ncon una CNN 2+1D',72,209,1120,210,66,C.white,true);text(s,'Progetto crop_bw\nSintesi del metodo e dei risultati',76,466,1080,90,30,'#D4E5EF');text(s,'15 settembre 2026',76,637,800,32,22,'#B7CEDC');s.speakerNotes.textFrame.setText('Sintesi del progetto presente nella cartella crop_bw, analisi nue SND@LHC. Stato dei file locali al 15 settembre 2026. Riferimenti principali: README.md, CNN_DATALOADER.md, cnn21d/model.py, configurazioni fitlt5 del 10 settembre e relativi report full. Le prestazioni mostrate sono valutazioni interne sui campioni disponibili.');}
// 2
{
const s=slide('Obiettivo e flusso di analisi','Fonti: README.md, ROOT_TO_HDF5.md, CNN_DATALOADER.md, scanning/scan_cnn21d_volumes.py, scanning/scan_cnn21d_loop2.sh. Il task riconosce strutture di segnale e stima la direzione, poi raggruppa finestre adiacenti in candidati da ispezionare.');
text(s,'Individuare candidati sciame nei volumi e stimarne la direzione',64,151,1100,75,34,C.blue,true);
para(s,'01   Mappe ROOT','Conteggi grezzi per piano e stima del fondo.',260);
para(s,'02   Dataset e CNN','Crop HDF5, classificazione del segnale e regressione delle slope.',374);
para(s,'03   Scansione dei volumi','Finestre sovrapposte, raggruppamento dei candidati e proiezioni XZ/YZ.',488);
}
// 3
{
const s=slide('Dataset e controlli di qualità','Fonti: README.md, data_preparation/root_to_hdf5.py, cnn_dataset_signal_p24_l5_fitlt5_hard.selection_summary.json, configurazioni fitlt5 e summary.json. Dataset fisico: 15872 esempi di cui 8656 hard e 7216 signal (differenza). I subset effettivamente usati dal training sono 10186/2377/1996 e non esauriscono necessariamente il dataset per via dei filtri e delle frazioni di classe. Non sommare i due conteggi come se avessero la stessa selezione. La vecchia descrizione 20x20 in ROOT_TO_HDF5.md è superata dal produttore corrente 32x32.');
text(s,'57 × 32 × 32',64,175,580,85,64,C.blue,true);text(s,'volume grezzo per esempio',68,261,570,45,28);
text(s,'57 × 20 × 20',680,175,530,85,64,C.blue,true);text(s,'crop all’ingresso della rete',684,261,530,45,28);
text(s,'15.872 esempi nel dataset fitlt5\n8.656 hard negatives e 7.216 signal',64,358,1140,95,32,C.ink,true);
text(s,'Separazione per evento di segnale o cella di fondo.\nEsclusione della cella 58 e controllo del supporto ai bordi.\nI conteggi originali restano disponibili senza smoothing.',64,479,1140,119,26);
foot(s,'Train / validation / test effettivi: 10.186 / 2.377 / 1.996 esempi');
}
// 4
{
const s=slide('CNN 2+1D a due teste','Fonti: cnn21d/model.py, cnn21d/losses.py, training/configs/cnn21d_signal_p24_l5_fitlt5_mu_residual.yaml, full/summary.json. Ogni blocco applica Conv3D(1,3,3), GroupNorm, SiLU, Conv3D(3,1,1), GroupNorm, SiLU e MaxPool. I canali sono 16,32,64,128. Le feature globali concatenano media e massimo. Le teste hanno layer nascosto di dimensione 64. Parametri: 196243.');
text(s,'Convoluzioni separate nello spazio e lungo z',64,155,1120,50,32,C.blue,true);
text(s,'1 × 3 × 3',65,242,500,80,57,C.ink,true);text(s,'Struttura del singolo piano XY',69,326,510,80,28);
text(s,'3 × 1 × 1',686,242,500,80,57,C.ink,true);text(s,'Evoluzione tra piani consecutivi',690,326,510,80,28);
para(s,'Due uscite complementari','Score di segnale e componenti della direzione (sx, sy).',446);
foot(s,'4 blocchi, canali 16 / 32 / 64 / 128, 196.243 parametri');
}
// 5
{
const s=slide('Training e trasformazioni coerenti','Fonti: CNN_DATALOADER.md, training/cnn_dataset.py, training/configs/cnn21d_signal_p24_l5_fitlt5_mu_residual.yaml. Figura QA preesistente, illustrativa delle trasformazioni, non un nuovo risultato del run fitlt5. Jitter massimo 5 pixel, rotazioni di 90 gradi e riflessioni coerenti con slope. Layer-drop train: probabilità 12% un piano, 3% due consecutivi. Hard selezionati con fit coerente theta<5 mrad. Nei run fitlt5 la BCE usa signal sopra 8 mrad e hard negativi, esclude signal a piccolo angolo dalla BCE ma li include nella regressione. Checkpoint classificazione sulla AUPRC e regressione sul MAE bilanciato per bin angolari.');
await img(s,'visualization/presentations/project_overview/inputs/qa_02_signal_hdf5_010602.png',55,160,758,403);
text(s,'Jitter ±5 pixel\nRotazioni e riflessioni\nSlope trasformate insieme',858,190,352,155,26,C.ink,true);
text(s,'Signal >8 mrad per la classificazione. Tutti i signal per la regressione.',858,374,352,132,25);
foot(s,'Checkpoint separati per classificazione e regressione. Augmentation solo nel training.');
}
// 6
{
const s=slide('Prestazioni sui crop di test',src('mu_residual')+' '+src('mpv_z_residual')+' '+src('mpv_z_sigma')+' AUPRC, precision e recall dal checkpoint di classificazione a soglia 0.5. MAE theta dal checkpoint di regressione su 879 signal. La classificazione è valutata su 613 positivi e 1117 hard: 266 signal a basso angolo del test sono esclusi dalla BCE. Il MAE è errore assoluto medio, distinto dalla deviazione standard del residuo.');
text(s,'Confronto delle tre normalizzazioni fitlt5',64,150,1130,50,29,C.blue);
const names=['mu_residual','mpv_z_residual','mpv_z_sigma'];
table(s,[['Ingresso CNN','AUPRC','Precisione','Recall','MAE θ'],...names.map((n,i)=>[[ 'H − μ','H − MPV(z)','[H − MPV(z)] / σ(z)'][i],fmt(D[n].test.auprc,4),fmt(D[n].test.precision*100)+'%',fmt(D[n].test.recall*100)+'%',fmt(D[n].reg.theta_mae_mrad,2)+' mrad'])],[395,175,190,180,212],233,295);
foot(s,'Score ≥0,50. Classificazione: 613 signal e 1.117 hard. Regressione: 879 signal.');
}
// 7
{
const s=slide('Scansione dei volumi completi','Fonti: scanning/scan_cnn21d_volumes.py, scanning/scan_cnn21d_loop2.sh, scanning/produce_data_scan_shards.py, candidates/export_background_scan_candidates.py. Configurazione script recente: fitlt5_mu_residual, crop20, stride10, threshold0.90, regressione minima score0.90, due checkpoint separati. I cluster di finestre sono candidati dello scan. Lo score CNN non rappresenta una probabilità fisica calibrata di evento neutrino.');
para(s,'Finestre 20 × 20, passo 10 pixel','La scansione conserva tutti i 57 piani e valuta finestre sovrapposte.',160);
para(s,'Score e direzione','Il checkpoint di classificazione seleziona le finestre. Quello di regressione stima le slope.',295);
para(s,'Candidati da ispezionare','Il raggruppamento delle finestre produce candidati con coordinate, score e proiezioni XZ/YZ.',430);
foot(s,'Script recente sui dati: H − μ, score ≥0,90. Il candidato richiede validazione fisica.');
}
// 8
{
const s=slide('Efficienza dello scan per angolo',src('mu_residual')+' Fonte aggiuntiva: scan_t090/signal_scan_summary.json, operating_points threshold0.90. Matching con verità MC, recomputed_propagation_limit_40. I bin sono 5-10, 10-20, 20-50, 50-100 mrad con N=223,245,253,37, trovati=114,240,252,36. Nessun evento >=100mrad. Efficienza aggregata theta>10mrad=528/535. Nessun taglio di size.');
text(s,'H − μ, score ≥0,90, associazione alla verità MC',64,148,1150,46,27,C.blue);
const chart=s.charts.add('bar',{position:{left:60,top:219,width:830,height:359},categories:['5–10','10–20','20–50','50–100'],series:[{name:'Efficienza',values:[51.1,98.0,99.6,97.3],fill:C.blue}],barOptions:{direction:'column',grouping:'clustered',gapWidth:95},hasLegend:false,xAxis:{title:'θ vero [mrad]',textStyle:{fontSize:21,fill:C.ink}},yAxis:{min:0,max:110,majorUnit:25,numberFormatCode:'0.0"%"',textStyle:{fontSize:20},majorGridlines:{fill:'#DCE6EC',width:1}},dataLabels:{showValue:true,position:'outEnd',textStyle:{fontSize:23,bold:true,fill:C.ink}}}); applyPresentationChartFont(chart,{fontFamily:family});
text(s,'98,7%',940,256,280,80,58,C.blue,true);text(s,'528 / 535 eventi\ncon θ >10 mrad',944,350,270,105,28);text(s,'Il limite principale è\nla regione 5–10 mrad.',944,493,270,90,25);
foot(s,'Eventi per bin: 223 / 245 / 253 / 37. Nessun taglio sulla dimensione del candidato.');
}
// 9
{
const s=slide('Soglia dello score e candidati di fondo',src('mu_residual')+' Dati: score_size_operating_points/summary.json. Efficienza truth-matched su 535 signal theta>10 mrad. Il conteggio fondo è il numero di cluster su 323 mappe e 116603 finestre, non una probabilità per finestra o un tasso fisico normalizzato. La riduzione 184 a72 è60.87%.');
text(s,'H − μ: scansione su 323 mappe di fondo',64,150,1150,48,28,C.blue);
table(s,[['Score minimo','Signal associati','Efficienza','Cluster nel fondo'],...D.mu_residual.op.rows.map(r=>[fmt(r.threshold,2),r.signal_truth_matched+' / 535',fmt(r.signal_truth_efficiency*100,2)+'%',String(r.background_clusters)])],[250,315,247,340],221,316);
foot(s,'Da 0,90 a 0,99: 61% di cluster di fondo in meno e 0,75 punti di efficienza in meno.');
}
// 10
{
const s=slide('Normalizzazione e taglio sulla dimensione',src('mu_residual')+' '+src('mpv_z_residual')+' '+src('mpv_z_sigma')+' Stesso campione di535 signal theta>10 e323 mappe background. La size è largest_component_voxels>=8. I risultati suggeriscono un compromesso interno al benchmark, non una scelta definitiva per i dati di tutti i brick.');
text(s,'Stesso benchmark, score ≥0,90',64,149,1150,50,28,C.blue);
table(s,[['Ingresso CNN','Efficienza','Cluster nel fondo'],...['mu_residual','mpv_z_residual','mpv_z_sigma'].map((n,i)=>[['H − μ','H − MPV(z)','[H − MPV(z)] / σ(z)'][i],fmt(D[n].op.rows[0].signal_truth_efficiency*100,2)+'%',String(D[n].op.rows[0].background_clusters)])],[535,295,322],219,259);
text(s,'Il taglio sulla dimensione ha un costo elevato',64,511,1130,43,31,C.blue,true);
text(s,'Con H − μ e almeno 8 voxel: efficienza 84,9% e 124 cluster di fondo.',64,563,1130,66,27);
}
// 11
{
const s=slide('Ispezione di un candidato nel fondo','Fonte figura: runs/cnn21d_signal_p24_l5/fitlt5_mu_residual/full/scan_t090/background_scan_t090/xz_yz_new_vs_gap5_10/rank_001_event_00271_cell_271_score_0.999890_xz_yz.png. La figura identifica hard_b000021_hmu_new_vs_gap5_10, rank1, event271 cell271. Score esatto0.999890 dal nome file. Theta CNN8.7mrad e fit9.5mrad dal titolo. Questo esempio illustra l’ispezione di un candidato nel campione background e non costituisce conferma di segnale fisico.');
await img(s,'runs/cnn21d_signal_p24_l5/fitlt5_mu_residual/full/scan_t090/background_scan_t090/xz_yz_new_vs_gap5_10/rank_001_event_00271_cell_271_score_0.999890_xz_yz.png',50,164,1170,445);
foot(s,'Cella 271: score 0,99989, θ CNN 8,7 mrad, θ fit 9,5 mrad. Candidato da validare.');
}
// 12
{
const s=slide('Stato del progetto e prossimi passi','Sintesi basata sui file locali. Risultati: pipeline implementata, test crop AUPRC~0.999, scan H-mu theta>10 528/535 a0.90. Prossimi passi proposti dall’analisi: valutare robustezza per brick e regime di fondo, fissare il punto operativo su validation indipendente, audit truth matching e candidati non associati, validazione manuale/fit su dati reali. Non vi sono in queste fonti una misura finale di purezza fisica sui dati né un risultato di osservazione.');
para(s,'Pipeline operativa','Produzione dei dataset, training, scansione e diagnostica dei candidati sono disponibili.',162);
para(s,'Risultato principale','Nel benchmark H − μ lo scan recupera il 98,7% del segnale sopra 10 mrad a score 0,90.',296);
para(s,'Validazione da completare','Robustezza tra brick e livelli di fondo, piccoli angoli e verifica fisica dei candidati.',430);
foot(s,'Prossimo passo proposto: fissare il punto operativo su validazione indipendente e dati reali.');
}
await fs.writeFile(path.join(B,'deck.json'),JSON.stringify(p.toProto()));
await (await PresentationFile.exportPptx(p)).save(path.join(B,'candidate.pptx'));
for(let i=0;i<p.slides.items.length;i++){const b=await p.export({slide:p.slides.items[i],format:'png',scale:1});await fs.writeFile(path.join(B,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await b.arrayBuffer()));}
const result=await finalizePresentation({workspaceDir:W,candidatePath:path.join(B,'candidate.pptx'),finalPath:path.join(O,'Sintesi_progetto_CNN_SND_v2.pptx'),pythonExecutable:'/Users/fabioali/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3',integrityValidatorPath:path.join(S,'container_tools/inspect_presentation_package_integrity.py'),layoutValidatorPath:path.join(S,'container_tools/inspect_presentation_layout_geometry.py'),layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit','--require-native-table-slide','6','--require-native-table-slide','9','--require-native-table-slide','10'],fontPolicy:{basis:'design',families:[family]},explicitTotalSlideCount:12,requiredNativeTableOwnerSlides:[6,9,10],requiredNativeChartOwnerSlides:[8],materializeLiteralChartWorkbooks:true,verifyArtifactToolImport:true,receiptPath:path.join(B,'validation.json')});
console.log(result);
